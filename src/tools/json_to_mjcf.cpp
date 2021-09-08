#include <CLI/CLI.hpp>
#include <igl/PI.h>
#include <nlohmann/json.hpp>
#include <tinyxml2.h>

#include <logger.hpp>

int json_to_mjcf(
    const std::string json_path,
    int mode = 0) // mode0: bullet, mode1: mujoco
{
    if (mode > 1) {
        spdlog::error("Invalid mode in mjcf file generation!");
        exit(-1);
    }

    std::ifstream input(json_path);
    if (input.good()) {
        nlohmann::json scene = nlohmann::json::parse(input, nullptr, false);
        if (scene.is_discarded()) {
            spdlog::error("Invalid Json file");
            input.close();
            return -1;
        }

        std::string xmlFilePath =
            json_path.substr(0, json_path.find_last_of('.'))
            + ((mode == 0) ? ".xml" : "_mjc.xml");
        FILE* output = fopen(xmlFilePath.c_str(), "w");
        if (!output) {
            spdlog::error("failed to create file");
            input.close();
            return -1;
        }

        // transform format
        tinyxml2::XMLPrinter printer(output);
        printer.OpenElement("mujoco");

        std::set<std::string> meshFilePathSet;
        std::vector<std::string> meshNames;
        printer.OpenElement("asset");
        int rbId = 0;
        for (const auto& rbI : scene["rigid_body_problem"]["rigid_bodies"]) {
            std::string meshFilePath = rbI["mesh"].get<std::string>();
            meshNames.emplace_back(
                meshFilePath.substr(0, meshFilePath.find_last_of('.')));
            if (meshFilePathSet.find(meshFilePath) == meshFilePathSet.end()) {
                meshFilePathSet.insert(meshFilePath);

                printer.OpenElement("mesh");
                if (mode == 0) {
                    printer.PushAttribute("file", meshFilePath.c_str());
                } else {
                    printer.PushAttribute(
                        "file", (meshNames.back() + ".stl").c_str());
                }
                printer.PushAttribute("name", meshNames.back().c_str());
                printer.CloseElement(); // mesh
            }
            if (mode == 1
                && (rbI.find("dimensions") != rbI.end()
                    || rbI.find("scale") != rbI.end())) {
                printer.OpenElement("mesh");
                printer.PushAttribute(
                    "file", (meshNames.back() + ".stl").c_str());
                meshNames.back() += "_" + std::to_string(rbId);
                printer.PushAttribute("name", meshNames.back().c_str());
                if (rbI.find("scale") != rbI.end()) {
                    try {
                        double scale = rbI["scale"].get<double>();
                        std::string sizeStr(std::to_string(scale));
                        sizeStr += " " + sizeStr + " " + sizeStr;
                        printer.PushAttribute("scale", sizeStr.c_str());
                    } catch (...) {
                        std::vector<double> scale =
                            rbI["scale"].get<std::vector<double>>();
                        if (scale.size() != 3) {
                            spdlog::error("scale dimension error!");
                            exit(-1);
                        }
                        for (auto& i : scale) {
                            if (!i) {
                                i = 1;
                            }
                        }
                        std::string sizeStr(std::to_string(scale[0]));
                        sizeStr += " " + std::to_string(scale[1]) + " "
                            + std::to_string(scale[2]);
                        printer.PushAttribute("scale", sizeStr.c_str());
                    }
                } else if (rbI.find("dimensions") != rbI.end()) {
                    try {
                        double scale = rbI["dimensions"].get<double>();
                        std::string sizeStr(std::to_string(scale));
                        sizeStr += " " + sizeStr + " " + sizeStr;
                        printer.PushAttribute("scale", sizeStr.c_str());
                    } catch (...) {
                        std::vector<double> scale =
                            rbI["dimensions"].get<std::vector<double>>();
                        if (scale.size() != 3) {
                            spdlog::error("scale dimension error!");
                            exit(-1);
                        }
                        for (auto& i : scale) {
                            if (!i) {
                                i = 1;
                            }
                        }
                        std::string sizeStr(std::to_string(scale[0]));
                        sizeStr += " " + std::to_string(scale[1]) + " "
                            + std::to_string(scale[2]);
                        printer.PushAttribute("scale", sizeStr.c_str());
                    }
                }
                printer.CloseElement(); // mesh
            }

            ++rbId;
        }
        printer.CloseElement(); // asset

        double fricCoef = 0;
        if (scene["rigid_body_problem"].find("coefficient_friction")
            != scene["rigid_body_problem"].end()) {
            fricCoef = scene["rigid_body_problem"]["coefficient_friction"]
                           .get<double>();
        }

        printer.OpenElement("worldbody");
        std::stringstream ss;
        Eigen::Matrix3d swapYZ;
        swapYZ << -1, 0, 0, 0, 0, 1, 0, 1, 0;
        rbId = 0;
        for (const auto& rbI : scene["rigid_body_problem"]["rigid_bodies"]) {
            printer.OpenElement("body");

            if (rbI.find("position") != rbI.end()) {
                ss << -rbI["position"][0].get<double>() << " "
                   << rbI["position"][2] << " " << rbI["position"][1];
                printer.PushAttribute("pos", ss.str().c_str());
                ss.str(std::string());
            } else {
                printer.PushAttribute("pos", "0 0 0");
            }

            std::vector<double> rotation;
            if (rbI.find("rotation") != rbI.end()) {
                rotation = rbI["rotation"].get<std::vector<double>>();
            } else {
                rotation.resize(3, 0.0);
            }
            Eigen::Quaterniond q;
            q = swapYZ
                * Eigen::AngleAxisd(
                      rotation[2] * igl::PI / 180.0, Eigen::Vector3d::UnitZ())
                * Eigen::AngleAxisd(
                      rotation[1] * igl::PI / 180.0, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(
                      rotation[0] * igl::PI / 180.0, Eigen::Vector3d::UnitX());
            ss << q.w() << " " << q.x() << " " << q.y() << " " << q.z();
            printer.PushAttribute("quat", ss.str().c_str());
            ss.str(std::string());

            bool isFixed = false;
            if (rbI.find("is_dof_fixed") != rbI.end()) {
                try {
                    for (bool i :
                         rbI["is_dof_fixed"].get<std::vector<bool>>()) {
                        isFixed |= i;
                    }
                } catch (...) {
                    isFixed = rbI["is_dof_fixed"].get<bool>();
                }
                if (isFixed && mode == 0) {
                    printer.OpenElement("inertial");
                    printer.PushAttribute("mass", "0");
                    printer.CloseElement(); // inertial
                }
            }

            if (!isFixed && rbI.find("type") != rbI.end()) {
                if (rbI["type"].get<std::string>() == "static") {
                    isFixed = true;
                    if (mode == 0) {
                        printer.OpenElement("inertial");
                        printer.PushAttribute("mass", "0");
                        printer.CloseElement(); // inertial
                    }
                }
            }

            printer.OpenElement("geom");
            printer.PushAttribute("mesh", meshNames[rbId].c_str());
            printer.PushAttribute("type", "mesh");
            if (mode == 0) {
                if (rbI.find("scale") != rbI.end()) {
                    try {
                        double scale = rbI["scale"].get<double>();
                        std::string sizeStr(std::to_string(scale));
                        sizeStr += " " + sizeStr + " " + sizeStr;
                        printer.PushAttribute("size", sizeStr.c_str());
                    } catch (...) {
                        std::vector<double> scale =
                            rbI["scale"].get<std::vector<double>>();
                        if (scale.size() != 3) {
                            spdlog::error("scale dimension error!");
                            exit(-1);
                        }
                        for (auto& i : scale) {
                            if (!i) {
                                i = 1;
                            }
                        }
                        std::string sizeStr(std::to_string(scale[0]));
                        sizeStr += " " + std::to_string(scale[1]) + " "
                            + std::to_string(scale[2]);
                        printer.PushAttribute("size", sizeStr.c_str());
                    }
                } else if (rbI.find("dimensions") != rbI.end()) {
                    try {
                        double scale = rbI["dimensions"].get<double>();
                        std::string sizeStr(std::to_string(scale));
                        sizeStr += " " + sizeStr + " " + sizeStr;
                        printer.PushAttribute("size", sizeStr.c_str());
                    } catch (...) {
                        std::vector<double> scale =
                            rbI["dimensions"].get<std::vector<double>>();
                        if (scale.size() != 3) {
                            spdlog::error("scale dimension error!");
                            exit(-1);
                        }
                        for (auto& i : scale) {
                            if (!i) {
                                i = 1;
                            }
                        }
                        std::string sizeStr(std::to_string(scale[0]));
                        sizeStr += " " + std::to_string(scale[1]) + " "
                            + std::to_string(scale[2]);
                        printer.PushAttribute("size", sizeStr.c_str());
                    }
                }
            }
            if (rbI.find("density") != rbI.end()) {
                printer.PushAttribute(
                    "density",
                    std::to_string(rbI["density"].get<double>()).c_str());
            }
            if (mode == 0) {
                // Bullet
                printer.PushAttribute(
                    "friction", std::to_string(std::sqrt(fricCoef)).c_str());
            } else if (mode == 1) {
                // Mujoco
                printer.PushAttribute(
                    "friction", std::to_string(fricCoef).c_str());
            }
            printer.CloseElement(); // geom

            if (mode == 1 && !isFixed) {
                printer.OpenElement("freejoint");
                printer.CloseElement(); // freejoint
            }

            printer.CloseElement(); // body
            ++rbId;
        }
        printer.CloseElement(); // worldbody

        printer.CloseElement(); // mujoco

        tinyxml2::XMLDocument doc;
        doc.Print(&printer);

        fclose(output);
        input.close();
        return 0;
    } else {
        spdlog::error("Json file not found!");
    }
    return -1;
}

int main(int argc, char* argv[])
{
    ipc::rigid::set_logger_level(spdlog::level::info);
    CLI::App app("transform JSON to MJCF");

    std::string json_path;
    app.add_option(
           "json_path", json_path, "JSON file with input scene to convert")
        ->required();

    int mode = 0;
    app.add_option("mode", mode, "conversion mode 0: Bullet, 1: MuJoCo");

    CLI11_PARSE(app, argc, argv);

    return json_to_mjcf(json_path, mode);
}
