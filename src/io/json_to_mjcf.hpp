#pragma once

#include <nlohmann/json.hpp>
using nlohmann::json;

#include <tinyxml2.h>

#include <sys/stat.h> // for mkdir

int json_to_mjcf(
    const char* jsonFilePath, int mode = 0) // mode0: bullet, mode1: mujoco
{
    if (mode > 1) {
        spdlog::error("Invalid mode in mjcf file generation!");
        exit(-1);
    }

    std::ifstream input(jsonFilePath);
    if (input.good()) {
        json scene = json::parse(input, nullptr, false);
        if (scene.is_discarded()) {
            spdlog::error("Invalid Json file");
            input.close();
            return -1;
        }

        std::string jsonFilePathStr(jsonFilePath);
        std::string xmlFilePath =
            jsonFilePathStr.substr(0, jsonFilePathStr.find_last_of('.'))
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
            if (mode == 1 && (rbI.find("dimensions") != rbI.end() || rbI.find("scale") != rbI.end())) {
                printer.OpenElement("mesh");
                printer.PushAttribute("file", (meshNames.back() + ".stl").c_str());
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
                        for (auto &i : scale) {
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
                else if (rbI.find("dimensions") != rbI.end()) {
                    try {
                        double scale = rbI["dimensions"].get<double>();
                        std::string sizeStr(std::to_string(scale));
                        sizeStr += " " + sizeStr + " " + sizeStr;
                        printer.PushAttribute("scale", sizeStr.c_str());
                    }
                    catch(...) {
                        std::vector<double> scale = rbI["dimensions"].get<std::vector<double>>();
                        if (scale.size() != 3) {
                            spdlog::error("scale dimension error!");
                            exit(-1);
                        }
                        for (auto &i : scale) {
                            if (!i) {
                                i = 1;
                            }
                        }
                        std::string sizeStr(std::to_string(scale[0]));
                        sizeStr += " " + std::to_string(scale[1]) + " " + std::to_string(scale[2]);
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
                      rotation[2] * M_PI / 180.0, Eigen::Vector3d::UnitZ())
                * Eigen::AngleAxisd(
                      rotation[1] * M_PI / 180.0, Eigen::Vector3d::UnitY())
                * Eigen::AngleAxisd(
                      rotation[0] * M_PI / 180.0, Eigen::Vector3d::UnitX());
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
                        for (auto &i : scale) {
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
                else if (rbI.find("dimensions") != rbI.end()) {
                    try {
                        double scale = rbI["dimensions"].get<double>();
                        std::string sizeStr(std::to_string(scale));
                        sizeStr += " " + sizeStr + " " + sizeStr;
                        printer.PushAttribute("size", sizeStr.c_str());
                    }
                    catch(...) {
                        std::vector<double> scale = rbI["dimensions"].get<std::vector<double>>();
                        if (scale.size() != 3) {
                            spdlog::error("scale dimension error!");
                            exit(-1);
                        }
                        for (auto &i : scale) {
                            if (!i) {
                                i = 1;
                            }
                        }
                        std::string sizeStr(std::to_string(scale[0]));
                        sizeStr += " " + std::to_string(scale[1]) + " " + std::to_string(scale[2]);
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
                printer.PushAttribute("friction", std::to_string(std::sqrt(fricCoef)).c_str());
            }
            else if (mode == 1) {
                // Mujoco
                printer.PushAttribute("friction", std::to_string(fricCoef).c_str());
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
    }
    else {
        spdlog::error("Json file not found!");
    }
    return -1;
}

int generate_bullet_results(
    const char* jsonFilePath, const char* bulletResultFilePath)
{
    std::ifstream input(jsonFilePath);
    std::ifstream bulletTransFile(bulletResultFilePath);
    if (input.good() && bulletTransFile.good()) {
        json scene = json::parse(input, nullptr, false);
        if (scene.is_discarded()) {
            spdlog::error("Invalid Json file");
            input.close();
            return -1;
        }
        mkdir("output", 0777);
        std::string bulletResultFilePathStr(bulletResultFilePath);
        int lastSlash = bulletResultFilePathStr.find_last_of('/');
        std::string folderPath = "output/"
            + bulletResultFilePathStr.substr(
                  lastSlash,
                  bulletResultFilePathStr.find_last_of('.') - lastSlash)
            + "/";
        mkdir(folderPath.c_str(), 0777);

        std::vector<Eigen::MatrixXd> V;
        std::vector<Eigen::MatrixXi> F;
        for (const auto& rbI : scene["rigid_body_problem"]["rigid_bodies"]) {
            std::string meshFilePath = rbI["mesh"].get<std::string>();
            V.resize(V.size() + 1);
            F.resize(F.size() + 1);
            igl::readOBJ("meshes/" + meshFilePath, V.back(), F.back());
            if (rbI.find("scale") != rbI.end()) {
                try {
                    double scale = rbI["scale"].get<double>();
                    V.back() *= scale;
                } catch (...) {
                    std::vector<double> scale =
                        rbI["scale"].get<std::vector<double>>();
                    if (scale.size() != 3) {
                        spdlog::error("scale dimension error!");
                        exit(-1);
                    }
                    V.back().col(0) *= scale[0];
                    V.back().col(1) *= scale[1];
                    V.back().col(2) *= scale[2];
                }
            }
            else if (rbI.find("dimensions") != rbI.end()) {
                try {
                    double scale = rbI["dimensions"].get<double>();
                    V.back() *= scale;
                }
                catch(...) {
                    std::vector<double> scale = rbI["dimensions"].get<std::vector<double>>();
                    if (scale.size() != 3) {
                        spdlog::error("scale dimension error!");
                        exit(-1);
                    }
                    V.back().col(0) *= scale[0];
                    V.back().col(1) *= scale[1];
                    V.back().col(2) *= scale[2];
                }
            }
        }
        input.close();

        Eigen::Matrix3d swapYZ;
        swapYZ << -1, 0, 0, 0, 0, 1, 0, 1, 0;
        int lastFrameI = -1, VStart = 0;
        FILE* objFile = NULL;
        while (!bulletTransFile.eof()) {
            int frameI, objI;
            bulletTransFile >> frameI >> objI;
            if (frameI == -1) {
                // timing output
                spdlog::info("timing: {:g} min", objI / 60.0);
                exit(0);
            }
            if (frameI != lastFrameI) {
                lastFrameI = frameI;
                VStart = 0;

                if (objFile) {
                    fclose(objFile);
                }
                objFile = fopen(
                    (folderPath + std::to_string(frameI) + ".obj").c_str(),
                    "w");
            }

            Eigen::Vector3d trans;
            bulletTransFile >> trans[0] >> trans[1] >> trans[2];

            double x, y, z, w;
            bulletTransFile >> x >> y >> z >> w;
            Eigen::Matrix3d rotMtr =
                Eigen::Quaterniond(w, x, y, z).toRotationMatrix();

            for (int vI = 0; vI < V[objI].rows(); ++vI) {
                Eigen::Vector3d v = swapYZ
                    * (rotMtr * V[objI].row(vI).transpose()
                       + trans); // center of mass?
                fprintf(objFile, "v %le %le %le\n", v[0], v[1], v[2]);
            }
            fprintf(objFile, "g obj%d\n", objI);
            for (int fI = 0; fI < F[objI].rows(); ++fI) {
                fprintf(
                    objFile, "f %d %d %d\n", F[objI].row(fI)[0] + 1 + VStart,
                    F[objI].row(fI)[1] + 1 + VStart,
                    F[objI].row(fI)[2] + 1 + VStart);
            }
            VStart += V[objI].rows();
        }
        if (objFile) {
            fclose(objFile);
        }
        bulletTransFile.close();
        return 0;
    } else {
        if (!input.good()) {
            spdlog::error("json file open error!");
        }
        if (!bulletTransFile.good()) {
            spdlog::error("bullet file open error!");
        }
    }
    return -1;
}
