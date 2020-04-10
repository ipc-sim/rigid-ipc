#pragma once

#include <nlohmann/json.hpp>
using nlohmann::json;

#include <tinyxml2.h>

int json_to_mjcf(const char* jsonFilePath) 
{
    std::ifstream input(jsonFilePath);
    if (input.good()) {
        json scene = json::parse(input, nullptr, false);
        if (scene.is_discarded()) {
            spdlog::error("Invalid Json file");
            input.close();
            return -1;
        }

        std::string jsonFilePathStr(jsonFilePath);
        std::string xmlFilePath = jsonFilePathStr.substr(0, jsonFilePathStr.find_last_of('.')) + ".xml";
        FILE *output = fopen(xmlFilePath.c_str(), "w");
        if (!output) {
            spdlog::error("failed to create file");
            input.close();
            return -1;
        }

        // transform format
        tinyxml2::XMLPrinter printer( output );
        printer.OpenElement("mujoco");

        std::vector<std::string> meshNames;
        printer.OpenElement("asset");
        for (const auto& rbI : scene["rigid_body_problem"]["rigid_bodies"]) {
            printer.OpenElement("mesh");

            std::string meshFilePath = rbI["mesh"].get<std::string>();
            printer.PushAttribute("file", meshFilePath.c_str());
            meshNames.emplace_back(meshFilePath.substr(0, meshFilePath.find_last_of('.')));
            printer.PushAttribute("name", meshNames.back().c_str());
            
            printer.CloseElement(); // mesh
        }
        printer.CloseElement(); // asset

        printer.OpenElement("worldbody");
        std::stringstream ss;
        Eigen::Matrix3d swapYZ;
        swapYZ << -1, 0, 0,
            0, 0, 1,
            0, 1, 0;
        int rbId = 0;
        for (const auto& rbI : scene["rigid_body_problem"]["rigid_bodies"]) {
            printer.OpenElement("body");

            ss << -rbI["position"][0].get<double>() << " " << rbI["position"][2] << " " << rbI["position"][1];
            printer.PushAttribute("pos", ss.str().c_str());
            ss.str(std::string());

            std::vector<double> rotation;
            if (rbI.find("rotation") != rbI.end()) {
                rotation = rbI["rotation"].get<std::vector<double>>();
            }
            else {
                rotation.resize(3, 0.0);
            }
            Eigen::Quaterniond q;
            q = Eigen::AngleAxisd(rotation[0] * M_PI / 180.0, -Eigen::Vector3d::UnitX())
                    * Eigen::AngleAxisd(rotation[1] * M_PI / 180.0, Eigen::Vector3d::UnitZ())
                    * Eigen::AngleAxisd(rotation[2] * M_PI / 180.0, Eigen::Vector3d::UnitY()) 
                    * swapYZ;
            ss << q.w() << " " << q.x() << " " << q.y() << " " << q.z();
            printer.PushAttribute("quat", ss.str().c_str());
            ss.str(std::string());

            printer.OpenElement("geom");
            printer.PushAttribute("mesh", meshNames[rbId].c_str());
            printer.PushAttribute("type", "mesh");
            printer.CloseElement(); // geom

            if (rbI.find("is_dof_fixed") != rbI.end()) {
                bool isFixed = false;
                for (bool i : rbI["is_dof_fixed"].get<std::vector<bool>>()) {
                    isFixed |= i;
                }
                if (isFixed) {
                    printer.OpenElement("inertial");
                    printer.PushAttribute("mass", "0");
                    printer.CloseElement(); // inertial
                }
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
    return -1;
}

int generate_bullet_results(
    const char* jsonFilePath, 
    const char* bulletResultFilePath) 
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

        std::vector<Eigen::MatrixXd> V;
        std::vector<Eigen::MatrixXi> F;
        for (const auto& rbI : scene["rigid_body_problem"]["rigid_bodies"]) {
            std::string meshFilePath = rbI["mesh"].get<std::string>();
            V.resize(V.size() + 1);
            F.resize(F.size() + 1);
            igl::readOBJ("meshes/" + meshFilePath, V.back(), F.back());
        }
        input.close();

        Eigen::Matrix3d swapYZ;
        swapYZ << -1, 0, 0,
            0, 0, 1,
            0, 1, 0;
        int lastFrameI = -1, VStart = 0;
        FILE *objFile = NULL;
        while (!bulletTransFile.eof()) {
            int frameI, objI;
            bulletTransFile >> frameI >> objI;
            if (frameI != lastFrameI) {
                lastFrameI = frameI;
                VStart = 0;

                if (objFile) {
                    fclose(objFile);
                }
                objFile = fopen((std::to_string(frameI) + ".obj").c_str(), "w");
            }

            Eigen::Vector3d trans;
            bulletTransFile >> trans[0] >> trans[1] >> trans[2];

            double x, y, z, w;
            bulletTransFile >> x >> y >> z >> w;
            Eigen::Matrix3d rotMtr = Eigen::Quaterniond(w, x, y, z).toRotationMatrix();

            for (int vI = 0; vI < V[objI].rows(); ++vI) {
                Eigen::Vector3d v = swapYZ * (rotMtr * V[objI].row(vI).transpose() + trans); // center of mass?
                fprintf(objFile, "v %le %le %le\n", v[0], v[1], v[2]);
            }
            fprintf(objFile, "g obj%d\n", objI);
            for (int fI = 0; fI < F[objI].rows(); ++fI) {
                fprintf(objFile, "f %d %d %d\n", 
                    F[objI].row(fI)[0] + 1 + VStart, 
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
    }
    return -1;
}