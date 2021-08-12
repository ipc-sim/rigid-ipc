#include <CLI/CLI.hpp>
#include <ghc/fs_std.hpp> // filesystem
#include <igl/readOBJ.h>
#include <nlohmann/json.hpp>

#include <logger.hpp>

int generate_bullet_results(
    const std::string& json_path, const std::string& bullet_result_file_path)
{
    std::ifstream input(json_path);
    std::ifstream bulletTransFile(bullet_result_file_path);
    if (input.good() && bulletTransFile.good()) {
        nlohmann::json scene = nlohmann::json::parse(input, nullptr, false);
        if (scene.is_discarded()) {
            spdlog::error("Invalid Json file");
            input.close();
            return -1;
        }
        fs::create_directories(fs::path("output"));
        int lastSlash = bullet_result_file_path.find_last_of('/');
        std::string folderPath = "output/"
            + bullet_result_file_path.substr(
                  lastSlash,
                  bullet_result_file_path.find_last_of('.') - lastSlash)
            + "/";
        fs::create_directories(fs::path(folderPath));

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
            } else if (rbI.find("dimensions") != rbI.end()) {
                try {
                    double scale = rbI["dimensions"].get<double>();
                    V.back() *= scale;
                } catch (...) {
                    std::vector<double> scale =
                        rbI["dimensions"].get<std::vector<double>>();
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

int main(int argc, char* argv[])
{
    ipc::rigid::set_logger_level(spdlog::level::info);
    CLI::App app("generating bullet output geometry from transformation");

    std::string json_path;
    app.add_option("json_path", json_path, "JSON file with input scene")
        ->required();

    std::string bullet_result_path;
    app.add_option(
           "bullet_result_path", bullet_result_path, "bullet output file path")
        ->required();

    CLI11_PARSE(app, argc, argv);

    return generate_bullet_results(json_path, bullet_result_path);
}
