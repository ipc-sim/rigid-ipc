#include "io/write_scene.hpp"

#include <fstream>
#include <iomanip> // std::setw

#include <io/serialize_json.hpp>
namespace ccd {
namespace io {

    void write_scene(const std::string filename,
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXd& displacements)
    {
        using nlohmann::json;
        json scene = write_scene(vertices, edges, displacements);

        std::ofstream o(filename);
        o << std::setw(4) << scene << std::endl;
    }

    nlohmann::json write_scene(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXd& displacements)
    {
        using nlohmann::json;
        json scene;
        scene["vertices"] = to_json<double>(vertices);
        scene["displacements"] = to_json<double>(displacements);
        scene["edges"] = to_json<int>(edges);
        return scene;
    }
} // namespace io
} // namespace ccd
