#include "read_scene.hpp"

#include <fstream>

namespace ccd {
namespace io {
    void read_scene(const std::string filename, Eigen::MatrixX2d& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixX2d& displacements)
    {
        using nlohmann::json;

        std::ifstream input(filename);
        json scene;
        input >> scene;

        return read_scene(scene, vertices, edges, displacements);
    }

    void read_scene_from_str(const std::string str, Eigen::MatrixX2d& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixX2d& displacements)
    {
        using nlohmann::json;
        json scene = json::parse(str.c_str());
        return read_scene(scene, vertices, edges, displacements);
    }

    void read_scene(const nlohmann::json scene, Eigen::MatrixX2d& vertices, Eigen::MatrixX2i& edges, Eigen::MatrixX2d& displacements)
    {
        vec_vec_d vec_vertices = scene["vertices"].get<vec_vec_d>();
        size_t num_vertices = vec_vertices.size();
        vertices.resize(long(num_vertices), 2);

        for (size_t i = 0; i < num_vertices; ++i) {
            vertices.row(int(i)) << vec_vertices[i][0], vec_vertices[i][1];
        }

        vec_vec_i vec_edges = scene["edges"].get<vec_vec_i>();
        size_t num_edges = vec_edges.size();
        edges.resize(long(num_edges), 2);

        for (size_t i = 0; i < num_edges; ++i) {
            edges.row(int(i)) << vec_edges[i][0], vec_edges[i][1];
        }

        vec_vec_d vec_displacements = scene["displacements"].get<vec_vec_d>();
        size_t num_displacements = vec_displacements.size();
        assert(num_displacements == num_vertices);
        displacements.resize(long(num_displacements), 2);

        for (size_t i = 0; i < num_displacements; ++i) {
            displacements.row(int(i)) << vec_displacements[i][0], vec_displacements[i][1];
        }
    }

}
}
