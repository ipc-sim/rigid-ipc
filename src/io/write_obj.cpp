// Modified version of read_obj from libigl to include writing polyline elements
// as codimensional edges.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "write_obj.hpp"

#include <cassert>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

#include <ghc/fs_std.hpp> // filesystem

#include <logger.hpp>
#include <physics/rigid_body_problem.hpp>

namespace ipc::rigid {

static const Eigen::IOFormat vertices_format(
    Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "v ", "", "", "\n");
static const Eigen::IOFormat faces_format(
    Eigen::FullPrecision, Eigen::DontAlignCols, " ", "\n", "f ", "", "", "\n");

bool write_obj(
    const std::string str,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const std::vector<size_t>& codim_edges_to_edges,
    const Eigen::MatrixXi& F)
{
    using namespace std;
    using namespace Eigen;
    ofstream s(str);
    if (!s.is_open()) {
        spdlog::error("IOError: write_obj() could not open {}", str);
        return false;
    }
    s << V.format(vertices_format) << (F.array() + 1).format(faces_format);
    for (const size_t& ei : codim_edges_to_edges) {
        s << fmt::format("l {:d} {:d}\n", E(ei, 0) + 1, E(ei, 1) + 1);
    }
    return true;
}

bool write_obj(
    const std::string str, const SimulationProblem& problem, bool write_mtl)
{
    std::ofstream s(str);
    if (!s.is_open()) {
        spdlog::error("IOError: write_obj() could not open {}", str);
        return false;
    }

    auto p = fs::path(str);
    std::string ps_name =
        (p.parent_path() / ("points-" + p.filename().string())).string();
    std::ofstream ps(ps_name);
    if (!s.is_open()) {
        spdlog::error("IOError: write_obj() could not open {}", ps_name);
        return false;
    }

    s << "mtllib mat.mtl\n";
    ps << "mtllib mat.mtl\n";

    size_t start_vi = 1;
    for (int i = 0; i < problem.num_bodies(); i++) {
        Eigen::MatrixXd V = problem.vertices(i);
        const auto& F = problem.faces(i);
        const auto& E = problem.edges(i);
        if (F.rows() == 0 && E.rows() == 0) {
            ps << fmt::format("o body{0:04d}\nusemtl body{0:04d}\n", i);
            ps << V.format(vertices_format);
        } else {
            s << fmt::format("o body{0:04d}\nusemtl body{0:04d}\n", i);
            s << V.format(vertices_format)
              << (F.array() + start_vi).format(faces_format);
            for (const size_t& ei : problem.codim_edges_to_edges(i)) {
                s << fmt::format(
                    "l {:d} {:d}\n", E(ei, 0) + start_vi, E(ei, 1) + start_vi);
            }
            start_vi += V.rows();
        }
    }

    if (!write_mtl) {
        return true;
    }

    // 0 is static, 1 is kinematic, 2 is dynamic
    Eigen::VectorXi body_types(problem.num_bodies());
    if (problem.is_rb_problem()) {
        const auto& bodies =
            dynamic_cast<const RigidBodyProblem&>(problem).m_assembler;
        for (int i = 0; i < problem.num_bodies(); i++) {
            body_types[i] = int(bodies[i].type);
        }
    } else {
        body_types.setConstant(2);
    }

    // for (int i = 0; i < problem.num_bodies(); i++) {
    //     std::ofstream mtl_ofs(fmt::format("body{:04d}.mtl", i));
    //     if (!mtl_ofs.is_open()) {
    //         spdlog::error(
    //             "IOError: write_obj() could not open body{:04d}.mtl", i);
    //         return false;
    //     }
    // }

    return true;
}

} // namespace ipc::rigid
