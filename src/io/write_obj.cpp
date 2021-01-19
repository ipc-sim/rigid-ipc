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

#include <logger.hpp>

namespace ccd {
namespace io {

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
        s << V.format(IOFormat(
                 FullPrecision, DontAlignCols, " ", "\n", "v ", "", "", "\n"))
          << (F.array() + 1)
                 .format(IOFormat(
                     FullPrecision, DontAlignCols, " ", "\n", "f ", "", "",
                     "\n"));
        for (const size_t& ei : codim_edges_to_edges) {
            s << fmt::format("l {:d} {:d}\n", E(ei, 0) + 1, E(ei, 1) + 1);
        }
        return true;
    }

} // namespace io
} // namespace ccd
