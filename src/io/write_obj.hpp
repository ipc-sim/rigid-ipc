// Modified version of read_obj from libigl to include reading polyline elements
// as edges.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <vector>

#include <Eigen/Core>

namespace ccd {
namespace io {

    bool write_obj(
        const std::string str,
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const std::vector<size_t>& codim_edges_to_edges,
        const Eigen::MatrixXi& F);

} // namespace io
} // namespace ccd
