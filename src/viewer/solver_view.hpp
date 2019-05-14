#pragma once

#include <opt/barrier_newton_solver.hpp>
#include <opt/ipopt_solver.hpp>
#include <opt/qp_solver.hpp>
#include <opt/ncp_solver.hpp>
#include <opt/nlopt_solver.hpp>

namespace ccd {
void solver_menu(ccd::opt::NLOptSolver& solver);
void solver_menu(ccd::opt::QPSolver& solver);
void solver_menu(ccd::opt::BarrierNewtonSolver& solver);
void solver_menu(ccd::opt::NCPSolver& ncp_solver);
void solver_menu(ccd::opt::IpoptSolver& solver);

} // namespace ccd
