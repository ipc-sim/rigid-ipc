#pragma once

#include <solvers/barrier_newton_solver.hpp>
#include <solvers/ipopt_solver.hpp>
#include <solvers/qp_solver.hpp>
#include <solvers/ncp_solver.hpp>
#include <solvers/nlopt_solver.hpp>

namespace ccd {
void solver_menu(ccd::opt::NLOptSolver& solver);
void solver_menu(ccd::opt::QPSolver& solver);
void solver_menu(ccd::opt::BarrierNewtonSolver& solver);
void solver_menu(ccd::opt::NCPSolver& ncp_solver);
void solver_menu(ccd::opt::IpoptSolver& solver);

} // namespace ccd
