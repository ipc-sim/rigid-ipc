#pragma once

#include <solvers/barrier_solver.hpp>
#include <solvers/ipopt_solver.hpp>
#include <solvers/ncp_solver.hpp>
#include <solvers/nlopt_solver.hpp>
#include <solvers/qp_solver.hpp>

namespace ccd {
void solver_menu(ccd::opt::NLOptSolver& solver);
void solver_menu(ccd::opt::QPSolver& solver);
void solver_menu(ccd::opt::BarrierSolver& solver);
void solver_menu(ccd::opt::NewtonSolver& solver);
void solver_menu(ccd::opt::BFGSSolver& solver);
void solver_menu(ccd::opt::GradientDescentSolver& solver);
void solver_menu(ccd::opt::NCPSolver& ncp_solver);
#ifdef BUILD_WITH_IPOPT
void solver_menu(ccd::opt::IpoptSolver& solver);
#endif

} // namespace ccd
