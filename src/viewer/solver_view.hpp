#pragma once

#include <opt/linearized_constraint_solver.hpp>
#include <opt/barrier_newton_solver.hpp>

namespace ccd {
void linearized_constraint_solver_view(ccd::opt::LinearizedCstrSolver& solver);
void barrier_newton_solver_view(ccd::opt::BarrierNewtonSolver& solver);
}
