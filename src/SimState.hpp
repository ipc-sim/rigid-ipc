#pragma once

#include <ccd/collision_detection.hpp>

#include <opt/rigid_body_problem.hpp>
#include <opt/barrier_constraint.hpp>
#include <solvers/barrier_solver.hpp>

#include <physics/rigid_body_system.hpp>

namespace ccd {

class SimState {
public:
    SimState();

    void load_scene(const std::string& filename);

    void simulation_step();
    Eigen::MatrixXd solve_collision();

    // CCD
    // ----------------------------------------------
    opt::RigidBodyProblem2 m_problem;
    opt::BarrierSolver m_ccd_solver;
    opt::BarrierConstraint m_ccd_constraint;

    double m_timestep_size;

    bool m_step_had_collision; ///< last step had a collision
    bool m_step_has_collision;///< last step failed to solve collisions
    int m_num_simulation_steps; ///< counts simulation steps
};

} // namespace ccd
