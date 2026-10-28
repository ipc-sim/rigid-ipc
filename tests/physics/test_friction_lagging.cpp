#include <catch2/catch.hpp>

#include <nlohmann/json.hpp>

#include <SimState.hpp>
#include <physics/rigid_body_problem.hpp>

using namespace ipc::rigid;

// A cube resting on a slab, held by static friction. The Newton solver stops
// on its default velocity tolerance, which leaves a momentum balance far above
// the requested tolerance. Lagging friction "until converged" must still end
// each step instead of re-solving forever.
TEST_CASE(
    "Lagging friction to an unreachable tolerance still finishes each step",
    "[friction][RB-Problem]")
{
    const double gap = 1.1e-3; // just outside the barrier activation distance
    SimState state;
    REQUIRE(state.init(
        nlohmann::json {
            { "timestep", 0.01 },
            { "friction_constraints",
              { { "iterations", -1 }, { "momentum_balance_tol", 1e-12 } } },
            { "rigid_body_problem",
              { { "coefficient_friction", 0.5 },
                { "gravity", { 0, -9.81, 0 } },
                { "rigid_bodies",
                  { { { "mesh", "cube.obj" },
                      { "position", { 0, -0.5, 0 } },
                      { "scale", { 10, 1, 10 } },
                      { "type", "static" } },
                    { { "mesh", "cube.obj" },
                      { "position", { 0, 0.5 + gap, 0 } } } } } } },
        }));
    auto problem =
        std::dynamic_pointer_cast<RigidBodyProblem>(state.problem_ptr);
    REQUIRE(problem != nullptr);

    const Eigen::VectorXd start = problem->m_assembler[1].pose.position;
    for (int i = 0; i < 10; i++) {
        state.simulation_step();
    }

    // The cube settles onto the slab without sliding sideways.
    const Eigen::VectorXd end = problem->m_assembler[1].pose.position;
    CHECK(std::abs(end.x() - start.x()) < 1e-6);
    CHECK(std::abs(end.z() - start.z()) < 1e-6);
    CHECK(std::abs(end.y() - start.y()) < 2 * gap);
}
