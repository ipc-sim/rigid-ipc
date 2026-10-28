#include <catch2/catch.hpp>

#include <nlohmann/json.hpp>

#include <SimState.hpp>
#include <physics/rigid_body_problem.hpp>

using namespace ipc::rigid;

namespace {

// Two unit cubes slide side by side across a static slab at the same initial
// speed, with the friction pairs given.
nlohmann::json sliding_cubes_scene(const nlohmann::json& friction_pairs)
{
    const double gap = 1.1e-3; // just outside the barrier activation distance
    return nlohmann::json {
        { "timestep", 0.01 },
        { "ipc_solver", { { "velocity_conv_tol", 1e-5 } } },
        { "friction_constraints",
          { { "iterations", -1 }, { "momentum_balance_tol", 1e-4 } } },
        { "rigid_body_problem",
          { { "coefficient_friction", 0.3 },
            { "friction_pairs", friction_pairs },
            { "gravity", { 0, -9.81, 0 } },
            { "rigid_bodies",
              { { { "name", "slab" },
                  { "mesh", "cube.obj" },
                  { "position", { 0, -0.5, 0 } },
                  { "scale", { 10, 1, 10 } },
                  { "type", "static" } },
                { { "name", "free" },
                  { "mesh", "cube.obj" },
                  { "position", { -2, 0.5 + gap, -2 } },
                  { "linear_velocity", { 1, 0, 0 } } },
                { { "name", "gripped" },
                  { "mesh", "cube.obj" },
                  { "position", { -2, 0.5 + gap, 2 } },
                  { "linear_velocity", { 1, 0, 0 } } } } } } },
    };
}

} // namespace

TEST_CASE(
    "Friction pairs override the global coefficient between two bodies",
    "[friction][RB-Problem]")
{
    SimState state;
    REQUIRE(state.init(sliding_cubes_scene(
        { { { "bodies", { "free", "slab" } },
            { "coefficient_friction", 0.0 } } })));
    auto problem =
        std::dynamic_pointer_cast<RigidBodyProblem>(state.problem_ptr);
    REQUIRE(problem != nullptr);

    // Run well past the 0.34 s the gripped cube needs to stop at μ = 0.3,
    // then measure both cubes' speeds over one more step.
    const double dt = 0.01;
    for (int i = 0; i < 70; i++) {
        state.simulation_step();
    }
    const Eigen::VectorXd free_start = problem->m_assembler[1].pose.position;
    const Eigen::VectorXd gripped_start = problem->m_assembler[2].pose.position;
    state.simulation_step();
    const Eigen::VectorXd free_velocity =
        (problem->m_assembler[1].pose.position - free_start) / dt;
    const Eigen::VectorXd gripped_velocity =
        (problem->m_assembler[2].pose.position - gripped_start) / dt;

    CHECK(free_velocity.x() == Approx(1.0).epsilon(1e-2));
    CHECK(gripped_velocity.norm() < 1e-2);
}

TEST_CASE("Friction pairs reject unknown bodies", "[friction][RB-Problem]")
{
    SimState state;
    CHECK_FALSE(state.init(sliding_cubes_scene(
        { { { "bodies", { "free", "nonexistent" } },
            { "coefficient_friction", 0.0 } } })));
}
