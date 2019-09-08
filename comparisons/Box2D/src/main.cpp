/**
 * @brief Read a JSON fixture file and simulate the fixture in Box2D.
 * @author Zachary Ferguson
 */

#include <Box2D/Box2D.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <fmt/format.h>
#include <nlohmann/json.hpp>
#include <spdlog/spdlog.h>

typedef std::array<double, 2> array2d;
typedef std::array<double, 3> array3d;
typedef std::array<size_t, 2> array2i;
typedef std::array<bool, 3> array3b;
namespace Eigen {
typedef Translation<double, 2> Translation2Dd;
}

/// @brief Convert a b2Vec2 to an Eigen::Vector2d.
inline Eigen::Vector2d b2Vec2_to_Vector2d(const b2Vec2& v)
{
    return Eigen::Vector2d(v.x, v.y);
}

/// @brief Convert an Eigen::Vector2d to a std::array<double, 2>.
inline array2d Vector2d_to_array2d(const Eigen::Vector2d& v)
{
    return { { v.x(), v.y() } };
}

/// @brief Convert an Eigen::Vector2d to a std::array<double, 2>.
inline b2Vec2 Vector2d_to_b2Vec2(const Eigen::Vector2d& v)
{
    return b2Vec2(v.x(), v.y());
}

/// @brief Convert a b2Vec2 to a std::array<double, 2>.
inline array2d b2Vec2_to_array2d(const b2Vec2& v) { return { { v.x, v.y } }; }

/// @brief Save the edges of the bodies.
void save_edges(const b2World& world, std::vector<array2i>& edges)
{
    size_t starting_vertex_index = 0;

    for (const b2Body* body = world.GetBodyList(); body != nullptr;
         body = body->GetNext()) {
        for (const b2Fixture* fixture = body->GetFixtureList();
             fixture != nullptr; fixture = fixture->GetNext()) {
            // Assume all shapes are polygons
            const b2Shape* shape = fixture->GetShape();
            assert(shape->GetType() == b2Shape::Type::e_polygon);
            const b2PolygonShape* polygon
                = dynamic_cast<const b2PolygonShape*>(shape);

            for (size_t i = 0; i < polygon->m_count; i++) {
                edges.push_back({ { starting_vertex_index + i,
                    starting_vertex_index + (i + 1) % polygon->m_count } });
            }

            starting_vertex_index += polygon->m_count;
        }
    }
}

/// @brief Save the group_id of the vertices.
void save_group_id(const b2World& world, std::vector<size_t>& group_ids)
{
    int group_id = 0;
    for (const b2Body* body = world.GetBodyList(); body != nullptr;
         body = body->GetNext()) {
        if (!(body->GetFixtureList())) {
            continue;
        }
        // Save the group_id per vertex
        for (const b2Fixture* fixture = body->GetFixtureList();
             fixture != nullptr; fixture = fixture->GetNext()) {
            // Assume all shapes are polygons
            const b2Shape* shape = fixture->GetShape();
            assert(shape->GetType() == b2Shape::Type::e_polygon);
            const b2PolygonShape* polygon
                = dynamic_cast<const b2PolygonShape*>(shape);

            for (size_t i = 0; i < polygon->m_count; i++) {
                group_ids.push_back(group_id);
            }
        }
        group_id++;
    }
}

/// @brief Save the state of the simulation at the current step.
void save_step(const b2World& world,
    std::vector<std::vector<array2d>>& vertices_sequence,
    std::vector<nlohmann::json>& state_sequence)
{
    // State variables
    std::vector<array2d> vertices;
    double angular_momentum = 0;
    b2Vec2 linear_momentum(0, 0);
    double kinetic_energy = 0;
    double potential_energy = 0;
    std::vector<nlohmann::json> body_states;

    // Save gravity to compute potential_energy
    const Eigen::Vector2d gravity = b2Vec2_to_Vector2d(world.GetGravity());

    for (const b2Body* body = world.GetBodyList(); body != nullptr;
         body = body->GetNext()) {

        if (!(body->GetFixtureList())) {
            continue;
        }

        Eigen::Vector2d position = b2Vec2_to_Vector2d(body->GetPosition());
        double theta = body->GetAngle();
        Eigen::Transform<double, 2, Eigen::Affine> T
            = Eigen::Translation2Dd(position) * Eigen::Rotation2Dd(theta);

        // Save the vertices
        for (const b2Fixture* fixture = body->GetFixtureList();
             fixture != nullptr; fixture = fixture->GetNext()) {
            // Assume all shapes are polygons
            const b2Shape* shape = fixture->GetShape();
            assert(shape->GetType() == b2Shape::Type::e_polygon);
            const b2PolygonShape* polygon
                = dynamic_cast<const b2PolygonShape*>(shape);

            for (size_t i = 0; i < polygon->m_count; i++) {
                // Transform the vertices to the world coordinates
                vertices.push_back(Vector2d_to_array2d(
                    (T * b2Vec2_to_Vector2d(polygon->m_vertices[i]))));
            }
        }

        // Save the state
        // Skip fixed bodys from the energy computation
        if (body->GetType() == b2_dynamicBody) {
            linear_momentum += body->GetMass() * body->GetLinearVelocity();
            angular_momentum += body->GetInertia() * body->GetAngularVelocity();
            kinetic_energy += 0.5 * body->GetMass()
                    * body->GetLinearVelocity().LengthSquared()
                + 0.5 * body->GetInertia() * body->GetAngularVelocity()
                    * body->GetAngularVelocity();
            Eigen::Vector2d world_com
                = b2Vec2_to_Vector2d(body->GetWorldCenter());
            potential_energy += -body->GetMass() * gravity.dot(world_com);
        }
        nlohmann::json body_state;
        body_state["position"] = { position.x(), position.y(), theta };
        body_state["velocity"] = { body->GetLinearVelocity().x,
            body->GetLinearVelocity().y, body->GetAngularVelocity() };
        body_states.push_back(body_state);
    }

    nlohmann::json state;
    state["angular_momentum"] = angular_momentum;
    state["linear_momentum"] = b2Vec2_to_array2d(linear_momentum);
    state["kinetic_energy"] = kinetic_energy;
    state["potential_energy"] = potential_energy;
    state["rigid_bodies"] = body_states;
    state["min_distance"] = nullptr;
    vertices_sequence.push_back(vertices);
    state_sequence.push_back(state);
}

void load_rigid_bodies(const nlohmann::json& body_args,
    const double friction_coeff,
    const double restitution_coeff,
    b2World& world)
{
    for (auto& jrb : body_args) {
        nlohmann::json args = R"({
              "polygons":[],
              "density":1.0,
              "is_dof_fixed":[false,false,false],
              "position":[0.0,0.0],
              "theta":0.0,
              "velocity":[0.0,0.0,0.0]
              })"_json;
        args.merge_patch(jrb);

        // Define the ground
        b2BodyDef body_def;

        // Is the body fixed?
        array3b is_dof_fixed = args["is_dof_fixed"].get<array3b>();
        bool is_fixed = (is_dof_fixed[0] || is_dof_fixed[1]) && is_dof_fixed[2];
        body_def.type = is_fixed ? b2_staticBody : b2_dynamicBody;
        body_def.fixedRotation = is_dof_fixed[2];

        body_def.bullet = body_def.type == b2_dynamicBody;
        body_def.allowSleep = false;
        body_def.awake = true;

        // Set the position and rotation
        array2d position = args["position"].get<array2d>();
        body_def.position.Set(position[0], position[1]);
        body_def.angle = args["theta"].get<double>() * M_PI / 180.0;
        // Set the linear and angular velocity
        array3d velocity = args["velocity"].get<array3d>(); // vx, vy, vω
        body_def.linearVelocity.Set(velocity[0], velocity[1]);
        body_def.angularVelocity = velocity[2];
        body_def.linearDamping = 0;
        body_def.angularDamping = 0;

        // Create the ground
        b2Body* body = world.CreateBody(&body_def);

        // Add revolute join to constrain position but not rotation
        if ((is_dof_fixed[0] || is_dof_fixed[1]) && !is_dof_fixed[2]) {
            b2BodyDef empty_body_def;
            empty_body_def.position.Set(position[0], position[1]);
            empty_body_def.type = b2_staticBody;
            b2Body* empty_body = world.CreateBody(&empty_body_def);
            b2RevoluteJointDef jointDef;
            jointDef.Initialize(empty_body, body, empty_body->GetWorldCenter());
            world.CreateJoint(&jointDef);
        }

        double density = args["density"].get<double>();

        // Example JSON Polygon:
        // "polygons": [
        //     [[0, 0], [1, 0], [1, 1], [0, 1]],
        //     [[-1, 0], [0, 0], [0, 1], [-1, 1]]
        // ]
        std::vector<std::vector<array2d>> polygons
            = args["polygons"].get<std::vector<std::vector<array2d>>>();
        std::vector<b2Vec2> points;
        points.reserve(b2_maxPolygonVertices);
        for (const auto& polygon : polygons) {
            if (polygon.size() > b2_maxPolygonVertices) {
                spdlog::error(
                    "input {}-gon exceeds maximum number of vertices in Box2D ({})",
                    polygon.size(), b2_maxPolygonVertices);
                exit(1);
            }

            for (const auto& point : polygon) {
                points.push_back(b2Vec2(point[0], point[1]));
            }

            // Define the shape
            b2PolygonShape shape;
            shape.Set(points.data(), points.size());
            // Add the shape fixtures to the line body
            b2Fixture* fixture = body->CreateFixture(&shape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);

            points.clear();
        }

        // Compute the mass, center of mass, and moment of intertia
        // {
        //     Eigen::MatrixXd vertices;
        //     from_json<double>(args["vertices"], vertices);
        //     Eigen::MatrixXi edges;
        //     from_json<int>(args["edges"], edges);
        //
        //     Eigen::Vector2d com = center_of_mass(vertices, edges);
        //     vertices.rowwise() -= com.transpose();
        //
        //     b2MassData mass_data;
        //     body->GetMassData(&mass_data);
        //     // mass_data.center = Vector2d_to_b2Vec2(com);
        //
        //     Eigen::VectorXd masses;
        //     from_json<double>(args["masses"], masses);
        //     if (masses.size() != vertices.rows()) {
        //         mass_vector(vertices, edges, masses);
        //     }
        //     mass_data.mass = masses.sum();
        //     mass_data.I = moment_of_inertia(vertices, masses);
        //     body->SetMassData(&mass_data);
        // }

        spdlog::info(
            "mass={:g} inertia={:g}", body->GetMass(), body->GetInertia());
    }
}

void load_scene(nlohmann::json& args,
    b2World& world,
    double& timestep,
    int& num_steps,
    int& num_iterations)
{
    {
        nlohmann::json full_args = R"({
            "timestep_size": 1e-2,
            "max_iterations":1000,
            "rigid_body_problem":{
                "coefficient_restitution":0.0,
                "coefficient_friction":0.0,
                "gravity":[0.0,0.0,0.0],
                "rigid_bodies": []
            }
        })"_json;

        full_args.merge_patch(args);
        args = full_args;
    }

    timestep = args["timestep_size"].get<double>();
    num_steps = args["max_iterations"].get<double>();
    num_iterations = 10000;

    double friction_coeff
        = args["rigid_body_problem"]["coefficient_friction"].get<double>();
    double restitution_coeff = std::max(
        args["rigid_body_problem"]["coefficient_restitution"].get<double>(),
        0.0);

    array3d gravity_array
        = args["rigid_body_problem"]["gravity"].get<array3d>();
    world.SetGravity(b2Vec2(gravity_array[0], gravity_array[1]));

    load_rigid_bodies(args["rigid_body_problem"]["rigid_bodies"],
        friction_coeff, restitution_coeff, world);
}

nlohmann::json load_simulation_args(const std::string& filename)
{
    std::ifstream input(filename);
    if (input.good()) {
        nlohmann::json scene = nlohmann::json::parse(input, nullptr, false);
        if (!scene.is_discarded()) {
            if (scene.find("args") != scene.end()) {
                if (scene["args"] != nullptr) {
                    return scene["args"];
                }
            } else {
                return scene;
            }
        }
    }
    throw fmt::format("Invalid simulation file: {}", filename);
}

void save_simulation_results(const std::string& filename,
    const nlohmann::json& args,
    const std::vector<array2i>& edges,
    const std::vector<size_t>& group_id,
    const std::vector<std::vector<array2d>>& vertices_sequence,
    const std::vector<nlohmann::json>& state_sequence)
{
    nlohmann::json results;
    results["args"] = args;
    results["active_args"] = nlohmann::json();
    results["animation"] = nlohmann::json();
    results["animation"]["edges"] = edges;
    results["animation"]["group_id"] = group_id;
    results["animation"]["vertices_sequence"] = vertices_sequence;
    results["animation"]["state_sequence"] = state_sequence;

    std::ofstream o(filename);
    o << std::setw(4) << results << std::endl;
    spdlog::info("simulation results saved to {}", filename);
}

struct CLIArgs {
    std::string scene_path = "";
    std::string output_dir = "";
    int num_steps = -1;
    bool is_log_trace;
};

int parse_command_line_args(int argc, char** argv, CLIArgs& cli_args)
{
    CLI::App app { "run headless comparison" };
    app.add_option("scene_path,-s,--scene-path", cli_args.scene_path,
           "JSON file with input scene")
        ->required();
    app.add_option("output_dir,-o,--output-path", cli_args.output_dir,
           "directory for results")
        ->required();

    app.add_option("--num-steps", cli_args.num_steps, "number of time-steps");
    app.add_flag("--trace", cli_args.is_log_trace, "log everything");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }
    return 0;
}

int main(int argc, char** argv)
{
    spdlog::set_level(spdlog::level::info);

    // Parse the command line arguments
    CLIArgs cli_args;
    int error_code = parse_command_line_args(argc, argv, cli_args);
    if (error_code) {
        return error_code;
    }

    if (cli_args.is_log_trace) {
        spdlog::set_level(spdlog::level::trace);
    }

    // Simulation variables
    b2World world(b2Vec2(0, -9.81));
    double timestep;
    int num_steps;
    int num_iterations;

    nlohmann::json json_args;
    try {
        json_args = load_simulation_args(cli_args.scene_path);
    } catch (const std::string& err) {
        spdlog::error(err);
        exit(1);
    }

    // Results
    std::vector<array2i> edges;
    std::vector<size_t> group_id;
    std::vector<std::vector<array2d>> vertices_sequence;
    std::vector<nlohmann::json> state_sequence;

    load_scene(json_args, world, timestep, num_steps, num_iterations);

    if (cli_args.num_steps >= 0) {
        num_steps = cli_args.num_steps;
    }

    save_edges(world, edges);
    save_group_id(world, group_id);
    save_step(world, vertices_sequence, state_sequence);

    spdlog::info("Running {} steps", num_steps);
    spdlog::info("Starting simulation {}", cli_args.scene_path);

    // Main simulation loop
    for (int i = 0; i < num_steps; i++) {
        world.Step(timestep, num_iterations, num_iterations);
        save_step(world, vertices_sequence, state_sequence);
        spdlog::info("finished step {}", i);
    }

    std::string fout = fmt::format("{}/sim.json", cli_args.output_dir);
    save_simulation_results(
        fout, json_args, edges, group_id, vertices_sequence, state_sequence);

    spdlog::info(
        "To postprocess run:\n `python tools/results_to_vtk_files.py {} {}`",
        fout, cli_args.output_dir);
    return 0;
}
