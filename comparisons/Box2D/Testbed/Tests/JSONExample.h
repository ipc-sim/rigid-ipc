#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "json.hpp"

typedef std::array<double, 2> array2d;
typedef std::array<double, 3> array3d;
typedef std::array<size_t, 2> array2i;
typedef std::array<bool, 3> array3b;

const static char* INPUT_FILENAME
    = "../../../fixtures/chain-cross-section.json";

class JSONExample : public Test {
protected:
    nlohmann::json args;
    std::vector<array2i> edges;
    std::vector<std::vector<array2d>> vertices_sequence;
    std::vector<nlohmann::json> state_sequence;

    double timestep;
    long iterations;
    double friction_coeff;
    double restitution_coeff;

    void SaveEdges()
    {
        edges.clear();
        size_t starting_vertex_index = 0;

        for (const b2Body* body = m_world->GetBodyList(); body != nullptr;
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

    static inline b2Vec2 MatrixMultiplication(const b2Mat22& M, const b2Vec2& v)
    {
        return b2Vec2(M.ex.x * v.x + M.ey.x * v.y, M.ex.y * v.x + M.ey.y * v.y);
    }

    static inline double DotProduct(const b2Vec2& u, const b2Vec2& v)
    {
        return u.x * v.x + u.y * v.y;
    }

    static inline b2Mat22 RotationMatrix(const double& theta)
    {
        return b2Mat22(cos(theta), -sin(theta), sin(theta), cos(theta));
    }

    static inline b2Vec2 TransformVector(
        const b2Mat22& R, const b2Vec2& t, const b2Vec2& v)
    {
        return MatrixMultiplication(R, v) + t;
    }

    static inline array2d b2Vec2ToArray(const b2Vec2& v)
    {
        return { { v.x, v.y } };
    }

    void SaveStep()
    {
        std::vector<array2d> vertices;
        double angular_momentum = 0;
        array2d linear_momentum = { { 0, 0 } };
        double kinetic_energy = 0;
        double potential_energy = 0;
        std::vector<nlohmann::json> body_states;

        // TODO: There is an extra rigid body for some reason
        for (const b2Body* body = m_world->GetBodyList(); body != nullptr;
             body = body->GetNext()) {

            const b2Mat22 R = RotationMatrix(body->GetAngle());
            const b2Vec2 pos = body->GetPosition();

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
                    vertices.push_back(b2Vec2ToArray(
                        TransformVector(R, pos, polygon->m_vertices[i])));
                }
            }

            // Save the state
            b2Vec2 body_linear_momentum
                = body->GetMass() * body->GetLinearVelocity();
            linear_momentum[0] += body_linear_momentum.x;
            linear_momentum[1] += body_linear_momentum.y;
            angular_momentum += body->GetInertia() * body->GetAngularVelocity();
            kinetic_energy += 0.5 * body->GetMass()
                    * body->GetLinearVelocity().LengthSquared()
                + 0.5 * body->GetInertia() * body->GetAngularVelocity()
                    * body->GetAngularVelocity();
            potential_energy += -body->GetMass()
                * DotProduct(m_world->GetGravity(), body->GetPosition());

            nlohmann::json body_state;
            body_state["position"] = { body->GetPosition().x,
                body->GetPosition().y, body->GetAngle() };
            body_state["velocity"] = { body->GetLinearVelocity().x,
                body->GetLinearVelocity().y, body->GetAngularVelocity() };
            body_states.push_back(body_state);
        }

        nlohmann::json state;
        state["angular_momentum"] = angular_momentum;
        state["linear_momentum"] = linear_momentum;
        state["kinetic_energy"] = kinetic_energy;
        state["potential_energy"] = potential_energy;
        state["rigid_bodies"] = body_states;
        state["min_distance"] = nullptr;
        vertices_sequence.push_back(vertices);
        state_sequence.push_back(state);
    }

    void LoadRigidBodies(nlohmann::json body_args)
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
            b2BodyDef bodyDef;

            // Is the body fixed?
            array3b is_dof_fixed = args["is_dof_fixed"].get<array3b>();
            bool is_fixed
                = (is_dof_fixed[0] || is_dof_fixed[1]) && is_dof_fixed[2];
            bodyDef.type = is_fixed ? b2_staticBody : b2_dynamicBody;
            bodyDef.fixedRotation = is_dof_fixed[2];

            bodyDef.bullet = bodyDef.type == b2_dynamicBody;
            bodyDef.allowSleep = false;
            bodyDef.awake = true;

            // Set the position and rotation
            array2d position = args["position"].get<array2d>();
            bodyDef.position.Set(position[0], position[1]);
            bodyDef.angle = args["theta"].get<double>() * M_PI / 180.0;
            // Set the linear and angular velocity
            array3d velocity = args["velocity"].get<array3d>(); // vx, vy, vω
            bodyDef.linearVelocity.Set(velocity[0], velocity[1]);
            bodyDef.angularVelocity = velocity[2];
            bodyDef.linearDamping = 0;
            bodyDef.angularDamping = 0;

            // Create the ground
            b2Body* body = m_world->CreateBody(&bodyDef);

            if ((is_dof_fixed[0] || is_dof_fixed[1]) && !is_dof_fixed[2]) {
                b2BodyDef empty_body_def;
                empty_body_def.position.Set(position[0], position[1]);
                empty_body_def.type = b2_staticBody;
                b2Body* empty_body = m_world->CreateBody(&empty_body_def);
                b2RevoluteJointDef jointDef;
                jointDef.Initialize(
                    empty_body, body, empty_body->GetWorldCenter());
                m_world->CreateJoint(&jointDef);
            }

            double density = args["density"].get<double>();

            std::vector<std::vector<array2d>> polygons
                = args["polygons"].get<std::vector<std::vector<array2d>>>();
            for (const auto& polygon : polygons) {
                // Define the ground shape
                b2PolygonShape shape;

                // Example JSON Polygon:
                // "polygons": [
                //     [[0, 0], [1, 0], [1, 1], [0, 1]],
                //     [[-1, 0], [0, 0], [0, 1], [-1, 1]]
                // ]
                std::vector<b2Vec2> points;
                for (const auto& point : polygon) {
                    points.push_back(b2Vec2(point[0], point[1]));
                }

                // Add the shape fixtures to the line body
                shape.Set(points.data(), points.size());
                b2Fixture* fixture = body->CreateFixture(&shape, density);
                fixture->SetFriction(friction_coeff);
                fixture->SetRestitution(restitution_coeff);
            }
        }
    }

    void LoadScene(nlohmann::json args_in)
    {
        // clang-format off
        args = R"({
            "timestep_size": 1e-2,
            "max_iterations":300,
            "rigid_body_problem":{
                "coefficient_restitution":0.0,
                "coefficient_friction":0.0,
                "gravity":[0.0,0.0,0.0],
                "rigid_bodies": []
            }
        })"_json;
        // clang-format on

        args.merge_patch(args_in);
        iterations = 10000;
        timestep = args["timestep_size"].get<double>();

        friction_coeff
            = args["rigid_body_problem"]["coefficient_friction"].get<double>();
        restitution_coeff = std::max(
            args["rigid_body_problem"]["coefficient_restitution"].get<double>(),
            0.0);

        array3d gravity_array
            = args["rigid_body_problem"]["gravity"].get<array3d>();
        m_world->SetGravity(b2Vec2(gravity_array[0], gravity_array[1]));

        LoadRigidBodies(args["rigid_body_problem"]["rigid_bodies"]);
    }

public:
    JSONExample(std::string filename)
    {
        std::ifstream input(filename);
        if (input.good()) {
            nlohmann::json scene = nlohmann::json::parse(input, nullptr, false);
            if (scene.is_discarded()) {
                std::cerr << "Invalid Json file" << std::endl;
            }
            if (scene.find("args") != scene.end()) {
                if (scene["args"] != nullptr) {
                    LoadScene(scene["args"]);
                }
            } else {
                LoadScene(scene);
            }
        } else {
            std::cerr << "Invalid Sim file" << std::endl;
        }
        SaveEdges();
        SaveStep();
    }

    static Test* Create() { return new JSONExample(INPUT_FILENAME); }

    virtual ~JSONExample()
    {
        nlohmann::json results;
        results["args"] = args;
        results["active_args"] = nlohmann::json();
        results["animation"] = nlohmann::json();
        results["animation"]["vertices_sequence"] = vertices_sequence;
        results["animation"]["state_sequence"] = state_sequence;
        results["animation"]["edges"] = edges;

        std::ofstream o("../../../results/comparisons/Box2D/results.json");
        o << std::setw(4) << results << std::endl;
    }

    virtual void Step(Settings* settings) override
    {
        static bool firstCall = true;
        static int counter = 0;
        if (firstCall) {
            settings->enableSleep = false;
            settings->hz = 1 / timestep;
            settings->velocityIterations = iterations;
            settings->positionIterations = iterations;
            settings->drawContactPoints = true;
            settings->drawContactImpulse = true;
            settings->pause = true;
            firstCall = false;
        }
        // if (counter >= 1000) {
        //     settings->pause = true;
        // }

        Test::Step(settings);

        if (!settings->pause || settings->singleStep) {
            SaveStep();
            counter++;
        }
    }
};
