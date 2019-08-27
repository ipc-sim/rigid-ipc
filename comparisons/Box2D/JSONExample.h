#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "json.hpp"

class JSONExample : public Test {
protected:
    float displacement_scale;
    std::vector<std::array<size_t, 2>> edges;
    std::vector<std::vector<std::array<double, 2>>> vertices_sequence;
    std::vector<nlohmann::json> state_sequence;

    double timestep;
    long iterations;

    void save_edges()
    {
        edges.clear();
        size_t starting_vertex_index = 0;
        const b2Body* body = m_world->GetBodyList();
        while (body != nullptr) {
            const b2Fixture* fixture = body->GetFixtureList();
            while (fixture != nullptr) {
                const b2Shape* shape = fixture->GetShape();
                // assert(shape->GetType() == Type::e_polygon);
                const b2PolygonShape* polygon
                    = dynamic_cast<const b2PolygonShape*>(shape);
                for (size_t i = 0; i < polygon->m_count; i++) {
                    edges.push_back({ { starting_vertex_index + i,
                        starting_vertex_index + (i + 1) % polygon->m_count } });
                }
                starting_vertex_index += polygon->m_count;
                fixture = fixture->GetNext();
            }
            body = body->GetNext();
        }
    }

    b2Vec2 b2Mat22_times_b2Vec2(b2Mat22 M, b2Vec2 v)
    {
        return b2Vec2(M.ex.x * v.x + M.ey.x * v.y, M.ex.y * v.x + M.ey.y * v.y);
    }

    double b2Vec2_dot_b2Vec2(b2Vec2 u, b2Vec2 v)
    {
        return u.x * v.x + u.y * v.y;
    }

    void save_step()
    {
        std::vector<std::array<double, 2>> vertices;
        double angular_momentum = 0;
        std::array<double, 2> linear_momentum = { { 0, 0 } };
        double kinetic_energy = 0;
        double potential_energy = 0;
        std::vector<nlohmann::json> body_states;

        const b2Body* body = m_world->GetBodyList();
        while (body != nullptr) {
            double theta = body->GetAngle();
            b2Mat22 R(cos(theta), -sin(theta), sin(theta), cos(theta));
            // Save the vertices
            const b2Fixture* fixture = body->GetFixtureList();
            while (fixture != nullptr) {
                const b2Shape* shape = fixture->GetShape();
                // assert(shape->GetType() == Type::e_polygon);
                const b2PolygonShape* polygon
                    = dynamic_cast<const b2PolygonShape*>(shape);
                for (size_t i = 0; i < polygon->m_count; i++) {
                    b2Vec2 vertex
                        = b2Mat22_times_b2Vec2(R, polygon->m_vertices[i])
                        + body->GetPosition();
                    vertices.push_back({ { vertex.x, vertex.y } });
                }
                fixture = fixture->GetNext();
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
                * b2Vec2_dot_b2Vec2(m_world->GetGravity(), body->GetPosition());

            nlohmann::json body_state;
            body_state["position"] = { body->GetPosition().x,
                body->GetPosition().y, body->GetAngle() };
            body_state["velocity"] = { body->GetLinearVelocity().x,
                body->GetLinearVelocity().y, body->GetAngularVelocity() };
            body_states.push_back(body_state);

            body = body->GetNext();
        }

        nlohmann::json state;
        state["angular_momentum"] = angular_momentum;
        state["linear_momentum"] = linear_momentum;
        state["kinetic_energy"] = kinetic_energy;
        state["potential_energy"] = potential_energy;
        state["rigid_bodies"] = body_states;
        vertices_sequence.push_back(vertices);
        state_sequence.push_back(state);
    }

public:
    JSONExample(std::string filename)
    {
        iterations = 300;
        timestep = 1e-2;

        m_world->SetGravity(b2Vec2(0, -9.8));

        float friction_coeff = 0.0;
        float restitution_coeff = 1.0;
        float density = 1.0f;

        {
            // Define the ground
            b2BodyDef groundBodyDef;
            groundBodyDef.type = b2_kinematicBody;
            groundBodyDef.position.Set(-10, 30);
            groundBodyDef.angularVelocity = 1;

            // Create the ground
            b2Body* ground = m_world->CreateBody(&groundBodyDef);

            // Define the ground shape
            b2PolygonShape groundShape;

            // Add the shape fixtures to the line body
            groundShape.SetAsBox(20, 1);
            b2Fixture* fixture = ground->CreateFixture(&groundShape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }

        {
            // Define the ground
            b2BodyDef boxBodyDef;
            boxBodyDef.type = b2_dynamicBody;
            boxBodyDef.position.Set(-10, 50);

            // Create the ground
            b2Body* box = m_world->CreateBody(&boxBodyDef);

            // Define the ground shape
            b2PolygonShape boxShape;

            // Add the shape fixtures to the line body
            boxShape.SetAsBox(1, 1);
            b2Fixture* fixture = box->CreateFixture(&boxShape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }

        save_edges();
        save_step();
    }

    static Test* Create() { return new JSONExample("example.json"); }

    virtual ~JSONExample()
    {
        nlohmann::json results;
        results["args"] = nlohmann::json();
        results["active_args"] = nlohmann::json();
        results["animation"] = nlohmann::json();
        results["animation"]["vertices_sequence"] = vertices_sequence;
        results["animation"]["state_sequence"] = state_sequence;
        results["animation"]["edges"] = edges;

        std::ofstream o("../../../results/simulations/Box2D/results.json");
        o << std::setw(4) << results << std::endl;
    }

    virtual void Step(Settings* settings) override
    {
        static bool firstCall = true;
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
        Test::Step(settings);

        if (!settings->pause || settings->singleStep) {
            save_step();
        }
    }
};
