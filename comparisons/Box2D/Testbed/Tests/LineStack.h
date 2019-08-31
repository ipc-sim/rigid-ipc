#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

class LineStack : public Test {
protected:
    int num_lines;
    float displacement_scale;

public:
    LineStack(int num_lines = 20, float displacement_scale = 1000.0)
        : num_lines(num_lines)
        , displacement_scale(displacement_scale)
    {
        m_world->SetGravity(b2Vec2(0, 0));

        // Define the lines.
        b2BodyDef lineBodyDef;
        lineBodyDef.type = b2_dynamicBody;
        lineBodyDef.bullet = true;

        float friction_coeff = 0.0;
        float restitution_coeff = 0.0;
        float density = 1.0f;
        for (int i = 0; i < num_lines; i++) {
            lineBodyDef.position.Set(-10, 30 * (1.0f - i / float(num_lines)));

            // Call the body factory which allocates memory for the ground body
            // from a pool and creates the line (also from a pool).
            // The body is also added to the world.
            b2Body* line = m_world->CreateBody(&lineBodyDef);

            // Define the line shape
            b2PolygonShape lineShape;

            // Add the shape fixtures to the line body
            lineShape.SetAsBox(20, 0.01);
            b2Fixture* fixture = line->CreateFixture(&lineShape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }

        {
            // lineBodyDef.angle = 3.1415926f * 0.1;
            lineBodyDef.position.Set(-10.0f, 35.0f);
            lineBodyDef.linearVelocity = b2Vec2(0, -1 * displacement_scale);

            // Call the body factory which allocates memory for the ground body
            // from a pool and creates the line (also from a pool).
            // The body is also added to the world.
            b2Body* line = m_world->CreateBody(&lineBodyDef);

            // Define the line shape
            b2EdgeShape lineShape;

            // Add the shape fixtures to the line body
            lineShape.Set(b2Vec2(-5, 0), b2Vec2(5, 0));
            b2Fixture* fixture = line->CreateFixture(&lineShape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }
    }

    static Test* Create() { return new LineStack(); }

    virtual void Step(Settings* settings) override
    {
        static bool firstCall = true;
        if (firstCall) {
            settings->enableSleep = false;
            settings->hz = 120; // 1 / 1e-3;
            settings->velocityIterations = 300;
            settings->positionIterations = 300;
            settings->drawContactPoints = true;
            settings->drawContactImpulse = true;
            firstCall = false;
        }
        Test::Step(settings);
    }
};
