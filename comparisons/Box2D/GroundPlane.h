#pragma once

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

class GroundPlane : public Test {
protected:
    float displacement_scale;

public:
    GroundPlane(float displacement_scale = 100000.0)
        : displacement_scale(displacement_scale)
    {
        m_world->SetGravity(b2Vec2(0, 0));

        float friction_coeff = 0.0;
        float restitution_coeff = 0.0;
        float density = 1.0f;

        {
            // Define the ground
            b2BodyDef groundBodyDef;
            groundBodyDef.type = b2_staticBody;
            groundBodyDef.position.Set(-10, 30);

            // Create the ground
            b2Body* ground = m_world->CreateBody(&groundBodyDef);

            // Define the ground shape
            b2PolygonShape groundShape;

            // Add the shape fixtures to the line body
            groundShape.SetAsBox(20, 0.01);
            b2Fixture* fixture = ground->CreateFixture(&groundShape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }

        {
            // Define the line
            b2BodyDef lineBodyDef;
            lineBodyDef.type = b2_dynamicBody;
            lineBodyDef.bullet = true;
            lineBodyDef.position.Set(0.0f, 35.0f);
            lineBodyDef.linearVelocity = b2Vec2(0, -displacement_scale);

            b2Body* line = m_world->CreateBody(&lineBodyDef);

            // Define the line shape
            b2EdgeShape lineShape;

            // Add the shape fixtures to the line body
            lineShape.Set(b2Vec2(-15, 0), b2Vec2(-5, 0));
            b2Fixture* fixture = line->CreateFixture(&lineShape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }
    }

    static Test* Create() { return new GroundPlane(); }

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
