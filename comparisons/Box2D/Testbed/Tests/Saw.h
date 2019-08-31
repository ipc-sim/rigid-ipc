#pragma once

#include <cmath>
#include <iostream>
#include <vector>

class Saw : public Test {
protected:
    int num_lines;
    float displacement_scale;

public:
    Saw()
    {
        m_world->SetGravity(b2Vec2(0, 0));

        float friction_coeff = 0.0;
        float restitution_coeff = 0.0;
        float density = 1.0f;

        // Define the box.
        {
            b2BodyDef boxBodyDef;
            boxBodyDef.type = b2_dynamicBody;
            boxBodyDef.bullet = true;
            boxBodyDef.position.Set(1, 2.5);
            b2Body* box = m_world->CreateBody(&boxBodyDef);
            box->SetLinearVelocity(b2Vec2(1, -10));
            b2PolygonShape boxShape;
            boxShape.SetAsBox(2, 0.5);
            b2Fixture* boxFixture = box->CreateFixture(&boxShape, density);
            boxFixture->SetFriction(friction_coeff);
            boxFixture->SetRestitution(restitution_coeff);
        }

        {
            b2BodyDef sawBodyDef;
            sawBodyDef.type = b2_staticBody;
            sawBodyDef.bullet = true;
            sawBodyDef.position.Set(0, 0);
            b2Body* saw = m_world->CreateBody(&sawBodyDef);
            b2ChainShape sawShape;
            std::vector<b2Vec2> sawPoints;
            const double delta_x = 1.0;
            for (double x = 0; x < 10; x += delta_x) {
                sawPoints.push_back(b2Vec2(x, 1));
                sawPoints.push_back(b2Vec2(x + delta_x, 0));
            }
            sawShape.CreateChain(&sawPoints.front(), sawPoints.size());
            b2Fixture* sawFixture = saw->CreateFixture(&sawShape, density);
            sawFixture->SetFriction(friction_coeff);
            sawFixture->SetRestitution(restitution_coeff);
        }
    }

    static Test* Create() { return new Saw(); }

    virtual void Step(Settings* settings) override
    {
        static bool firstCall = true;
        if (firstCall) {
            settings->enableSleep = false;
            settings->hz = 1 / 1e-1;
            settings->velocityIterations = 300;
            settings->positionIterations = 300;
            settings->drawContactPoints = true;
            settings->drawContactImpulse = true;
            firstCall = false;
        }
        Test::Step(settings);
    }
};
