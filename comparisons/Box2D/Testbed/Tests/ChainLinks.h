#ifndef CHAINLINKS_H
#define CHAINLINKS_H

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

class ChainLinks : public Test {
protected:
    float scale;
    float thickness;

public:
    static void createLinkVertexEdgeArrays(std::vector<b2Vec2>& vertices,
        std::vector<std::array<int, 2>>& edges,
        float thickness)
    {
        auto f1 = [](float x) -> float { return -1.25 * x + 1.625; };
        auto f2 = [](float x) -> float { return 1.25 * x + 0.375; };
        vertices = std::vector<b2Vec2>(10);
        vertices[0].Set(-1.5f * f1(thickness / 2), 0.0f);
        vertices[1].Set(0.0f, 0.0f);
        vertices[2].Set(1.5f * f1(thickness / 2), 0.0f);
        vertices[3].Set(-3.0f, -3.0f);
        vertices[4].Set(0.0f, -3.0f);
        vertices[5].Set(3.0f, -3.0f);
        vertices[6].Set(-3.0f, -6.0f);
        vertices[7].Set(-1.5f * f2(thickness / 2), -6.0f);
        vertices[8].Set(1.5f * f2(thickness / 2), -6.0f);
        vertices[9].Set(3.0f, -6.0f);

        edges = std::vector<std::array<int, 2>>(9);
        edges[0] = { 0, 1 };
        edges[1] = { 1, 2 };
        edges[2] = { 1, 4 };
        edges[3] = { 3, 4 };
        edges[4] = { 4, 5 };
        edges[5] = { 3, 6 };
        edges[6] = { 5, 9 };
        edges[7] = { 6, 7 };
        edges[8] = { 8, 9 };
    }

    virtual void createLink(b2World* world, float x, float y, b2BodyType type)
    {
        // Define the link.
        b2BodyDef linkBodyDef;
        linkBodyDef.type = type;
        linkBodyDef.bullet = type == b2_dynamicBody;
        linkBodyDef.position.Set(x, y);
        linkBodyDef.angle = 3.14159f / 2.0f;
        // Call the body factory which allocates memory for the ground body
        // from a pool and creates the link (also from a pool).
        // The body is also added to the world.
        b2Body* link = world->CreateBody(&linkBodyDef);

        // Define the link shape
        b2CircleShape circle;
        circle.m_radius = thickness / 2.0f;
        b2PolygonShape box;

        float density = 1.0f;

        std::vector<b2Vec2> vertices;
        std::vector<std::array<int, 2>> edges;
        createLinkVertexEdgeArrays(vertices, edges, thickness);

        // Add the shape fixtures to the link body
        float friction_coeff = 0.0;
        float restitution_coeff = 0.0;
        for (const auto& vertex : vertices) {
            circle.m_p = scale * vertex;
            b2Fixture* fixture = link->CreateFixture(&circle, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }
        for (const auto& edge : edges) {
            b2Vec2 v0 = scale * vertices[edge[0]];
            b2Vec2 v1 = scale * vertices[edge[1]];
            b2Vec2 center = 0.5 * (v1 + v0);
            float hx = (v1 - v0).Length() / 2;
            float hy = thickness / 2.0f;
            float angle = atan2((v1 - v0).y, (v1 - v0).x);
            box.SetAsBox(hx, hy, center, angle);
            b2Fixture* fixture = link->CreateFixture(&box, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }
    }

    ChainLinks(float thickness = 1.0, float scale = 1.0)
        : thickness(thickness)
        , scale(scale)
    {
        m_world->SetGravity(b2Vec2(0, -9.81));
        float y = 15;
        float x = -30;
        for (int i = 0; i < 30; i++) {
            createLink(m_world, x, y, i ? b2_dynamicBody : b2_staticBody);
            x += 4.5 * scale;
        }
    }

    static Test* Create() { return new ChainLinks(); }
    static Test* CreateThin() { return new ChainLinks(0.1); }

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

class SimpleChainLinks : public Test {
protected:
    float scale;

public:
    virtual void createLink(b2World* world, float x, float y, b2BodyType type)
    {
        // Define the link.
        b2BodyDef linkBodyDef;
        linkBodyDef.type = type;
        linkBodyDef.bullet = type == b2_dynamicBody;
        linkBodyDef.position.Set(x, y);
        linkBodyDef.angle = 3.14159f * 3.0f / 4.0f;
        // Call the body factory which allocates memory for the ground body
        // from a pool and creates the link (also from a pool).
        // The body is also added to the world.
        b2Body* link = world->CreateBody(&linkBodyDef);

        // Define the link shape
        b2EdgeShape edgeShape;

        float density = 1.0f;

        std::vector<b2Vec2> vertices;
        std::vector<std::array<int, 2>> edges;
        ChainLinks::createLinkVertexEdgeArrays(vertices, edges, 0.1);

        // Add the shape fixtures to the link body
        float friction_coeff = 0.0;
        float restitution_coeff = 0.0;
        for (const auto& edge : edges) {
            b2Vec2 v0 = scale * vertices[edge[0]];
            b2Vec2 v1 = scale * vertices[edge[1]];
            edgeShape.Set(v0, v1);
            b2Fixture* fixture = link->CreateFixture(&edgeShape, density);
            fixture->SetFriction(friction_coeff);
            fixture->SetRestitution(restitution_coeff);
        }
    }

    SimpleChainLinks(float scale = 1.0)
        : scale(scale)
    {
        m_world->SetGravity(b2Vec2(0, -9.81));
        float y = 15;
        float x = -30;
        for (int i = 0; i < 10; i++) {
            createLink(m_world, x, y, i ? b2_dynamicBody : b2_staticBody);
            x += 3 * scale;
            y += 3 * scale;
        }
    }

    static Test* Create() { return new SimpleChainLinks(); }

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

#endif
