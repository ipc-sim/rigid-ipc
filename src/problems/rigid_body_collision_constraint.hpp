#pragma once

#include <array>

#include <ipc/collision_constraint.hpp>

#include <physics/rigid_body_assembler.hpp>

namespace ccd {
namespace opt {

    struct RigidBodyVertexVertexConstraint {
        long vertex0_body_id;
        long vertex1_body_id;
        long vertex0_local_id;
        long vertex1_local_id;

        long multiplicity;

        RigidBodyVertexVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::VertexVertexConstraint& vv_constraint);

        inline std::array<long, 2> body_ids()
        {
            return { { vertex0_body_id, vertex1_body_id } };
        }
    };

    struct RigidBodyEdgeVertexConstraint {
        long vertex_body_id;
        long edge_body_id;
        long vertex_local_id;
        long edge_vertex0_local_id;
        long edge_vertex1_local_id;

        long multiplicity;

        RigidBodyEdgeVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::EdgeVertexConstraint& ev_constraint);

        inline std::array<long, 2> body_ids()
        {
            return { { vertex_body_id, edge_body_id } };
        }
    };

    struct RigidBodyEdgeEdgeConstraint {
        long edge0_body_id;
        long edge1_body_id;
        long edge0_vertex0_local_id;
        long edge0_vertex1_local_id;
        long edge1_vertex0_local_id;
        long edge1_vertex1_local_id;

        const long multiplicity = 1;

        double eps_x;

        RigidBodyEdgeEdgeConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::EdgeEdgeConstraint& ee_constraint);

        inline std::array<long, 2> body_ids()
        {
            return { { edge0_body_id, edge1_body_id } };
        }
    };

    struct RigidBodyFaceVertexConstraint {
        long vertex_body_id;
        long face_body_id;
        long vertex_local_id;
        long face_vertex0_local_id;
        long face_vertex1_local_id;
        long face_vertex2_local_id;

        const long multiplicity = 1;

        RigidBodyFaceVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::FaceVertexConstraint& fv_constraint);

        inline std::array<long, 2> body_ids()
        {
            return { { vertex_body_id, face_body_id } };
        }
    };

} // namespace opt
} // namespace ccd
