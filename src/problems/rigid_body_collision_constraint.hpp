#pragma once

#include <array>

#include <ipc/collision_constraint.hpp>
#include <ipc/friction/friction_constraint.hpp>

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
            long vertex0_index,
            long vertex1_index,
            long multiplicity);
        RigidBodyVertexVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::VertexVertexConstraint& vv_constraint)
            : RigidBodyVertexVertexConstraint(
                  bodies,
                  vv_constraint.vertex0_index,
                  vv_constraint.vertex1_index,
                  vv_constraint.multiplicity)
        {
        }
        RigidBodyVertexVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::VertexVertexFrictionConstraint& vv_constraint)
            : RigidBodyVertexVertexConstraint(
                  bodies,
                  vv_constraint.vertex0_index,
                  vv_constraint.vertex1_index,
                  vv_constraint.multiplicity)
        {
        }

        std::array<long, 2> body_ids() const
        {
            return { { vertex0_body_id, vertex1_body_id } };
        }

        static std::vector<uint8_t> vertex_local_body_ids()
        {
            return { { 0, 1 } };
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
            long edge_index,
            long vertex_index,
            long multiplicity);
        RigidBodyEdgeVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::EdgeVertexConstraint& ev_constraint)
            : RigidBodyEdgeVertexConstraint(
                  bodies,
                  ev_constraint.edge_index,
                  ev_constraint.vertex_index,
                  ev_constraint.multiplicity)
        {
        }
        RigidBodyEdgeVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::EdgeVertexFrictionConstraint& ev_constraint)
            : RigidBodyEdgeVertexConstraint(
                  bodies,
                  ev_constraint.edge_index,
                  ev_constraint.vertex_index,
                  ev_constraint.multiplicity)
        {
        }

        std::array<long, 2> body_ids() const
        {
            return { { vertex_body_id, edge_body_id } };
        }

        static std::vector<uint8_t> vertex_local_body_ids()
        {
            return { { 0, 1, 1 } };
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
            long edge0_index,
            long edge1_index,
            double eps_x);
        RigidBodyEdgeEdgeConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::EdgeEdgeConstraint& ee_constraint)
            : RigidBodyEdgeEdgeConstraint(
                  bodies,
                  ee_constraint.edge0_index,
                  ee_constraint.edge1_index,
                  ee_constraint.eps_x)
        {
        }
        RigidBodyEdgeEdgeConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::EdgeEdgeFrictionConstraint& ee_constraint)
            : RigidBodyEdgeEdgeConstraint(
                  bodies,
                  ee_constraint.edge0_index,
                  ee_constraint.edge1_index,
                  -1)
        {
        }

        std::array<long, 2> body_ids() const
        {
            return { { edge0_body_id, edge1_body_id } };
        }

        static std::vector<uint8_t> vertex_local_body_ids()
        {
            return { { 0, 0, 1, 1 } };
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
            long face_index,
            long vertex_index);
        RigidBodyFaceVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::FaceVertexConstraint& fv_constraint)
            : RigidBodyFaceVertexConstraint(
                  bodies, fv_constraint.face_index, fv_constraint.vertex_index)
        {
        }
        RigidBodyFaceVertexConstraint(
            const physics::RigidBodyAssembler& bodies,
            const ipc::FaceVertexFrictionConstraint& fv_constraint)
            : RigidBodyFaceVertexConstraint(
                  bodies, fv_constraint.face_index, fv_constraint.vertex_index)
        {
        }

        std::array<long, 2> body_ids() const
        {
            return { { vertex_body_id, face_body_id } };
        }

        static std::vector<uint8_t> vertex_local_body_ids()
        {
            return { { 0, 1, 1, 1 } };
        }
    };

    template <typename Constraints>
    std::vector<uint8_t>
    vertex_local_body_ids(const Constraints& constraints, size_t ci)
    {
        if (ci < constraints.vv_constraints.size()) {
            return RigidBodyVertexVertexConstraint::vertex_local_body_ids();
        }
        ci -= constraints.vv_constraints.size();
        if (ci < constraints.ev_constraints.size()) {
            return RigidBodyEdgeVertexConstraint::vertex_local_body_ids();
        }
        ci -= constraints.ev_constraints.size();
        if (ci < constraints.ee_constraints.size()) {
            return RigidBodyEdgeEdgeConstraint::vertex_local_body_ids();
        }
        ci -= constraints.ee_constraints.size();
        if (ci < constraints.fv_constraints.size()) {
            return RigidBodyFaceVertexConstraint::vertex_local_body_ids();
        }
        assert(false);
        throw "Invalid constraint index!";
    }

    template <typename Constraints>
    std::array<long, 2> body_ids(
        const physics::RigidBodyAssembler& bodies,
        const Constraints& constraints,
        size_t ci)
    {
        if (ci < constraints.vv_constraints.size()) {
            return RigidBodyVertexVertexConstraint(
                       bodies, constraints.vv_constraints[ci])
                .body_ids();
        }
        ci -= constraints.vv_constraints.size();
        if (ci < constraints.ev_constraints.size()) {
            return RigidBodyEdgeVertexConstraint(
                       bodies, constraints.ev_constraints[ci])
                .body_ids();
        }
        ci -= constraints.ev_constraints.size();
        if (ci < constraints.ee_constraints.size()) {
            return RigidBodyEdgeEdgeConstraint(
                       bodies, constraints.ee_constraints[ci])
                .body_ids();
        }
        ci -= constraints.ee_constraints.size();
        if (ci < constraints.fv_constraints.size()) {
            return RigidBodyFaceVertexConstraint(
                       bodies, constraints.fv_constraints[ci])
                .body_ids();
        }
        assert(false);
        throw "Invalid constraint index!";
    }

} // namespace opt
} // namespace ccd
