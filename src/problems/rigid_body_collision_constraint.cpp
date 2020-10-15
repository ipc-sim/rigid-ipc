#include "rigid_body_collision_constraint.hpp"

namespace ccd {
namespace opt {

    RigidBodyVertexVertexConstraint::RigidBodyVertexVertexConstraint(
        const physics::RigidBodyAssembler& bodies,
        const ipc::VertexVertexConstraint& vv_constraint)
    {
        bodies.global_to_local_vertex(
            vv_constraint.vertex0_index, vertex0_body_id, vertex0_local_id);
        bodies.global_to_local_vertex(
            vv_constraint.vertex1_index, vertex1_body_id, vertex1_local_id);
        multiplicity = vv_constraint.multiplicity;
    }

    RigidBodyEdgeVertexConstraint::RigidBodyEdgeVertexConstraint(
        const physics::RigidBodyAssembler& bodies,
        const ipc::EdgeVertexConstraint& ev_constraint)
    {
        const auto& edge = bodies.m_edges.row(ev_constraint.edge_index);
        bodies.global_to_local_vertex(
            ev_constraint.vertex_index, vertex_body_id, vertex_local_id);
        bodies.global_to_local_vertex(
            edge(0), edge_body_id, edge_vertex0_local_id);
        bodies.global_to_local_vertex(
            edge(1), edge_body_id, edge_vertex1_local_id);
        multiplicity = ev_constraint.multiplicity;
    }

    RigidBodyEdgeEdgeConstraint::RigidBodyEdgeEdgeConstraint(
        const physics::RigidBodyAssembler& bodies,
        const ipc::EdgeEdgeConstraint& ee_constraint)
    {
        const auto& edge0 = bodies.m_edges.row(ee_constraint.edge0_index);
        const auto& edge1 = bodies.m_edges.row(ee_constraint.edge1_index);
        bodies.global_to_local_vertex(
            edge0(0), edge0_body_id, edge0_vertex0_local_id);
        bodies.global_to_local_vertex(
            edge0(1), edge0_body_id, edge0_vertex1_local_id);
        bodies.global_to_local_vertex(
            edge1(0), edge1_body_id, edge1_vertex0_local_id);
        bodies.global_to_local_vertex(
            edge1(1), edge1_body_id, edge1_vertex1_local_id);
        eps_x = ee_constraint.eps_x;
    }

    RigidBodyFaceVertexConstraint::RigidBodyFaceVertexConstraint(
        const physics::RigidBodyAssembler& bodies,
        const ipc::FaceVertexConstraint& fv_constraint)
    {
        const auto& face = bodies.m_faces.row(fv_constraint.face_index);
        bodies.global_to_local_vertex(
            fv_constraint.vertex_index, vertex_body_id, vertex_local_id);
        bodies.global_to_local_vertex(
            face(0), face_body_id, face_vertex0_local_id);
        bodies.global_to_local_vertex(
            face(1), face_body_id, face_vertex1_local_id);
        bodies.global_to_local_vertex(
            face(2), face_body_id, face_vertex2_local_id);
    }

} // namespace opt
} // namespace ccd
