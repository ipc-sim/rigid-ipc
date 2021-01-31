#include "rigid_body_collision_constraint.hpp"

namespace ipc::rigid {

RigidBodyVertexVertexConstraint::RigidBodyVertexVertexConstraint(
    const RigidBodyAssembler& bodies,
    long vertex0_index,
    long vertex1_index,
    long multiplicity)
    : multiplicity(multiplicity)
{
    bodies.global_to_local_vertex(
        vertex0_index, vertex0_body_id, vertex0_local_id);
    bodies.global_to_local_vertex(
        vertex1_index, vertex1_body_id, vertex1_local_id);
}

RigidBodyEdgeVertexConstraint::RigidBodyEdgeVertexConstraint(
    const RigidBodyAssembler& bodies,
    long edge_index,
    long vertex_index,
    long multiplicity)
    : multiplicity(multiplicity)
{
    const auto& edge = bodies.m_edges.row(edge_index);
    bodies.global_to_local_vertex(
        vertex_index, vertex_body_id, vertex_local_id);
    bodies.global_to_local_vertex(edge(0), edge_body_id, edge_vertex0_local_id);
    bodies.global_to_local_vertex(edge(1), edge_body_id, edge_vertex1_local_id);
}

RigidBodyEdgeEdgeConstraint::RigidBodyEdgeEdgeConstraint(
    const RigidBodyAssembler& bodies,
    long edge0_index,
    long edge1_index,
    double eps_x)
    : eps_x(eps_x)
{
    const auto& edge0 = bodies.m_edges.row(edge0_index);
    const auto& edge1 = bodies.m_edges.row(edge1_index);
    bodies.global_to_local_vertex(
        edge0(0), edge0_body_id, edge0_vertex0_local_id);
    bodies.global_to_local_vertex(
        edge0(1), edge0_body_id, edge0_vertex1_local_id);
    bodies.global_to_local_vertex(
        edge1(0), edge1_body_id, edge1_vertex0_local_id);
    bodies.global_to_local_vertex(
        edge1(1), edge1_body_id, edge1_vertex1_local_id);
}

RigidBodyFaceVertexConstraint::RigidBodyFaceVertexConstraint(
    const RigidBodyAssembler& bodies,
    long face_index,
    long vertex_index)
{
    const auto& face = bodies.m_faces.row(face_index);
    bodies.global_to_local_vertex(
        vertex_index, vertex_body_id, vertex_local_id);
    bodies.global_to_local_vertex(face(0), face_body_id, face_vertex0_local_id);
    bodies.global_to_local_vertex(face(1), face_body_id, face_vertex1_local_id);
    bodies.global_to_local_vertex(face(2), face_body_id, face_vertex2_local_id);
}

} // namespace ipc::rigid
