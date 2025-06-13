#include "rigid_body_bvh.hpp"

namespace ipc::rigid {

void detect_body_pair_collision_candidates_from_aabbs(
    const RigidBodyAssembler& bodies,
    const std::vector<AABB>& bodyA_vertex_aabbs,
    const int bodyA_id,
    const int bodyB_id,
    const int collision_types,
    Candidates& candidates,
    const double inflation_radius)
{
    bool build_ev = collision_types & CollisionType::EDGE_VERTEX;
    bool build_ee = collision_types & CollisionType::EDGE_EDGE;
    bool build_fv = collision_types & CollisionType::FACE_VERTEX;
    auto add_ev = [&](size_t eai, size_t vbi) {
        if (build_ev) {
            candidates.ev_candidates.emplace_back(
                bodies.m_body_edge_id[bodyA_id] + eai,
                bodies.m_body_vertex_id[bodyB_id] + vbi);
        }
    };
    auto add_ve = [&](size_t vai, size_t ebi) {
        if (build_ev) {
            candidates.ev_candidates.emplace_back(
                bodies.m_body_edge_id[bodyB_id] + ebi,
                bodies.m_body_vertex_id[bodyA_id] + vai);
        }
    };
    auto add_ee = [&](size_t eai, size_t ebi) {
        if (build_ee) {
            candidates.ee_candidates.emplace_back(
                bodies.m_body_edge_id[bodyA_id] + eai,
                bodies.m_body_edge_id[bodyB_id] + ebi);
        }
    };
    auto add_fv = [&](size_t fai, size_t vbi) {
        if (build_fv) {
            candidates.fv_candidates.emplace_back(
                bodies.m_body_face_id[bodyA_id] + fai,
                bodies.m_body_vertex_id[bodyB_id] + vbi);
        }
    };
    auto add_vf = [&](size_t vai, size_t fbi) {
        if (build_fv) {
            candidates.fv_candidates.emplace_back(
                bodies.m_body_face_id[bodyB_id] + fbi,
                bodies.m_body_vertex_id[bodyA_id] + vai);
        }
    };

    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];

    const std::vector<AABB> bodyB_vertex_aabbs =
        vertex_aabbs(bodyB.vertices, inflation_radius);
    const Eigen::MatrixXi &EA = bodyA.edges, &EB = bodyB.edges,
                          &FA = bodyA.faces, &FB = bodyB.faces;

    const auto& selectorA = bodyA.mesh_selector;
    const auto& selectorB = bodyB.mesh_selector;

    auto bodyA_edge_aabb = [&](size_t ei) {
        return AABB(
            bodyA_vertex_aabbs[EA(ei, 0)], bodyA_vertex_aabbs[EA(ei, 1)]);
    };
    auto bodyB_edge_aabb = [&](size_t ei) {
        return AABB(
            bodyB_vertex_aabbs[EB(ei, 0)], bodyB_vertex_aabbs[EB(ei, 1)]);
    };
    auto bodyA_face_aabb = [&](size_t fi) {
        return AABB(
            bodyA_vertex_aabbs[FA(fi, 0)], bodyA_vertex_aabbs[FA(fi, 1)],
            bodyA_vertex_aabbs[FA(fi, 2)]);
    };
    auto bodyB_face_aabb = [&](size_t fi) {
        return AABB(
            bodyB_vertex_aabbs[FB(fi, 0)], bodyB_vertex_aabbs[FB(fi, 1)],
            bodyB_vertex_aabbs[FB(fi, 2)]);
    };

    ///////////////////////////////////////////////////////////////////////////
    // query (f, *)
    for (size_t fa_id = 0; fa_id < FA.rows(); fa_id++) {
        // all (f_e, *v) and (f_v, *e) are not needed because faces are only 3D
        assert(!build_ev);

        // Construct a bbox of bodyA's face
        AABB fa_aabb = bodyA_face_aabb(fa_id);

        std::vector<unsigned int> ids;
        bodyB.bvh.intersect_box(
            // Grow the box by inflation_radius because the BVH is not grown
            fa_aabb.getMin().array() - inflation_radius,
            fa_aabb.getMax().array() + inflation_radius, //
            ids);

        for (const auto& id : ids) {
            if (id < bodyB.num_codim_vertices()) {

                // (f, cv) - no need to do a AABB check
                add_fv(fa_id, selectorB.codim_vertices_to_vertices(id));
                // ignore (f_e, cv) and (f_v, cv)

            } else if (
                id < bodyB.num_codim_vertices() + bodyB.num_codim_edges()) {

                size_t eb_id = selectorB.codim_edges_to_edges(
                    id - bodyB.num_codim_vertices());

                // (f, ce_v)
                for (int vi = 0; vi < EB.cols(); vi++) {
                    size_t vb_id = EB(eb_id, vi);
                    if (selectorB.vertex_to_edge(vb_id) == eb_id
                        && AABB::are_overlapping(
                               fa_aabb, bodyB_vertex_aabbs[vb_id])) {
                        add_fv(fa_id, vb_id);
                    }
                }

                // (f_e, ce)
                AABB eb_aabb = bodyB_edge_aabb(eb_id);
                for (int ei = 0; ei < FA.cols(); ei++) {
                    size_t ea_id = selectorA.face_to_edge(fa_id, ei);
                    if (selectorA.edge_to_face(ea_id) == fa_id) {
                        AABB ea_aabb = bodyA_edge_aabb(ea_id);
                        if (AABB::are_overlapping(ea_aabb, eb_aabb)) {
                            add_ee(ea_id, eb_id);
                        }
                    }
                }

                // ignore (f, ce), (f_v, ce), (f_v, ce_v), and (f_e, ce_v)

            } else {

                size_t fb_id =
                    id - bodyB.num_codim_vertices() - bodyB.num_codim_edges();

                AABB fb_aabb = bodyB_face_aabb(fb_id);
                for (int f_vi = 0; f_vi < FA.cols(); f_vi++) {
                    // (f_v, f)
                    long va_id = FA(fa_id, f_vi);
                    if (selectorA.vertex_to_face(va_id) == fa_id) {
                        if (AABB::are_overlapping(
                                bodyA_vertex_aabbs[va_id], fb_aabb)) {
                            // Convert the local ids to the global ones
                            add_vf(va_id, fb_id);
                        }
                    }

                    // (f, f_v)
                    long vb_id = FB(fb_id, f_vi);
                    if (selectorB.vertex_to_face(vb_id) == fb_id) {
                        if (AABB::are_overlapping(
                                fa_aabb, bodyB_vertex_aabbs[vb_id])) {
                            // Convert the local ids to the global ones
                            add_fv(fa_id, vb_id);
                        }
                    }
                }

                for (int fa_ei = 0; fa_ei < FA.cols(); fa_ei++) {
                    long ea_id = selectorA.face_to_edge(fa_id, fa_ei);

                    if (selectorA.edge_to_face(ea_id) != fa_id) {
                        continue;
                    }

                    AABB ea_aabb = bodyA_edge_aabb(ea_id);

                    for (int fb_ei = 0; fb_ei < FB.cols(); fb_ei++) {
                        long eb_id = selectorB.face_to_edge(fb_id, fb_ei);

                        if (selectorB.edge_to_face(eb_id) != fb_id) {
                            continue;
                        }

                        AABB eb_aabb = bodyB_edge_aabb(eb_id);

                        if (AABB::are_overlapping(ea_aabb, eb_aabb)) {
                            // Convert the local ids to the global ones
                            add_ee(ea_id, eb_id);
                        }
                    }
                }

                // ignore (f, f), (f, f_e), (f_v, f_v), (f_v, f_e), (f_e, f),
                // (f_e, f_v)
            }
        }
    }

    // query (ce, *)
    for (long co_ea_id = 0; co_ea_id < bodyA.num_codim_edges(); co_ea_id++) {
        size_t ea_id = selectorA.codim_edges_to_edges(co_ea_id);

        AABB ea_aabb = bodyA_edge_aabb(ea_id);

        std::vector<unsigned int> ids;
        bodyB.bvh.intersect_box(
            // Grow the box by inflation_radius because the BVH is not grown
            ea_aabb.getMin().array() - inflation_radius,
            ea_aabb.getMax().array() + inflation_radius, //
            ids);

        for (const auto& id : ids) {
            if (id < bodyB.num_codim_vertices()) {
                size_t vb_id = selectorB.codim_vertices_to_vertices(id);

                // (ce, cv)
                add_ev(ea_id, vb_id);

            } else if (
                id < bodyB.num_codim_edges() + bodyB.num_codim_vertices()) {
                size_t eb_id = selectorB.codim_edges_to_edges(id);

                // (ce, ce)
                add_ee(ea_id, eb_id);

                for (int vi = 0; vi < EB.cols(); vi++) {
                    // (ce, ce_v)
                    size_t vb_id = EB(eb_id, vi);
                    if (selectorB.vertex_to_edge(vb_id) == eb_id) {
                        add_ev(ea_id, vb_id);
                    }

                    // (ce_v, ce)
                    size_t va_id = EA(ea_id, vi);
                    if (selectorA.vertex_to_edge(va_id) == ea_id) {
                        add_ve(va_id, eb_id);
                    }
                }

                // (ce_v, ce_v) is not needed

            } else {
                // (ce, f*)
                size_t fb_id =
                    id - bodyB.num_codim_vertices() - bodyB.num_codim_edges();

                // (ce_v, f_v) is not needed
                // (ce_v, f_e) is not needed because in 3D
                // (ce, f_v) is not needed because in 3D
                assert(!build_ev);

                // (ce_v, f)
                for (int vi = 0; vi < EA.cols(); vi++) {
                    size_t va_id = EA(ea_id, vi);
                    if (selectorA.vertex_to_edge(va_id) == ea_id) {
                        add_vf(va_id, fb_id);
                    }
                }

                // (ce, f_e)
                for (int ei = 0; ei < FB.cols(); ei++) {
                    size_t eb_id = selectorB.face_to_edge(fb_id, ei);
                    if (selectorB.edge_to_face(eb_id) == fb_id) {
                        add_ee(ea_id, eb_id);
                    }
                }

                // (ce, f) is not needed
            }
        }
    }

    // query (cv, *)
    for (long co_va_id = 0; co_va_id < bodyA.num_codim_vertices(); co_va_id++) {
        size_t va_id = selectorA.codim_vertices_to_vertices(co_va_id);

        AABB va_aabb = bodyA_vertex_aabbs[va_id];

        std::vector<unsigned int> ids;
        bodyB.bvh.intersect_box(
            // Grow the box by inflation_radius because the BVH is not grown
            va_aabb.getMin().array() - inflation_radius,
            va_aabb.getMax().array() + inflation_radius, //
            ids);

        for (const auto& id : ids) {
            if (id < bodyB.num_codim_vertices()) {
                // (cv, cv) is not needed
            } else if (
                id < bodyB.num_codim_vertices() + bodyB.num_codim_edges()) {

                size_t eb_id = selectorB.codim_edges_to_edges(
                    id - bodyB.num_codim_vertices());

                // (cv, ce)
                add_ev(eb_id, va_id);

                // (cv, ce_v) is not needed
            } else {
                // (cv, f)
                size_t fb_id =
                    id - bodyB.num_codim_vertices() - bodyB.num_codim_edges();
                add_vf(va_id, fb_id);

                // (cv, f_e) is not needed because in 3D
                assert(!build_ev);

                // (cv, f_v) is not needed
            }
        }
        // for (int fb_id = 0; fb_id < FB.rows(); fb_id++) {
        //     if (AABB::are_overlapping(va_aabb, bodyB_face_aabb(fb_id))) {
        //         add_vf(va_id, fb_id);
        //     }
        // }
    }
}

void detect_body_pair_intersection_candidates_from_aabbs(
    const RigidBodyAssembler& bodies,
    const std::vector<AABB>& bodyA_vertex_aabbs,
    const int bodyA_id,
    const int bodyB_id,
    std::vector<EdgeFaceCandidate>& ef_candidates,
    const double inflation_radius)
{
    auto add_ef = [&](size_t eai, size_t fbi) {
        ef_candidates.emplace_back(
            bodies.m_body_edge_id[bodyA_id] + eai,
            bodies.m_body_face_id[bodyB_id] + fbi);
    };
    auto add_fe = [&](size_t fai, size_t ebi) {
        ef_candidates.emplace_back(
            bodies.m_body_edge_id[bodyB_id] + ebi,
            bodies.m_body_face_id[bodyA_id] + fai);
    };

    const RigidBody& bodyA = bodies[bodyA_id];
    const RigidBody& bodyB = bodies[bodyB_id];

    const std::vector<AABB> bodyB_vertex_aabbs =
        vertex_aabbs(bodyB.vertices, inflation_radius);
    const Eigen::MatrixXi &EA = bodyA.edges, &EB = bodyB.edges,
                          &FA = bodyA.faces, &FB = bodyB.faces;

    const auto& selectorA = bodyA.mesh_selector;
    const auto& selectorB = bodyB.mesh_selector;

    auto bodyA_edge_aabb = [&](size_t ei) {
        return AABB(
            bodyA_vertex_aabbs[EA(ei, 0)], bodyA_vertex_aabbs[EA(ei, 1)]);
    };
    auto bodyB_edge_aabb = [&](size_t ei) {
        return AABB(
            bodyB_vertex_aabbs[EB(ei, 0)], bodyB_vertex_aabbs[EB(ei, 1)]);
    };
    auto bodyA_face_aabb = [&](size_t fi) {
        return AABB(
            bodyA_vertex_aabbs[FA(fi, 0)], bodyA_vertex_aabbs[FA(fi, 1)],
            bodyA_vertex_aabbs[FA(fi, 2)]);
    };
    auto bodyB_face_aabb = [&](size_t fi) {
        return AABB(
            bodyB_vertex_aabbs[FB(fi, 0)], bodyB_vertex_aabbs[FB(fi, 1)],
            bodyB_vertex_aabbs[FB(fi, 2)]);
    };

    ///////////////////////////////////////////////////////////////////////////
    // query (f, *)
    for (size_t fa_id = 0; fa_id < FA.rows(); fa_id++) {
        // Construct a bbox of bodyA's face
        AABB fa_aabb = bodyA_face_aabb(fa_id);

        std::vector<unsigned int> ids;
        bodyB.bvh.intersect_box(
            // Grow the box by inflation_radius because the BVH is not grown
            fa_aabb.getMin().array() - inflation_radius,
            fa_aabb.getMax().array() + inflation_radius, //
            ids);

        for (const auto& id : ids) {
            if (id < bodyB.num_codim_vertices()) {
                // ignore (f, cv)
            } else if (
                id < bodyB.num_codim_edges() + bodyB.num_codim_vertices()) {

                // (f, ce)
                size_t eb_id = selectorB.codim_edges_to_edges(
                    id - bodyB.num_codim_vertices());
                add_fe(fa_id, eb_id);

            } else {
                size_t fb_id =
                    id - bodyB.num_codim_vertices() - bodyB.num_codim_edges();

                AABB fb_aabb = bodyB_face_aabb(fb_id);

                for (int ei = 0; ei < FA.cols(); ei++) {
                    long ea_id = bodyA.mesh_selector.face_to_edge(fa_id, ei);
                    if (selectorA.edge_to_face(ea_id) == fa_id) {
                        AABB ea_aabb = bodyA_edge_aabb(ea_id);
                        if (AABB::are_overlapping(ea_aabb, fb_aabb)) {
                            add_ef(ea_id, fb_id);
                        }
                    }

                    long eb_id = bodyB.mesh_selector.face_to_edge(fb_id, ei);
                    if (bodyB.mesh_selector.edge_to_face(eb_id) == fb_id) {
                        AABB eb_aabb = bodyB_edge_aabb(eb_id);
                        if (AABB::are_overlapping(fa_aabb, eb_aabb)) {
                            add_fe(fa_id, eb_id);
                        }
                    }
                }
            }
        }
    }

    // query (ce, *)
    for (long co_ea_id = 0; co_ea_id < bodyA.num_codim_edges(); co_ea_id++) {
        size_t ea_id = selectorA.codim_edges_to_edges(co_ea_id);

        AABB ea_aabb = bodyA_edge_aabb(ea_id);

        std::vector<unsigned int> ids;
        bodyB.bvh.intersect_box(
            // Grow the box by inflation_radius because the BVH is not grown
            ea_aabb.getMin().array() - inflation_radius,
            ea_aabb.getMax().array() + inflation_radius, //
            ids);

        for (const auto& id : ids) {
            if (id < bodyB.num_codim_vertices()) {
                // ignore (ce, cv)
            } else if (
                id < bodyB.num_codim_edges() + bodyB.num_codim_vertices()) {
                // ignore (ce, cv)
            } else {
                size_t fb_id =
                    id - bodyB.num_codim_vertices() - bodyB.num_codim_edges();
                // (ce, f)
                add_ef(ea_id, fb_id);
            }
        }
    }

    // no need to query (cv, *)
}

} // namespace ipc::rigid
