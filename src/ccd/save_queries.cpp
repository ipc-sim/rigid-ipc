#include "save_queries.hpp"

#include <io/serialize_json.hpp>

namespace ipc::rigid {

struct CombinedQueryWriter {

    virtual ~CombinedQueryWriter() { save_queries(); }

    void append_query(const nlohmann::json& query)
    {
        std::scoped_lock lock(append_query_mutex);
        queries.push_back(query);
        if (queries.size() >= 1e4) {
            save_queries();
            queries.clear();
        }
    }

    void save_queries()
    {
        static int file_counter = 0;
        if (queries.size()) {
            nlohmann::json queries_json;
            queries_json["queries"] = queries;
            std::ofstream(
                fmt::format("ccd-queries-{:03d}.json", file_counter++))
                << queries_json.dump();
        }
    }

private:
    std::mutex append_query_mutex;
    std::vector<nlohmann::json> queries;
};

static CombinedQueryWriter combined_query_writer;

void save_ccd_candidate(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeVertexCandidate& ev_candidate)
{
    nlohmann::json query;
    long bodyA_id, vertex_id, bodyB_id, edge_id;
    bodies.global_to_local_vertex(
        ev_candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_edge(ev_candidate.edge_index, bodyB_id, edge_id);
    const auto& bodyA = bodies[bodyA_id];
    const auto& bodyB = bodies[bodyB_id];
    query["type"] = "ev";

    query["edge"] = nlohmann::json();
    query["edge"]["vertex0"] =
        to_json(bodyB.vertices.row(bodyB.edges(edge_id, 0)).transpose());
    query["edge"]["vertex1"] =
        to_json(bodyB.vertices.row(bodyB.edges(edge_id, 1)).transpose());
    query["edge"]["pose_t0"] = nlohmann::json();
    query["edge"]["pose_t0"]["position"] =
        to_json(poses_t0[bodyB_id].position);
    query["edge"]["pose_t0"]["rotation"] =
        to_json(poses_t0[bodyB_id].rotation);
    query["edge"]["pose_t1"] = nlohmann::json();
    query["edge"]["pose_t1"]["position"] =
        to_json(poses_t1[bodyB_id].position);
    query["edge"]["pose_t1"]["rotation"] =
        to_json(poses_t1[bodyB_id].rotation);

    query["vertex"] = nlohmann::json();
    query["vertex"]["vertex"] =
        to_json(bodyA.vertices.row(vertex_id).transpose());
    query["vertex"]["pose_t0"]["position"] =
        to_json(poses_t0[bodyA_id].position);
    query["vertex"]["pose_t0"]["rotation"] =
        to_json(poses_t0[bodyA_id].rotation);
    query["vertex"]["pose_t1"] = nlohmann::json();
    query["vertex"]["pose_t1"]["position"] =
        to_json(poses_t1[bodyA_id].position);
    query["vertex"]["pose_t1"]["rotation"] =
        to_json(poses_t1[bodyA_id].rotation);

    combined_query_writer.append_query(query);
}

void save_ccd_candidate(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const FaceVertexCandidate& fv_candidate)
{
    nlohmann::json query;
    long bodyA_id, vertex_id, bodyB_id, face_id;
    bodies.global_to_local_vertex(
        fv_candidate.vertex_index, bodyA_id, vertex_id);
    bodies.global_to_local_face(fv_candidate.face_index, bodyB_id, face_id);
    const auto& bodyA = bodies[bodyA_id];
    const auto& bodyB = bodies[bodyB_id];
    query["type"] = "fv";

    query["face"] = nlohmann::json();
    query["face"]["vertex0"] =
        to_json(bodyB.vertices.row(bodyB.faces(face_id, 0)).transpose());
    query["face"]["vertex1"] =
        to_json(bodyB.vertices.row(bodyB.faces(face_id, 1)).transpose());
    query["face"]["vertex2"] =
        to_json(bodyB.vertices.row(bodyB.faces(face_id, 2)).transpose());
    query["face"]["pose_t0"] = nlohmann::json();
    query["face"]["pose_t0"]["position"] =
        to_json(poses_t0[bodyB_id].position);
    query["face"]["pose_t0"]["rotation"] =
        to_json(poses_t0[bodyB_id].rotation);
    query["face"]["pose_t1"] = nlohmann::json();
    query["face"]["pose_t1"]["position"] =
        to_json(poses_t1[bodyB_id].position);
    query["face"]["pose_t1"]["rotation"] =
        to_json(poses_t1[bodyB_id].rotation);

    query["vertex"] = nlohmann::json();
    query["vertex"]["vertex"] =
        to_json(bodyA.vertices.row(vertex_id).transpose());
    query["vertex"]["pose_t0"]["position"] =
        to_json(poses_t0[bodyA_id].position);
    query["vertex"]["pose_t0"]["rotation"] =
        to_json(poses_t0[bodyA_id].rotation);
    query["vertex"]["pose_t1"] = nlohmann::json();
    query["vertex"]["pose_t1"]["position"] =
        to_json(poses_t1[bodyA_id].position);
    query["vertex"]["pose_t1"]["rotation"] =
        to_json(poses_t1[bodyA_id].rotation);

    combined_query_writer.append_query(query);
}

void save_ccd_candidate(
    const RigidBodyAssembler& bodies,
    const PosesD& poses_t0,
    const PosesD& poses_t1,
    const EdgeEdgeCandidate& ee_candidate)
{
    nlohmann::json query;
    long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
    bodies.global_to_local_edge(ee_candidate.edge0_index, bodyA_id, edgeA_id);
    bodies.global_to_local_edge(ee_candidate.edge1_index, bodyB_id, edgeB_id);
    const auto& bodyA = bodies[bodyA_id];
    const auto& bodyB = bodies[bodyB_id];
    query["type"] = "ee";

    query["edge0"] = nlohmann::json();
    query["edge0"]["vertex0"] =
        to_json(bodyA.vertices.row(bodyA.edges(edgeA_id, 0)).transpose());
    query["edge0"]["vertex1"] =
        to_json(bodyA.vertices.row(bodyA.edges(edgeA_id, 1)).transpose());
    query["edge0"]["pose_t0"] = nlohmann::json();
    query["edge0"]["pose_t0"]["position"] =
        to_json(poses_t0[bodyA_id].position);
    query["edge0"]["pose_t0"]["rotation"] =
        to_json(poses_t0[bodyA_id].rotation);
    query["edge0"]["pose_t1"] = nlohmann::json();
    query["edge0"]["pose_t1"]["position"] =
        to_json(poses_t1[bodyA_id].position);
    query["edge0"]["pose_t1"]["rotation"] =
        to_json(poses_t1[bodyA_id].rotation);

    query["edge1"] = nlohmann::json();
    query["edge1"]["vertex0"] =
        to_json(bodyB.vertices.row(bodyB.edges(edgeB_id, 0)).transpose());
    query["edge1"]["vertex1"] =
        to_json(bodyB.vertices.row(bodyB.edges(edgeB_id, 1)).transpose());
    query["edge1"]["pose_t0"] = nlohmann::json();
    query["edge1"]["pose_t0"]["position"] =
        to_json(poses_t0[bodyB_id].position);
    query["edge1"]["pose_t0"]["rotation"] =
        to_json(poses_t0[bodyB_id].rotation);
    query["edge1"]["pose_t1"] = nlohmann::json();
    query["edge1"]["pose_t1"]["position"] =
        to_json(poses_t1[bodyB_id].position);
    query["edge1"]["pose_t1"]["rotation"] =
        to_json(poses_t1[bodyB_id].rotation);

    combined_query_writer.append_query(query);
}

} // namespace ipc::rigid
