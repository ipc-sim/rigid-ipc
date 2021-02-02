#include "write_gltf.hpp"

#define TINYGLTF_IMPLEMENTATION
#define TINYGLTF_NO_INCLUDE_JSON
#define TINYGLTF_NO_STB_IMAGE
#define TINYGLTF_NO_STB_IMAGE_WRITE
#include <nlohmann/json.hpp>
#include <tiny_gltf.h>

namespace ipc::rigid {

bool write_gltf(
    const std::string& filename,
    const RigidBodyAssembler& bodies,
    const std::vector<PosesD>& poses,
    double timestep,
    bool embed_buffers,
    bool write_binary,
    bool prettyPrint)
{
    typedef float Float;
    int float_component_type;
    if (std::is_same<Float, float>::value) {
        float_component_type = TINYGLTF_COMPONENT_TYPE_FLOAT;
    } else {
        // assert(std::is_same<Float, double>::value);
        float_component_type = TINYGLTF_COMPONENT_TYPE_DOUBLE;
    }

    using namespace tinygltf;

    assert(bodies.dim() == 3);
    size_t num_bodies = bodies.num_bodies();
    assert(poses.size() > 0);
    size_t num_steps = poses.size();

    Model model;
    model.defaultScene = 0;

    model.asset.version = "2.0";
    model.asset.generator = "RigidIPC";

    Scene& scene = model.scenes.emplace_back();
    scene.name = "RigidIPCSimulation";
    scene.nodes.resize(num_bodies);
    std::iota(scene.nodes.begin(), scene.nodes.end(), 0);

    model.nodes.resize(num_bodies);
    model.animations.resize(1);
    Animation& animation = model.animations[0];
    animation.name = "Simulation";
    animation.channels.resize(2 * num_bodies);
    animation.samplers.resize(2 * num_bodies);
    model.meshes.resize(num_bodies);
    // Four accessors per bidy: vertices, faces, translations, and rotations
    model.accessors.resize(4 * num_bodies + 1);
    {
        Accessor& accessor = model.accessors[2 * num_bodies];
        accessor.name = "Times";
        accessor.bufferView = 2 * num_bodies;
        accessor.componentType = float_component_type;
        accessor.count = num_steps;
        accessor.minValues.push_back(timestep);
        accessor.maxValues.push_back(num_steps * timestep);
        accessor.type = TINYGLTF_TYPE_SCALAR;
    }

    for (int i = 0; i < num_bodies; i++) {
        Node& node = model.nodes[i];
        node.mesh = i;
        node.name = bodies[i].name;
        node.translation = { { poses[0][i].position.x(),
                               poses[0][i].position.y(),
                               poses[0][i].position.z() } };
        Eigen::Quaternion<double> q = poses[0][i].construct_quaternion();
        node.rotation = { { q.x(), q.y(), q.z(), q.w() } };

        animation.channels[2 * i + 0].sampler = 2 * i;
        animation.channels[2 * i + 0].target_node = i;
        animation.channels[2 * i + 0].target_path = "translation";
        animation.channels[2 * i + 1].sampler = 2 * i + 1;
        animation.channels[2 * i + 1].target_node = i;
        animation.channels[2 * i + 1].target_path = "rotation";

        animation.samplers[2 * i + 0].input = 2 * num_bodies;
        animation.samplers[2 * i + 0].output = 2 * num_bodies + 2 * i + 1;
        animation.samplers[2 * i + 0].interpolation = "LINEAR";
        animation.samplers[2 * i + 1].input = 2 * num_bodies;
        animation.samplers[2 * i + 1].output = 2 * num_bodies + 2 * i + 2;
        animation.samplers[2 * i + 1].interpolation = "LINEAR";

        Mesh& mesh = model.meshes[i];
        mesh.name = bodies[i].name;
        Primitive& primitive = mesh.primitives.emplace_back();
        primitive.attributes["POSITION"] = 2 * i;
        primitive.indices = 2 * i + 1;
        primitive.mode = TINYGLTF_MODE_TRIANGLES;

        Accessor* accessor = &model.accessors[2 * i];
        accessor->name = bodies[i].name + "Vertices";
        accessor->bufferView = 2 * i;
        accessor->componentType = float_component_type;
        accessor->count = bodies[i].num_vertices();
        // accessor->max = ...;
        // accessor->min = ...;
        accessor->type = TINYGLTF_TYPE_VEC3;

        accessor = &model.accessors[2 * i + 1];
        accessor->name = bodies[i].name + "Faces";
        accessor->bufferView = 2 * i + 1;
        accessor->componentType = TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT;
        accessor->count = 3 * bodies[i].num_faces();
        accessor->type = TINYGLTF_TYPE_SCALAR;

        accessor = &model.accessors[2 * num_bodies + 2 * i + 1];
        accessor->name = bodies[i].name + "Translations";
        accessor->bufferView = 2 * num_bodies + 2 * i + 1;
        accessor->componentType = float_component_type;
        accessor->count = num_steps;
        accessor->type = TINYGLTF_TYPE_VEC3;

        accessor = &model.accessors[2 * num_bodies + 2 * i + 2];
        accessor->name = bodies[i].name + "Rotations";
        accessor->bufferView = 2 * num_bodies + 2 * i + 2;
        accessor->componentType = float_component_type;
        accessor->count = num_steps;
        accessor->type = TINYGLTF_TYPE_VEC4;

        BufferView& buffer_view = model.bufferViews[2 * i];
    }

    ///////////////////////////////////////////////////////////////////////////

    model.bufferViews.resize(4 * num_bodies + 1);
    size_t byte_offset = 0;
    for (int i = 0; i < num_bodies; i++) {
        BufferView* buffer_view = &model.bufferViews[2 * i];
        buffer_view->name = bodies[i].name + "Vertices";
        buffer_view->buffer = 0;
        buffer_view->byteLength = sizeof(Float) * bodies[i].vertices.size();
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;

        buffer_view = &model.bufferViews[2 * i + 1];
        buffer_view->name = bodies[i].name + "Faces";
        buffer_view->buffer = 0;
        buffer_view->byteLength = sizeof(unsigned int) * bodies[i].faces.size();
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;
    }
    {
        BufferView& buffer_view = model.bufferViews[2 * num_bodies];
        buffer_view.name = "Times";
        buffer_view.buffer = 0;
        buffer_view.byteLength = num_steps * sizeof(Float);
        buffer_view.byteOffset = byte_offset;
        byte_offset += buffer_view.byteLength;
    }
    for (int i = 0; i < num_bodies; i++) {
        BufferView* buffer_view =
            &model.bufferViews[2 * num_bodies + 2 * i + 1];
        buffer_view->name = bodies[i].name + "Translations";
        buffer_view->buffer = 0;
        buffer_view->byteLength = 3 * sizeof(Float) * num_steps;
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;

        buffer_view = &model.bufferViews[2 * num_bodies + 2 * i + 2];
        buffer_view->name = bodies[i].name + "Rotations";
        buffer_view->buffer = 0;
        buffer_view->byteLength = 4 * sizeof(Float) * num_steps;
        buffer_view->byteOffset = byte_offset;
        byte_offset += buffer_view->byteLength;
    }

    ///////////////////////////////////////////////////////////////////////////

    std::vector<unsigned char> byte_data(byte_offset);
    size_t byte_i = 0;
    for (int i = 0; i < num_bodies; i++) {
        Eigen::MatrixXd V = bodies[i].vertices;
        for (int r = 0; r < V.rows(); r++) {
            for (int c = 0; c < V.cols(); c++) {
                Float v = V(r, c);
                std::memcpy(&byte_data[byte_i], &v, sizeof(Float));
                byte_i += sizeof(Float);
            }
        }

        Eigen::MatrixXi F = bodies[i].faces;
        for (int r = 0; r < F.rows(); r++) {
            for (int c = 0; c < F.cols(); c++) {
                unsigned int fij = F(r, c);
                std::memcpy(&byte_data[byte_i], &fij, sizeof(unsigned int));
                byte_i += sizeof(unsigned int);
            }
        }
    }
    for (int i = 0; i < num_steps; i++) {
        Float t = i * timestep;
        std::memcpy(&byte_data[byte_i], &t, sizeof(Float));
        byte_i += sizeof(Float);
    }
    for (int i = 0; i < num_bodies; i++) {
        for (int j = 0; j < num_steps; j++) {
            Eigen::Vector3d p = poses[j][i].position;
            for (int d = 0; d < p.size(); d++) {
                Float pd = p[d];
                std::memcpy(&byte_data[byte_i], &pd, sizeof(Float));
                byte_i += sizeof(Float);
            }
        }

        for (int j = 0; j < num_steps; j++) {
            Eigen::Quaternion<double> quat = poses[j][i].construct_quaternion();
            Eigen::Vector4d q(quat.x(), quat.y(), quat.z(), quat.w());
            for (int d = 0; d < q.size(); d++) {
                Float qd = q[d];
                std::memcpy(&byte_data[byte_i], &qd, sizeof(Float));
                byte_i += sizeof(Float);
            }
        }
    }
    assert(byte_i == byte_data.size());

    model.buffers.emplace_back().data = byte_data;
    return TinyGLTF().WriteGltfSceneToFile(
        &model, filename,
        /*embedImages=*/true, embed_buffers, prettyPrint, write_binary);
}

} // namespace ipc::rigid
