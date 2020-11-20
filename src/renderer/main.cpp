#include <ctime>

#include <CLI/CLI.hpp>
#include <boost/filesystem.hpp>
#include <igl/Timer.h>
#include <nlohmann/json.hpp>
#include <tbb/parallel_for.h>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>
#include <physics/pose.hpp>
#include <physics/rigid_body_assembler.hpp>

#include "render_mesh.hpp"

using namespace swr;

struct SimRenderArgs {
    boost::filesystem::path sim_path;
    boost::filesystem::path output_path = boost::filesystem::path("sim.mp4");
    spdlog::level::level_enum loglevel = spdlog::level::level_enum::info;
};

SimRenderArgs parse_args(int argc, char* argv[])
{
    CLI::App app { "render simulation headless" };

    SimRenderArgs args;

    app.add_option(
           "sim_path,-i,-s,--sim-path", args.sim_path,
           "path to simulation JSON")
        ->required();
    app.add_option("-o,--output", args.output_path, "path to output render");
    app.add_option(
        "-l,--log,--loglevel", args.loglevel,
        "set log level 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, "
        "6=off");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        exit(app.exit(e));
    }

    return args;
}

bool read_json(const std::string& filename, nlohmann::json& json)
{
    std::ifstream input(filename);
    json = nlohmann::json::parse(input, nullptr, false);
    return !json.is_discarded();
}

int main(int argc, char* argv[])
{
    SimRenderArgs args = parse_args(argc, argv);

    ccd::logger::set_level(args.loglevel);

    ///////////////////////////////////////////////////////////////////////////
    // Create folder for PNG frames
    // Create the output directory if it does not exist
    if (args.output_path.has_parent_path()) {
        boost::filesystem::create_directories(args.output_path.parent_path());
    }
    boost::filesystem::path frames_dir = args.output_path.parent_path()
        / fmt::format("frames-{}", ccd::logger::now());
    boost::filesystem::create_directories(frames_dir);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Read the simulation json file
    nlohmann::json sim;
    if (!read_json(args.sim_path.string(), sim)) {
        spdlog::error("Invalid simulation JSON file");
        return 1;
    }

    auto state_sequence =
        sim["animation"]["state_sequence"].get<std::vector<nlohmann::json>>();
    if (state_sequence.size() == 0) {
        return 0;
    }

    std::vector<ccd::physics::RigidBody> rbs;
    ccd::io::read_rb_scene(sim["args"]["rigid_body_problem"], rbs);
    ccd::physics::RigidBodyAssembler bodies;
    bodies.init(rbs);

    const Eigen::MatrixXi &E = bodies.m_edges, &F = bodies.m_faces;
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Create a scene
    auto render_args_path = boost::filesystem::path(__FILE__).parent_path()
        / "render_settings.json";
    nlohmann::json render_args;
    if (!read_json(render_args_path.string(), render_args)) {
        spdlog::error(
            "Invalid render settings JSON file ({})",
            render_args_path.string());
        return 1;
    }
    Scene scene(render_args);
    scene.camera.align_camera_center(bodies.world_vertices(), F);
    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////
    // Per-vertex colors
    Eigen::VectorXi C(bodies.num_vertices());
    int start_i = 0;
    for (const auto& body : bodies.m_rbs) {
        C.segment(start_i, body.vertices.rows()).setConstant(int(body.type));
        start_i += body.vertices.rows();
    }
    ///////////////////////////////////////////////////////////////////////////

    tbb::parallel_for(size_t(0), state_sequence.size(), [&](size_t i) {
        ccd::physics::Poses<double> poses(bodies.num_bodies());
        assert(state_sequence[i]["rigid_bodies"].size() == bodies.num_bodies());
        for (int j = 0; j < bodies.num_bodies(); j++) {
            const auto& jrb = state_sequence[i]["rigid_bodies"][j];
            ccd::io::from_json(jrb["position"], poses[j].position);
            ccd::io::from_json(jrb["rotation"], poses[j].rotation);
        }

        std::string frame_name =
            (frames_dir / fmt::format("frame{:06d}.png", i)).string();

        igl::Timer render_timer;
        render_timer.start();
        bool wrote_frame = render_mesh(
            scene, bodies.world_vertices(poses), E, F, C, frame_name);
        render_timer.stop();

        if (wrote_frame) {
            spdlog::info(
                "Rendered frame {:d} to '{}' in {:g} seconds", //
                i, frame_name, render_timer.getElapsedTime());
        } else {
            spdlog::error("Unable to render frame {:d} to '{}'", i, frame_name);
        }
    });

    ///////////////////////////////////////////////////////////////////////////
    // Combine PNG frames to video
    int fps = render_args["fps"];
    if (fps <= 0) {
        fps = int(1 / sim["args"]["timestep"].get<double>());
    }
    std::string ffmpeg_cmd = fmt::format(
        "ffmpeg -hide_banner -loglevel warning -y -r {:d} -i {}/frame%06d.png "
        "-vcodec libx264 -crf 0 {}",
        fps, frames_dir.string(), args.output_path.string());
    spdlog::info("Combining frames using '{}'", ffmpeg_cmd);
    std::system(ffmpeg_cmd.c_str());
    ///////////////////////////////////////////////////////////////////////////
}
