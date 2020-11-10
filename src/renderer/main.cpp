#include <ctime>

#include <CLI/CLI.hpp>
#include <boost/filesystem.hpp>
#include <fmt/chrono.h>
#include <igl/Timer.h>
#include <nlohmann/json.hpp>
#include <tbb/parallel_for.h>

#include <io/read_rb_scene.hpp>
#include <io/serialize_json.hpp>
#include <logger.hpp>
#include <physics/pose.hpp>
#include <physics/rigid_body_assembler.hpp>
#include <renderer/render_mesh.hpp>

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

int main(int argc, char* argv[])
{
    SimRenderArgs args = parse_args(argc, argv);

    ccd::logger::set_level(args.loglevel);

    // Create the output directory if it does not exist
    if (args.output_path.has_parent_path()) {
        boost::filesystem::create_directories(args.output_path.parent_path());
    }
    boost::filesystem::path frames_dir = args.output_path.parent_path()
        / fmt::format("frames-{}", ccd::logger::now());
    boost::filesystem::create_directories(frames_dir);

    std::ifstream input(args.sim_path.string());
    nlohmann::json sim = nlohmann::json::parse(input, nullptr, false);
    if (sim.is_discarded()) {
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

    std::string render_args_fname =
        (boost::filesystem::path(__FILE__).parent_path()
         / "render_settings.json")
            .string();
    std::ifstream render_args_file(render_args_fname);
    nlohmann::json render_args =
        nlohmann::json::parse(render_args_file, nullptr, false);
    if (render_args.is_discarded()) {
        spdlog::error(
            "Invalid render settings JSON file ({})", render_args_fname);
        return 1;
    }
    Scene scene(render_args);

    Eigen::MatrixXd C(bodies.num_vertices(), 3);
    Eigen::Vector3d free_color, fixed_color;
    ccd::io::from_json(render_args["free_body_color"], free_color);   // #E74C3C
    ccd::io::from_json(render_args["fixed_body_color"], fixed_color); // #B3B3B3
    Eigen::VectorXb is_vertex_fixed = bodies.is_dof_fixed.rowwise().all();
    for (int i = 0; i < bodies.num_vertices(); i++) {
        C.row(i) = is_vertex_fixed(i) ? fixed_color : free_color;
    }

    Eigen::MatrixXi E = bodies.m_edges, F = bodies.m_faces;
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

    int fps = render_args["fps"];
    if (fps <= 0) {
        fps = int(1 / sim["args"]["timestep"].get<double>());
    }
    std::string ffmpeg_cmd = fmt::format(
        "ffmpeg -y -r {:d} -i {}/frame%06d.png -vcodec libx264 -crf 0 {}", fps,
        frames_dir.string(), args.output_path.string());
    spdlog::info("Combining frames using '{}'", ffmpeg_cmd);
    std::system(ffmpeg_cmd.c_str());
}
