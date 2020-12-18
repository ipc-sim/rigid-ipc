#include <CLI/CLI.hpp>
#include <boost/filesystem.hpp>

#include <SimState.hpp>
#include <logger.hpp>
#include <profiler.hpp>

int main(int argc, char* argv[])
{
    ccd::logger::set_level(spdlog::level::info);
    ccd::SimState sim;

    struct {
        std::string scene_path = "";
        std::string output_dir = "";
        std::string output_name = "sim.json";
        int num_steps = -1;
        int checkpoint_freq = -1;
        int loglevel = 2; // info
    } args;

    CLI::App app { "run headless simulation" };

    app.add_option(
           "scene_path,-i,-s,--scene-path", args.scene_path,
           "JSON file with input scene")
        ->required();

    app.add_option(
           "output_dir,-o,--output-path", args.output_dir,
           "directory for results")
        ->required();
    app.add_option(
        "-f,--output-name", args.output_name, "name for simulation file", true);
    app.add_option("--num-steps", args.num_steps, "number of time-steps");
    app.add_option(
        "--chkpt,--checkpoint-frequency", args.checkpoint_freq,
        "number of time-steps between checkpoints");
    app.add_option(
        "--log,--loglevel", args.loglevel,
        "set log level 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, "
        "6=off",
        true);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    ccd::logger::set_level(
        static_cast<spdlog::level::level_enum>(args.loglevel));

    // Create the output directory if it does not exist
    boost::filesystem::create_directories(
        boost::filesystem::path(args.output_dir));
    PROFILER_OUTDIR(args.output_dir)
    std::string fout = fmt::format("{}/{}", args.output_dir, args.output_name);

    sim.load_scene(args.scene_path);
    if (args.num_steps > 0) {
        sim.m_max_simulation_steps = args.num_steps;
    }

    if (args.checkpoint_freq > 0) {
        sim.m_checkpoint_frequency = args.checkpoint_freq;
    }

    sim.run_simulation(fout);
}
