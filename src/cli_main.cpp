#include <SimState.hpp>
#include <logger.hpp>
#include <profiler.hpp>

#include <CLI/CLI.hpp>

int main(int argc, char* argv[])
{
    spdlog::set_level(spdlog::level::info);
    ccd::SimState sim;

    struct {
        std::string scene_path = "";
        std::string output_dir = "";
        std::string output_name = "sim.json";
        int num_it = -1;
        bool is_log_trace = false;
        int loglevel = 2; // info
    } args;

    CLI::App app { "run headless simulation" };

    app.add_option("scene_path,-s,--scene-path", args.scene_path,
           "JSON file with input scene.")
        ->required();

    app.add_option("output_dir,-o,--output-path", args.output_dir,
           "directory for results.")
        ->required();
    app.add_option("-f,--output-name", args.output_name, "name for simulation file");
    app.add_option("--num-iterations", args.num_it, "number of iterations");
    app.add_flag("--trace", args.is_log_trace, "log everything information");
    app.add_option("--loglevel", args.loglevel, "set log level 0=trace, 1=debug, 2=info, 3=warn");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }
    if (args.is_log_trace) {
        spdlog::set_level(spdlog::level::trace);
    }
    else {
        spdlog::set_level(static_cast<spdlog::level::level_enum>(args.loglevel));
    }

    PROFILER_OUTDIR(args.output_dir)
    std::string fout = fmt::format("{}/{}", args.output_dir, args.output_name);

    sim.load_scene(args.scene_path);
    if (args.num_it > 0) {
        sim.m_max_simulation_steps = args.num_it;
    }

    sim.run_simulation(fout);
    spdlog::info(
        "To postprocess run:\n `python tools/results_to_vtk_files.py {} {}`",
        fout, args.output_dir);
}
