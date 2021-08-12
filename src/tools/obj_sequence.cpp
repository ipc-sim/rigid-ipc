#include <CLI/CLI.hpp>
#include <ghc/fs_std.hpp> // filesystem

#include <SimState.hpp>
#include <logger.hpp>

int main(int argc, char* argv[])
{
    using namespace ipc::rigid;

    CLI::App app("Generate OBJ sequence from simulation results.");

    std::string sim_path = "";
    app.add_option(
           "sim_path,-i,-s,--sim-path", sim_path,
           "JSON file with simulation results")
        ->required();

    std::string output_dir = "";
    app.add_option(
           "output_dir,-o,--output-dir", output_dir,
           "directory for OBJ sequence")
        ->required();

    spdlog::level::level_enum loglevel = spdlog::level::warn;
    app.add_option("--log,--loglevel", loglevel, "log level")
        ->default_val(loglevel)
        ->transform(CLI::CheckedTransformer(
            SPDLOG_LEVEL_NAMES_TO_LEVELS, CLI::ignore_case));

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    set_logger_level(loglevel);

    // Create the output directory if it does not exist
    fs::create_directories(fs::path(output_dir));

    SimState sim;
    bool success = sim.load_scene(sim_path);
    if (!success) {
        return app.exit(
            CLI::Error("load_sim_failed", "Unable to load simulation result!"));
    }
    sim.save_obj_sequence(output_dir);
}
