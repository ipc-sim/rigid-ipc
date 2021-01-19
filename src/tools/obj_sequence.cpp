#include <CLI/CLI.hpp>
#include <boost/filesystem.hpp>

#include <SimState.hpp>
#include <logger.hpp>

int main(int argc, char* argv[])
{
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

    int loglevel = 3; // info
    app.add_option(
        "--log,--loglevel", loglevel,
        "set log level 0=trace, 1=debug, 2=info, 3=warn, 4=error, 5=critical, "
        "6=off",
        true);

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return app.exit(e);
    }

    ccd::logger::set_level(static_cast<spdlog::level::level_enum>(loglevel));

    // Create the output directory if it does not exist
    boost::filesystem::create_directories(boost::filesystem::path(output_dir));

    ccd::SimState sim;
    bool success = sim.load_scene(sim_path);
    if (!success) {
        return app.exit(
            CLI::Error("load_sim_failed", "Unable to load simulation result!"));
    }
    sim.save_obj_sequence(output_dir);
}
