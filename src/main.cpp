// Launch simulation GUI

#include <CLI/CLI.hpp>

#include <tbb/global_control.h>
#include <tbb/task_scheduler_init.h>
#include <thread>

#include <logger.hpp>
#include <viewer/UISimState.hpp>

int main(int argc, char* argv[])
{
    using namespace ipc::rigid;
    set_logger_level(spdlog::level::info);

    CLI::App app("run simulation with viewer");

    std::string scene_path = "";
    app.add_option(
        "scene_path,-i,-s,--scene-path", scene_path,
        "JSON file with input scene");

    spdlog::level::level_enum loglevel = spdlog::level::info;
    app.add_option("--log,--loglevel", loglevel, "log level")
        ->default_val(loglevel)
        ->transform(CLI::CheckedTransformer(
            SPDLOG_LEVEL_NAMES_TO_LEVELS, CLI::ignore_case));

    int nthreads = tbb::task_scheduler_init::default_num_threads();
    app.add_option("--nthreads", nthreads, "maximum number of threads to use")
        ->default_val(nthreads);

    CLI11_PARSE(app, argc, argv);

    set_logger_level(loglevel);

    if (nthreads <= 0) {
        nthreads = tbb::task_scheduler_init::default_num_threads();
    }

    if (nthreads > tbb::task_scheduler_init::default_num_threads()) {
        spdlog::warn(
            "Attempting to use more threads than available ({:d} > {:d})!",
            nthreads, tbb::task_scheduler_init::default_num_threads());
    }

    tbb::global_control thread_limiter(
        tbb::global_control::max_allowed_parallelism, nthreads);

    UISimState ui;
    ui.launch(scene_path);
}
