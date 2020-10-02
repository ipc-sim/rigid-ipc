// Launch simulation GUI

#include <Eigen/Core>

#include <logger.hpp>
#include <viewer/UISimState.hpp>

#include <io/json_to_mjcf.hpp>

int main(int argc, char* argv[])
{
    spdlog::debug(
        "Using Eigen {}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION,
        EIGEN_MINOR_VERSION);

    if (argc > 1) {
        int mode = std::stoi(argv[1]);
        switch (mode) {
        case 1: { // json to MJCF
            spdlog::info("transform json to MJCF");

            if (argc < 3) {
                spdlog::error("need input json file path!");
                return -1;
            }

            return json_to_mjcf(argv[2]);
        }

        case 2: { // process bullet output
            spdlog::info(
                "generating bullet output geometry from transformation");

            if (argc < 4) {
                spdlog::error(
                    "need input json file and bullet output file path!");
                return -1;
            }

            return generate_bullet_results(argv[2], argv[3]);
        }

        case 0: // simulation
        default:
            break;
        }
    }

    // tbb::task_scheduler_init init(1);
    ccd::logger::set_level(spdlog::level::info);
    ccd::UISimState ui;
    ui.launch();
}
