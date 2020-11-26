// Launch simulation GUI

#include <Eigen/Core>
#include <tbb/global_control.h>

#include <logger.hpp>
#include <viewer/UISimState.hpp>

#include <io/json_to_mjcf.hpp>

bool is_number(const std::string& s)
{
    return !s.empty()
        && std::find_if(
               s.begin(), s.end(),
               [](unsigned char c) { return !std::isdigit(c); })
        == s.end();
}

int main(int argc, char* argv[])
{
    spdlog::debug(
        "Using Eigen {}.{}.{}", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION,
        EIGEN_MINOR_VERSION);

    bool is_arg1_number = argc > 1 && is_number(argv[1]);

    if (is_arg1_number) {
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

    // tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
    ccd::logger::set_level(spdlog::level::info);
    ccd::UISimState ui;
    ui.launch((argc > 1 && !is_arg1_number) ? argv[1] : "");
}
