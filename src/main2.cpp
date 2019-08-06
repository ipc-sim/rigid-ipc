#include <SimState.hpp>
#include <logger.hpp>

int main(int argc, char* argv[])
{
    spdlog::set_level(spdlog::level::info);
    ccd::SimState sim;

    // TODO: add proper CLI
    if (argc >= 3) {
        sim.load_scene(argv[1]);
        if (argc >= 4) {
            int max_iter = atoi(argv[3]);
            sim.m_max_simulation_steps = max_iter;
            spdlog::info("running {} iterations", max_iter);
        }
        sim.run_simulation(argv[2]);
    } else {
        spdlog::error("provided only {} arguments ", argc);
    }
}
