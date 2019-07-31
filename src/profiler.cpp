#include <logger.hpp>
#include <profiler.hpp>

#ifdef PROFILE_FUNCTIONS

namespace ccd {
namespace profiler {

    Profiler& Profiler::instance()
    {
        static Profiler profiler;
        return profiler;
    }
    std::shared_ptr<ProfilerPoint> Profiler::create_point(std::string name)
    {
        auto point = std::make_shared<ProfilerPoint>(name);
        points.push_back(point);
        return point;
    }
    std::shared_ptr<ProfilerPoint> Profiler::create_main_point(std::string name)
    {
        main = std::make_shared<ProfilerPoint>(name);
        return main;
    }
    void Profiler::log(const std::string& fin)
    {

        std::string parent_name
            = fmt::format("{}/log_{}", LOGS_OUTPUT_DIR, log::now());
        if (mkdir(parent_name.c_str(), ACCESSPERMS) == 0) {
            write_summary(parent_name, fin);
            for (auto& p : points) {
                write_point(parent_name, *p);
            }
        }
    }

    void Profiler::write_summary(
        const std::string& dout, const std::string& fin)
    {
        std::string filename = fmt::format("{}/summary.csv", dout);

        std::ofstream myfile;
        myfile.open(filename);
        myfile << fin << "\n";
        myfile
            << "section, total_time (sec), percentage_time, num_calls, avg_time (sec), num_success\n";

        double total_time = main->total_time();
        size_t num_calls = main->num_evaluations();
        myfile << fmt::format("{}, {:10e}, 100.00%, {}, {}, -1\n", main->name(),
            total_time, num_calls, total_time / num_calls);

        for (auto& p : points) {
            double p_time = p->total_time();
            size_t p_num_calls = p->num_evaluations();
            myfile << fmt::format("{}, {:10e}, {:2f}%, {}, {}, {} \n",
                p->name(), p_time, p_time / total_time * 100, p_num_calls,
                p_time / p_num_calls, p->num_success());
        }
        myfile.close();
    } // namespace profiler

    void Profiler::write_point(
        const std::string& dout, const ProfilerPoint& point)
    {
        std::string filename = fmt::format("{}/{}.csv", dout, point.name());

        std::ofstream myfile;

        auto& message_history = point.messages();
        auto& time_history = point.time();
        if (message_history.size() != time_history.size()) {
            return;
        }

        myfile.open(filename);
        myfile << "time (sec), message\n";
        if (message_history.size()) {
            for (size_t i = 0; i < time_history.size(); ++i) {
                myfile << fmt::format(
                    "{:10e},{}\n", time_history[i], message_history[i]);
            }
        }
        myfile.close();
    }
    ProfilerLog::ProfilerLog()
    {

        auto console_sink
            = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::debug);

        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
            fmt::format("{}/profiler.log", LOGS_OUTPUT_DIR), false);
        file_sink->set_level(spdlog::level::debug);

        spdlog::sinks_init_list sinks = { console_sink, file_sink };
        spdlog::logger logger_("multi_sink", sinks);
        logger_.set_level(spdlog::level::debug);

        this->logger = std::make_shared<spdlog::logger>("profile", sinks);
        this->logger->set_level(spdlog::level::debug);
        this->logger->debug("BEGIN SESSION");
    }

    spdlog::logger& log()
    {
        static ProfilerLog log;
        return *log.logger;
    }

} // namespace profiler
} // namespace ccd
#endif
