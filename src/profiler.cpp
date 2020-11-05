#include <logger.hpp>
#include <map>
#include <profiler.hpp>

#ifdef PROFILE_FUNCTIONS

namespace ccd {
namespace profiler {

    // -----------------------------------------------------------------
    // PROFILER POINT
    // -----------------------------------------------------------------

    ProfilerPoint::ProfilerPoint(const std::string name)
        : name_(name)
    {
    }
    void ProfilerPoint::begin(const std::vector<std::string>& stack)
    {
        stacks_.push_back(stack);
        timer.start();
    }
    void ProfilerPoint::end() { times_.push_back(timer.getElapsedTime()); }
    void ProfilerPoint::message(const std::string& m)
    {
        messages_.push_back(m);
    }
    void ProfilerPoint::clear()
    {
        times_.clear();
        stacks_.clear();
        messages_.clear();
    }
    double ProfilerPoint::total_time() const
    {
        return std::accumulate(times_.begin(), times_.end(), 0.0);
    }

    // -----------------------------------------------------------------
    // PROFILER MANAGER
    // -----------------------------------------------------------------

    Profiler& Profiler::instance()
    {
        static Profiler profiler;
        return profiler;
    }

    void Profiler::clear()
    {
        if (main != nullptr) {
            main->clear();
        }
        for (auto& p : points) {
            p->clear();
        }
        stack_.clear();
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

    void Profiler::push(const std::shared_ptr<ProfilerPoint>& point)
    {
        stack_.push_back(point->name());
    }

    void Profiler::pop() { stack_.pop_back(); }

    void Profiler::log(const std::string& fin)
    {

        std::string parent_name =
            fmt::format("{}/log_{}", dout, ccd::logger::now());
        if (mkdir(parent_name.c_str(), ACCESSPERMS) == 0) {
            write_summary(parent_name, fin);
            for (auto& p : points) {
                //                write_point_summary(parent_name, *p);
                write_point_details(parent_name, *p);
            }
        }
    }

    void
    Profiler::write_summary(const std::string& dout, const std::string& fin)
    {
        // create more detailed map of function calls:
        typedef std::tuple<int, double> T;
        typedef std::pair<std::string, T> M;
        std::map<std::string, T> num_events_per_path;

        for (auto& point : points) {
            /// add up the events of each stack
            auto& time_history = point->time();
            auto& stack_history = point->stacks();
            for (size_t i = 0; i < time_history.size(); ++i) {
                std::string stack = "";
                std::for_each(
                    stack_history[i].begin(), stack_history[i].end(),
                    [&](const std::string& piece) { stack += "->" + piece; });

                auto time = time_history[i];
                auto el = num_events_per_path.find(stack);
                if (el != num_events_per_path.end()) { // found
                    std::get<0>(el->second) += 1;
                    std::get<1>(el->second) += time;
                } else {
                    num_events_per_path.insert(M(stack, T(1, time)));
                }
            }
        }

        std::string filename = fmt::format("{}/summary.csv", dout);

        std::ofstream myfile;
        myfile.open(filename);
        myfile << fin << "\n";
        myfile << "section,total_time (sec),percentage_time,num_calls,avg_time "
                  "(sec)\n";

        double total_time = main->total_time();
        size_t num_calls = main->num_evaluations();
        myfile << fmt::format(
            "{},{:10e},100.00%,{},{}\n", main->name(), total_time, num_calls,
            total_time / num_calls);

        for (auto& p : points) {
            double p_time = p->total_time();
            size_t p_num_calls = p->num_evaluations();
            myfile << fmt::format(
                "{},{:10e},{:2f}%,{},{}\n", p->name(), p_time,
                p_time / total_time * 100, p_num_calls, p_time / p_num_calls);
        }
        myfile << "\n\n";

        for (auto& it : num_events_per_path) {
            int p_num_calls = std::get<0>(it.second);
            double p_time = std::get<1>(it.second);

            myfile << fmt::format(
                "{},{:10e},{:2f}%,{},{}\n", it.first, p_time,
                p_time / total_time * 100, p_num_calls, p_time / p_num_calls);
        }

        myfile.close();
    } // namespace profiler

    void Profiler::write_point_details(
        const std::string& dout, const ProfilerPoint& point)
    {
        std::string filename =
            fmt::format("{}/{}_details.csv", dout, point.name());

        std::ofstream myfile;

        auto& message_history = point.messages();
        auto& time_history = point.time();
        auto& stack_history = point.stacks();

        myfile.open(filename);
        myfile << "stack,time (sec),message\n";

        for (size_t i = 0; i < time_history.size(); ++i) {
            std::string stack = "";
            std::for_each(
                stack_history[i].begin(), stack_history[i].end(),
                [&](const std::string& piece) { stack += "->" + piece; });

            auto time = time_history[i];
            std::string message =
                message_history.size() > i ? message_history[i] : "";

            myfile << fmt::format("{},{:10e},{}\n", stack, time, message);
        }

        myfile.close();
    }
    void Profiler::write_point_summary(
        const std::string& dout, const ProfilerPoint& point)
    {

        typedef std::tuple<int, double> T;
        typedef std::pair<std::string, T> M;
        std::map<std::string, T> num_events_per_path;

        /// add up the events of each stack
        auto& time_history = point.time();
        auto& stack_history = point.stacks();
        for (size_t i = 0; i < time_history.size(); ++i) {
            std::string stack = "";
            std::for_each(
                stack_history[i].begin(), stack_history[i].end(),
                [&](const std::string& piece) { stack += "->" + piece; });

            auto time = time_history[i];
            auto el = num_events_per_path.find(stack);
            if (el != num_events_per_path.end()) { // found
                std::get<0>(el->second) += 1;
                std::get<1>(el->second) += time;
            } else {
                num_events_per_path.insert(M(stack, T(1, time)));
            }
        }
        std::string filename =
            fmt::format("{}/{}_summary.csv", dout, point.name());
        std::ofstream myfile;
        myfile.open(filename);
        myfile << "stack,total_time (sec),num_calls\n";
        for (auto& it : num_events_per_path) {
            myfile << fmt::format(
                "{},{:10e},{}\n", it.first, std::get<1>(it.second),
                std::get<0>(it.second));
        }
        myfile.close();
    }

    // -----------------------------------------------------------------
    // PROFILER LOG
    // -----------------------------------------------------------------
    ProfilerLog::ProfilerLog()
    {
        auto console_sink =
            std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
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
