#include <ghc/fs_std.hpp> // filesystem
#include <logger.hpp>
#include <map>
#include <profiler.hpp>
#include <utils/get_rss.hpp>

#ifdef RIGID_IPC_PROFILE_FUNCTIONS

namespace ipc::rigid {
namespace profiler {

    // -----------------------------------------------------------------
    // PROFILER POINT
    // -----------------------------------------------------------------

    ProfilerPoint::ProfilerPoint(const std::string name)
        : m_name(name)
    {
        clear();
    }

    void ProfilerPoint::begin()
    {
        if (is_active) {
            throw fmt::format("ProfilePoint {} is already active!", name());
        }
        timer.start();
        beginning_peak_rss = getPeakRSS();
        is_active = true;
    }

    void ProfilerPoint::end()
    {
        if (is_active) {
            timer.stop();
            m_total_time += timer.getElapsedTime();
            m_num_evaluations++;
            m_max_peak_rss_change = std::max(
                m_max_peak_rss_change, getPeakRSS() - beginning_peak_rss);
            is_active = false;
        }
    }

    void ProfilerPoint::message_header(const std::string& header)
    {
        m_message_header = header;
    }

    void ProfilerPoint::message(const std::string& m)
    {
        m_messages.push_back(m);
    }

    void ProfilerPoint::clear()
    {
        m_num_evaluations = 0;
        m_total_time = 0;
        m_max_peak_rss_change = 0;
        m_message_header = "";
        m_messages.clear();
        is_active = false;
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
        fs::path outpath(dout);
        outpath /= fmt::format("log-{}", current_time_string());
        fs::create_directories(outpath);

        write_summary(outpath.string(), fin);
        for (auto& p : points) {
            write_point_details(outpath.string(), *p);
        }
    }

    void
    Profiler::write_summary(const std::string& dout, const std::string& fin)
    {
        std::string filename = fmt::format("{}/summary.csv", dout);

        std::ofstream myfile;
        myfile.open(filename);
        myfile << fin << "\n";
        myfile << "section,total_time (sec),percentage_time,num_calls,"
                  "avg_time (sec),max peak RSS change (KB)\n";

        double total_time = main->total_time();
        size_t num_calls = main->num_evaluations();
        double max_peak_rss_change = main->max_peak_rss_change() / 1024.0;
        myfile << fmt::format(
            "{},{:10e},100.00%,{},{},{:g}\n", main->name(), total_time,
            num_calls, total_time / num_calls, max_peak_rss_change);

        for (auto& p : points) {
            double p_time = p->total_time();
            size_t p_num_calls = p->num_evaluations();
            double p_max_peak_rss_change = p->max_peak_rss_change() / 1024.0;
            myfile << fmt::format(
                "{},{:10e},{:2f}%,{},{},{:g}\n", p->name(), p_time,
                p_time / total_time * 100, p_num_calls, p_time / p_num_calls,
                p_max_peak_rss_change);
        }

        myfile.close();
    }

    void Profiler::write_point_details(
        const std::string& dout, const ProfilerPoint& point)
    {
        auto& messages = point.messages();
        if (messages.empty()) {
            return;
        }

        std::string point_name = point.name();
        std::replace(point_name.begin(), point_name.end(), ':', '_');
        std::string filename =
            fmt::format("{}/{}_details.csv", dout, point_name);

        std::ofstream myfile;
        myfile.open(filename);

        myfile << point.message_header() << "\n";
        for (const std::string& message : messages) {
            myfile << message << "\n";
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
            "logs/profiler.log", false);
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
} // namespace ipc::rigid
#endif
