#pragma once

#ifdef RIGID_IPC_PROFILE_FUNCTIONS

#include <igl/Timer.h>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace ipc::rigid {
namespace profiler {

    class ProfilerPoint {
    public:
        ProfilerPoint(const std::string name);

        void clear();

        void begin();
        void end();

        const std::string& name() const { return m_name; }

        size_t num_evaluations() const { return m_num_evaluations; }
        const double& total_time() const { return m_total_time; }
        const size_t& max_peak_rss_change() const
        {
            return m_max_peak_rss_change;
        }

        void message_header(const std::string& header);
        const std::string& message_header() const { return m_message_header; }
        void message(const std::string& m);
        const std::vector<std::string>& messages() const { return m_messages; }

    protected:
        std::string m_name;
        size_t m_num_evaluations;
        double m_total_time;
        size_t m_max_peak_rss_change;
        std::string m_message_header;
        std::vector<std::string> m_messages;

        igl::Timer timer;
        bool is_active;
        size_t beginning_peak_rss;
    };

    class Profiler {
    public:
        static Profiler& instance();
        void clear();

        void output_dir(const std::string& d) { dout = d; }
        std::shared_ptr<ProfilerPoint> create_point(std::string name);
        std::shared_ptr<ProfilerPoint> create_main_point(std::string name);

        void log(const std::string& fin = "");

    protected:
        std::string dout = "logs";
        Profiler() {}
        void write_summary(const std::string& dout, const std::string& fin);
        void write_point_details(
            const std::string& dout, const ProfilerPoint& point);

        std::shared_ptr<ProfilerPoint> main;
        std::vector<std::shared_ptr<ProfilerPoint>> points;
    };

    class ProfilerLog {
    public:
        ProfilerLog();
        std::shared_ptr<spdlog::logger> logger;
    };
    spdlog::logger& log();

} // namespace profiler
} // namespace ipc::rigid

#define PROFILE_MAIN_POINT(Description)                                        \
    static std::shared_ptr<profiler::ProfilerPoint> _PROFILER_POINT_ =         \
        profiler::Profiler::instance().create_main_point(Description);

#define NAMED_PROFILE_POINT(Description, Name)                                 \
    static std::shared_ptr<profiler::ProfilerPoint> _PROFILER_POINT_##Name =   \
        profiler::Profiler::instance().create_point(Description);

#define PROFILE_POINT(Description)                                             \
    static std::shared_ptr<profiler::ProfilerPoint> _PROFILER_POINT_ =         \
        profiler::Profiler::instance().create_point(Description);

#define PROFILE_START(Name) _PROFILER_POINT_##Name->begin();

#define PROFILE_END(Name) _PROFILER_POINT_##Name->end();

#define PROFILE_MESSAGE(Name, Header, Val)                                     \
    _PROFILER_POINT_##Name->message_header(Header);                            \
    _PROFILER_POINT_##Name->message(Val);
#define PROFILER_CLEAR() profiler::Profiler::instance().clear();
#define LOG_PROFILER(SceneFile) profiler::Profiler::instance().log(SceneFile);
#define PROFILER_OUTDIR(OutputDir)                                             \
    profiler::Profiler::instance().output_dir(OutputDir);
#else

#define PROFILE_MAIN_POINT(Description) ;
#define NAMED_PROFILE_POINT(Description, Name) ;
#define PROFILE_POINT(Description) ;
#define PROFILE_START(Name) ;
#define PROFILE_END(Name) ;
#define PROFILE_MESSAGE(Name, Header, Val) ;
#define PROFILER_CLEAR() ;
#define LOG_PROFILER(SceneFile) ;
#define PROFILER_OUTDIR(OutputDir) ;

#endif
