#pragma once

#ifdef PROFILE_FUNCTIONS

#include <igl/Timer.h>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

namespace ccd {
namespace profiler {

    class ProfilerPoint {
    public:
        ProfilerPoint(const std::string name);

        void clear();

        void begin() { timer.start(); }
        void begin(const std::vector<std::string>& stack);
        void end();

        void message(const std::string& m);
        size_t num_evaluations() const { return times_.size(); }
        double total_time() const;

        const std::string& name() const { return name_; }
        const std::vector<double>& time() const { return times_; }
        const std::vector<std::string>& messages() const { return messages_; }
        const std::vector<std::vector<std::string>>& stacks() const
        {
            return stacks_;
        }

    protected:
        std::vector<double> times_;
        std::vector<std::string> messages_;
        std::vector<std::vector<std::string>> stacks_;

        igl::Timer timer;
        std::string name_;
    };

    class Profiler {
    public:
        static Profiler& instance();
        void clear();

        void push(const std::shared_ptr<ProfilerPoint>& point);
        void pop();
        const std::vector<std::string>& stack() { return stack_; }
        std::shared_ptr<ProfilerPoint> create_point(std::string name);
        std::shared_ptr<ProfilerPoint> create_main_point(std::string name);

        void log(const std::string& fin = "");

    protected:
        Profiler() {}
        void write_summary(const std::string& dout, const std::string& fin);
        void write_point_summary(
            const std::string& dout, const ProfilerPoint& point);
        void write_point_details(
            const std::string& dout, const ProfilerPoint& point);

        std::vector<std::string> stack_;
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
} // namespace ccd

#define PROFILE_MAIN_POINT(Description)                                        \
    static std::shared_ptr<ccd::profiler::ProfilerPoint> _PROFILER_POINT_      \
        = ccd::profiler::Profiler::instance().create_main_point(Description);

#define NAMED_PROFILE_POINT(Description, Name)                                 \
    static std::shared_ptr<ccd::profiler::ProfilerPoint>                       \
        _PROFILER_POINT_##Name                                                 \
        = ccd::profiler::Profiler::instance().create_point(Description);

#define PROFILE_POINT(Description)                                             \
    static std::shared_ptr<ccd::profiler::ProfilerPoint> _PROFILER_POINT_      \
        = ccd::profiler::Profiler::instance().create_point(Description);

#define PROFILE_START(Name)                                                    \
    ccd::profiler::Profiler::instance().push(_PROFILER_POINT_##Name);          \
    _PROFILER_POINT_##Name->begin(ccd::profiler::Profiler::instance().stack());

#define PROFILE_END(Name)                                                      \
    _PROFILER_POINT_##Name->end();                                             \
    ccd::profiler::Profiler::instance().pop();

#define PROFILE_MESSAGE(Name, Val) _PROFILER_POINT_##Name->message(Val);
#define PROFILER_CLEAR() ccd::profiler::Profiler::instance().clear();
#define LOG_PROFILER(SceneFile)                                                \
    ccd::profiler::Profiler::instance().log(SceneFile);

#else

#define PROFILE_MAIN_POINT(Description) ;
#define NAMED_PROFILE_POINT(Description, Name) ;
#define PROFILE_POINT(Description) ;
#define PROFILE_START(Name) ;
#define PROFILE_END(Name) ;
#define PROFILE_MESSAGE(Name, Val) ;
#define PROFILER_CLEAR() ;
#define LOG_PROFILER(SceneFile) ;

#endif
