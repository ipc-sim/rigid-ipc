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
        ProfilerPoint(const std::string name)
            : name_(name)
        {
        }
        void begin() { timer.start(); }
        void end() { each_time.push_back(timer.getElapsedTime()); }
        void success(const bool val) { success_.push_back(val); }
        void message(const std::string& m) { messages_.push_back(m); }
        size_t num_evaluations() const { return each_time.size(); }
        double total_time() const
        {
            return std::accumulate(each_time.begin(), each_time.end(), 0.0);
        }
        long num_success() const
        {
            if (success_.size() > 0) {
                return std::count(success_.begin(), success_.end(), true);
            } else
                return -1;
        }
        const std::string& name() const { return name_; }
        const std::vector<double>& time() const { return each_time; }
        const std::vector<bool>& success() const { return success_; }
        const std::vector<std::string>& messages() const { return messages_; }

    protected:
        std::vector<double> each_time;
        std::vector<bool> success_;
        std::vector<std::string> messages_;
        igl::Timer timer;
        std::string name_;
    };

    class Profiler {
    public:
        static Profiler& instance();
        std::shared_ptr<ProfilerPoint> create_point(std::string name);
        std::shared_ptr<ProfilerPoint> create_main_point(std::string name);
        void log(const std::string& fin = "");

    protected:
        Profiler() {}
        void write_summary(const std::string& dout, const std::string& fin);
        void write_point(const std::string& dout, const ProfilerPoint& point);
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

#define PROFILE_START(Name) _PROFILER_POINT_##Name->begin();
#define PROFILE_END(Name) _PROFILER_POINT_##Name->end();
#define PROFILE_SUCCESS(Name, Val) _PROFILER_POINT_##Name->success(Val);
#define PROFILE_MESSAGE(Name, Val) _PROFILER_POINT_##Name->message(Val);

#define LOG_PROFILER(SceneFile)                                                \
    ccd::profiler::Profiler::instance().log(SceneFile);

#else

#define PROFILE_MAIN_POINT(Description) ;
#define NAMED_PROFILE_POINT(Description, Name) ;
#define PROFILE_POINT(Description) ;
#define PROFILE_START(Name) ;
#define PROFILE_END(Name) ;
#define PROFILE_SUCCESS(Name, Val) ;
#define LOG_PROFILER(SceneFile) ;

#endif
