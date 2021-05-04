// clang-format off
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/operators.h>
// clang-format on

#include <boost/filesystem.hpp>
#include <tbb/global_control.h>
#include <tbb/task_scheduler_init.h>
#include <thread>

#include <SimState.hpp>
#include <physics/rigid_body.hpp>
#include <physics/rigid_body_problem.hpp>
#include <io/read_rb_scene.hpp>
#include <logger.hpp>
#include <profiler.hpp>

namespace py = pybind11;
using namespace ipc;
using namespace ipc::rigid;

// static tbb::global_control thread_limiter = tbb::global_control(
//     tbb::global_control::max_allowed_parallelism,
//     tbb::task_scheduler_init::default_num_threads());

PYBIND11_MODULE(rigidipc, m)
{
    m.doc() = "Rigid IPC";

    py::enum_<spdlog::level::level_enum>(m, "LoggerLevel")
        .value("trace", spdlog::level::level_enum::trace)
        .value("debug", spdlog::level::level_enum::debug)
        .value("info", spdlog::level::level_enum::info)
        .value("warn", spdlog::level::level_enum::warn)
        .value("error", spdlog::level::level_enum::err)
        .value("critical", spdlog::level::level_enum::critical)
        .value("off", spdlog::level::level_enum::off)
        .export_values();

    m.def(
        "set_logger_level",
        [](spdlog::level::level_enum level) {
            set_logger_level(static_cast<spdlog::level::level_enum>(level));
        },
        "Set log level", py::arg("level"));

    // m.def(
    //     "set_num_threads",
    //     [](int nthreads) {
    //         if (nthreads <= 0) {
    //             nthreads = tbb::task_scheduler_init::default_num_threads();
    //         } else if (
    //             nthreads > tbb::task_scheduler_init::default_num_threads()) {
    //             spdlog::warn(
    //                 "Attempting to use more threads than available ({:d} > "
    //                 "{:d})!",
    //                 nthreads,
    //                 tbb::task_scheduler_init::default_num_threads());
    //             nthreads = tbb::task_scheduler_init::default_num_threads();
    //         }
    //         thread_limiter = tbb::global_control(
    //             tbb::global_control::max_allowed_parallelism, nthreads);
    //     },
    //     "maximum number of threads to use", py::arg("nthreads"));

    m.def(
        "set_profiler_output_directory",
        [](const std::string& out_dir) { PROFILER_OUTDIR(out_dir); },
        "Set the output directory for the profiler (if enabled through CMake)",
        py::arg("out_dir"));

    py::class_<PoseD>(m, "Pose")
        .def(py::init<>())
        .def(
            py::init<const VectorMax3d&, const VectorMax3d&>(),
            "Pose of a rigid body", py::arg("position"), py::arg("rotation"))
        .def(
            "__repr__",
            [](const PoseD& p) {
                return fmt::format(
                    "Pose(position={}, rotation={})", fmt_eigen(p.position),
                    fmt_eigen(p.rotation));
            })
        .def_property(
            "position",
            [](PoseD& self) -> VectorMax3d& { return self.position; },
            [](PoseD& self, const VectorMax3d& position) {
                self.position = position;
            },
            py::return_value_policy::reference)
        .def_property(
            "rotation",
            [](PoseD& self) -> VectorMax3d& { return self.rotation; },
            [](PoseD& self, const VectorMax3d& rotation) {
                self.rotation = rotation;
            },
            py::return_value_policy::reference);

    py::enum_<RigidBodyType>(m, "RigidBodyType")
        .value("STATIC", RigidBodyType::STATIC)
        .value("KINEMATIC", RigidBodyType::KINEMATIC)
        .value("DYNAMIC", RigidBodyType::DYNAMIC)
        .export_values();

    py::class_<RigidBody>(m, "RigidBody")
        .def(py::init<
             const Eigen::MatrixXd&, const Eigen::MatrixXi&,
             const Eigen::MatrixXi&, const PoseD&, const PoseD&, const PoseD&,
             const double, const VectorMax6b&, const bool, const int>())
        .def(py::init([](const std::string& mesh_filename) {
            std::vector<nlohmann::json> jrbs = { { { "mesh",
                                                     mesh_filename } } };
            nlohmann::json rb_scene;
            rb_scene["rigid_bodies"] = jrbs;
            std::vector<RigidBody> rbs;
            read_rb_scene(rb_scene, rbs);
            return rbs[0];
        }))
        .def_readwrite("name", &RigidBody::name)
        .def_readwrite("group_id", &RigidBody::group_id)
        .def_readwrite("type", &RigidBody::type)
        .def_readonly("vertices", &RigidBody::vertices)
        .def(
            "world_vertices",
            [](const RigidBody& self) { return self.world_vertices(); })
        .def_readonly("edges", &RigidBody::edges)
        .def_readonly("faces", &RigidBody::faces)
        .def_readwrite("pose", &RigidBody::pose)
        .def_readwrite("kinematic_poses", &RigidBody::kinematic_poses);

    py::class_<RigidBodyAssembler>(m, "RigidBodyAssembler")
        .def(
            "__getitem__",
            [](RigidBodyAssembler& self, size_t i) -> RigidBody& {
                return self[i];
            },
            py::return_value_policy::reference)
        .def(
            "__len__",
            [](const RigidBodyAssembler& self) { return self.num_bodies(); })
        .def(
            "__iter__",
            [](RigidBodyAssembler& self) {
                return py::make_iterator(self.m_rbs.begin(), self.m_rbs.end());
            },
            // Essential: keep object alive while iterator exists
            py::keep_alive<0, 1>(), py::return_value_policy::reference);

    py::class_<SimState>(m, "Simulation")
        .def(py::init<>())
        .def(
            "load_scene", &SimState::load_scene,
            "Load a simulation scene from a file (JSON).\n"
            "Optionally provide a JSON to patch the file.",
            py::arg("filename"), py::arg("patch") = "")
        .def(
            "bodies",
            [](const SimState& sim) -> RigidBodyAssembler& {
                return std::dynamic_pointer_cast<RigidBodyProblem>(
                           sim.problem_ptr)
                    ->m_assembler;
            },
            py::return_value_policy::reference)
        .def("step", &SimState::simulation_step, "Take a single step")
        .def(
            "run", &SimState::run_simulation, "Run the entire simulation",
            py::arg("fout"))
        .def(
            "save_obj_sequence", &SimState::save_obj_sequence,
            "Save the simulation as a sequence of OBJ files",
            py::arg("dir_name"))
        .def(
            "save_gltf", &SimState::save_gltf,
            "Save the simulation as a GLTF animation file", py::arg("filename"))
        .def(
            "save_simulation", &SimState::save_simulation,
            "Save the simulation as a JSON file", py::arg("filename"))
        .def_readwrite(
            "max_simulation_steps", &SimState::m_max_simulation_steps)
        .def_readwrite(
            "checkpoint_frequency", &SimState::m_checkpoint_frequency)
        .def_property(
            "timestep",
            [](const SimState& self) { return self.problem_ptr->timestep(); },
            [](const SimState& self, double timestep) {
                self.problem_ptr->timestep(timestep);
            },
            "Time step size");
}
