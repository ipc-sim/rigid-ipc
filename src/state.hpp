#ifndef CCD_STATE_HPP
#define CCD_STATE_HPP

#include <Eigen/Core>
#include <set>

#include <ccd/collision_detection.hpp>
#include <ccd/impact.hpp>
#include <ccd/prune_impacts.hpp>

#include <solvers/barrier_solver.hpp>
#include <solvers/ipopt_solver.hpp>
#include <solvers/ncp_solver.hpp>
#include <solvers/nlopt_solver.hpp>
#include <solvers/optimization_solver.hpp>
#include <solvers/qp_solver.hpp>

#include <opt/barrier_constraint.hpp>
#include <opt/collision_constraint.hpp>
#include <opt/particles_problem.hpp>
#include <opt/rigid_body_problem.hpp>
#include <opt/volume_constraint.hpp>

#include <physics/rigid_body_system.hpp>

namespace ccd {

enum class ConstraintType { VOLUME, BARRIER };
static const char* ConstraintNames[2] = { "VOLUME", "BARRIER" };

enum class OptimizationMethod {
    NLOPT,                  ///< @brief Container for various algorithms.
    IPOPT,                  ///< @brief Interior-Point Method (Ipopt).
    LINEARIZED_CONSTRAINTS, ///< @brief Linearize the constraints and solve the
                            ///< QP (OSQP/MOSEK).
    NCP,                    ///< @brief Nonlinear Complementarity Problem.
    BARRIER_SOLVER          ///< @brief Apply constraints as barriers.
};

static const char* OptimizationMethodNames[]
    = { "NLOPT", "IPOPT", "Linearized Const.", "NCP", "Barrier Solver" };
/**
 * @brief The State class keeps the full state of the UI and the collisions.
 */
class State {
public:
    static const int kDIM = 2;

    State();
    virtual ~State() = default;

    /// @brief #V,2 vertices positions
    Eigen::MatrixX2d vertices;
    /// @brief #E,2 vertices connnectivity
    Eigen::MatrixX2i edges;
    /// @brief #V,2 vertices displacements
    Eigen::MatrixX2d displacements;

    ccd::physics::RigidBodySystem rigid_body_system;
    bool is_rigid_bodies_mode;

    /// @brief The current number of pruned impacts
    int num_pruned_impacts;

    /// @brief #E,1 contact volume for each edge
    Eigen::VectorXd volumes;

    /// @brief #E,2*V contact gradient for each edge
    Eigen::MatrixXd volume_grad;

    std::string output_dir;

    ////////////////////////////////////////////////////////////////////////////
    // Optimization Fields

    opt::OptimizationResults opt_results;

    OptimizationMethod opt_method;

    opt::NCPSolver ncp_solver;
    opt::IpoptSolver ipopt_solver;
    opt::NLOptSolver nlopt_solver;
    opt::QPSolver qp_solver;
    opt::BarrierSolver barrier_solver;

    opt::VolumeConstraint volume_constraint;
    opt::BarrierConstraint barrier_constraint;
    ConstraintType constraint_function;

    opt::ParticlesDisplProblem particles_problem;
    opt::RigidBodyProblem rigid_body_problem;

    /// @brief reuse the current opt_displacements to initialize the
    /// optimization
    bool reuse_opt_displacements = false;

    ///@brief Optimization step history for displacements
    std::vector<Eigen::MatrixX2d> u_history;

    ///@brief Optimization step history for functional
    std::vector<double> f_history;

    ///@brief Optimization step history for constraints
    std::vector<double> gsum_history;
    std::vector<Eigen::VectorXd> g_history;
    std::vector<Eigen::MatrixXd> jac_g_history;
    std::vector<std::vector<long>> collision_history;

    ////////////////////////////////////////////////////////////////////////////
    // SCENE CRUD
    // ----------------------------------------------------------------------
    void load_scene(const std::string filename);
    void load_scene(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges, const Eigen::MatrixX2d& displacements);
    void save_scene(const std::string filename);
    void reset_scene();
    void fit_scene_to_canvas();

    void add_vertex(const Eigen::RowVector2d& vertex);
    void add_edges(const Eigen::MatrixX2i& edges);
    void duplicate_selected_vertices(Eigen::Vector2d delta_center_of_mass);
    void expand_is_dof_fixed();
    bool remove_free_vertices();

    void set_vertex_position(
        const int vertex_idx, const Eigen::RowVector2d& position);
    void move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta);
    void move_displacement(
        const int vertex_idx, const Eigen::RowVector2d& delta);

    ////////////////////////////////////////////////////////////////////////////
    // SCENE CCD
    // ----------------------------------------------------------------------
    void reset_scene_data();
    void run_ccd_pipeline();

    ////////////////////////////////////////////////////////////////////////////
    // SCENE OPT
    // ----------------------------------------------------------------------
    opt::CollisionConstraint& getCollisionConstraint();
    opt::OptimizationSolver& getOptimizationSolver();
    opt::OptimizationProblem& getOptimizationProblem();

    void reset_optimization_problem();
    void optimize_displacements(const std::string filename = "");

    void load_optimization(const std::string filename);
    void save_optimization(const std::string filename);
    bool record_optimization_step(
        const Eigen::VectorXd& x, const Eigen::MatrixX2d& Uk);
    void post_process_optimization();

    ////////////////////////////////////////////////////////////////////////////
    // UI
    // ----------------------------------------------------------------------
    Eigen::MatrixX2d get_vertex_at_time();
    Eigen::MatrixX2d get_volume_grad();
    // opt results
    double get_opt_functional();
    Eigen::MatrixX2d get_opt_displacements();
    Eigen::MatrixX2d get_opt_vertex_at_time();
    Eigen::MatrixX2d get_opt_volume_grad();
    Eigen::VectorXd get_opt_volume();
    std::vector<long> get_opt_collision_edges();

    /// @brief setup log
    int log_level;

    /// @brief Background rectangle to detect clicks
    double canvas_width, canvas_height;

    /// @brief We show the scene at time=`time` between 0 and 1

    float current_time;
    /// @brief Current user-selection of vertex and displacement points
    std::vector<int> selected_points, selected_displacements;

    /// @brief Use for any functionallity that requires showing info of only
    /// one volume
    int current_volume;

    /// when going to next edge, skip edges with no impact
    bool skip_no_impact_edge;

    /// scaling for drawing the gradient
    float grad_scaling;

    /// if true, gradient displayed comes from opt data.
    /// Else, it comes from the volume generated by the user-displacements.
    bool use_opt_gradient;

    // UI OPT
    // ----------------------------------------------------------------------
    ///@brief Time along the optimal displacements
    float current_opt_time;

    ///@brief we show the values of this iteration
    int current_opt_iteration;
    void reset_results();

    std::vector<std::vector<int>> create_adjacency_list();
    void find_connected_vertices(const long& vertex_id,
        const std::vector<std::vector<int>>& adjacency_list,
        std::unordered_set<int>& connected_vertices);

    void convert_connected_components_to_rigid_bodies();
    void update_fields_from_rigid_bodies();
    void update_displacements_from_rigid_bodies();
};

} // namespace ccd
#endif
