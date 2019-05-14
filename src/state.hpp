#ifndef CCD_STATE_HPP
#define CCD_STATE_HPP

#include <Eigen/Core>

#include <ccd/collision_detection.hpp>
#include <ccd/impact.hpp>
#include <ccd/prune_impacts.hpp>

#include <opt/ipopt_solver.hpp>
#include <opt/solver.hpp>

#include <opt/barrier_constraint.hpp>
#include <opt/displacement_opt.hpp>
#include <opt/volume_constraint.hpp>
#include <opt/collision_constraint.hpp>

namespace ccd {

enum class ConstraintType { VOLUME, BARRIER };
static const char* ConstraintNames[2] = { "VOLUME", "BARRIER" };

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

    /// @brief All edge-vertex contact
    EdgeVertexImpacts ev_impacts;

    /// @brief All edge-edge contact
    EdgeEdgeImpacts ee_impacts;

    /// @brief #E,1 indices of the edges' first impact
    Eigen::VectorXi edge_impact_map;

    /// @brief The current number of pruned impacts
    int num_pruned_impacts;

    /// @brief #E,1 contact volume for each edge
    Eigen::VectorXd volumes;

    /// @brief #E,2*V contact gradient for each edge
    Eigen::MatrixXd volume_grad;

    /// @brief method to use for contact detection
    DetectionMethod detection_method;

    /// @brief epsilon use on volume computation
    double volume_epsilon;

    std::string output_dir;

    ////////////////////////////////////////////////////////////////////////////
    // Optimization Fields

    /// @brief Optimization problem to solve
//    opt::OptimizationProblem opt_problem;

    /// @brief #V,2 optimized vertices displacements
    opt::OptimizationResults opt_results;

    /// @brief Settings for the problem solver
    opt::SolverSettings solver_settings;

    opt::NCPDisplacementOptimization ncp_displ_solver;
    opt::IpoptSolver ipopt_solver;

    opt::VolumeConstraint volume_constraint;
    opt::BarrierConstraint barrier_constraint;
    ConstraintType constraint_function;

    opt::ParticlesDisplProblem opt_problem;

    /// @brief if True, reuse the current opt_displacements for initial
    /// optimization
    bool reuse_opt_displacements = false;

    /// @brief if True, recompute collision set on each evaluation of the
    /// collision volume and gradient
    bool recompute_collision_set = false;

    ///@brief Optimization step history for displacements
    std::vector<Eigen::MatrixX2d> u_history;

    ///@brief Optimization step history for functional
    std::vector<double> f_history;

    ///@brief Optimization step history for constraints
    std::vector<double> gsum_history;
    std::vector<Eigen::VectorXd> g_history;

    std::vector<Eigen::MatrixXd> jac_g_history;

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

    void set_vertex_position(
        const int vertex_idx, const Eigen::RowVector2d& position);
    void move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta);
    void move_displacement(
        const int vertex_idx, const Eigen::RowVector2d& delta);

    ////////////////////////////////////////////////////////////////////////////
    // SCENE CCD
    // ----------------------------------------------------------------------
    void reset_impacts();
    void run_ccd_pipeline();

    ////////////////////////////////////////////////////////////////////////////
    // SCENE OPT
    // ----------------------------------------------------------------------
    opt::CollisionConstraint& getCollisionConstraint();
    void reset_optimization_problem();
    void optimize_displacements(const std::string filename = "");

    void load_optimization(const std::string filename);
    void save_optimization(const std::string filename);

    ////////////////////////////////////////////////////////////////////////////
    // SCENE OUTPUT
    // ----------------------------------------------------------------------
    void log_optimization_steps(const std::string filename,
        std::vector<Eigen::VectorXd>& it_x,
        std::vector<Eigen::VectorXd>& it_lambda, std::vector<double>& it_gamma);

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

    /// @brief Background rectangle to detect clicks
    double canvas_width, canvas_height;

    /// @brief We show the scene at time=`time` between 0 and 1

    float current_time;
    /// @brief Current user-selection of vertex and displacement points
    std::vector<int> selected_points, selected_displacements;

    /// @brief Use for any functionallity that requires showing only one ev
    /// impact
    int current_ev_impact;

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
    ///@brief Time along the optimal displacments
    float current_opt_time;

    ///@brief we show the values of this iteration
    int current_opt_iteration;
    void reset_results();
};

} // namespace ccd
#endif
