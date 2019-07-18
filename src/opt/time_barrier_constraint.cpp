#include "time_barrier_constraint.hpp"

#include <opt/barrier.hpp>

#include <profiler.hpp>

namespace ccd {
namespace opt {

    NLOHMANN_JSON_SERIALIZE_ENUM(InitialBarrierEpsilon,
        { { MIN_TOI, "min_toi" }, { MAX_TOI, "max_toi" },
            { CUSTOM, "custom" } })

    TimeBarrierConstraint::TimeBarrierConstraint()
        : TimeBarrierConstraint("time_barrier_constraint")
    {
    }

    TimeBarrierConstraint::TimeBarrierConstraint(const std::string& name)
        : CollisionConstraint(name)
        , barrier_epsilon(0.0)
        , custom_inital_epsilon(1.0)
        , initial_epsilon(InitialBarrierEpsilon::MIN_TOI)

    {
    }

    void TimeBarrierConstraint::settings(const nlohmann::json& json)
    {
        initial_epsilon = json["initial_epsilon"].get<InitialBarrierEpsilon>();
        custom_inital_epsilon = json["custom_initial_epsilon"].get<double>();
        CollisionConstraint::settings(json);
    }

    nlohmann::json TimeBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["initial_epsilon"] = initial_epsilon;
        json["custom_inital_epsilon"] = custom_inital_epsilon;

        return json;
    }

    void TimeBarrierConstraint::initialize(const Eigen::MatrixX2d& vertices,
        const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& Uk)
    {
        this->barrier_epsilon = 0.0; // Temporary until collisions are computed.
        CollisionConstraint::initialize(vertices, edges, Uk);
        resetBarrierEpsilon();
    }

    void TimeBarrierConstraint::detectCollisions(const Eigen::MatrixXd& Uk)
    {
        const double time_scale = 1 + 2 * this->barrier_epsilon;
        ccd::detect_edge_vertex_collisions(*vertices, time_scale * Uk, *edges,
            ev_impacts, detection_method, true);
        for (EdgeVertexImpact& ev_impact : ev_impacts) {
            ev_impact.time *= time_scale;
        }
        ccd::convert_edge_vertex_to_edge_edge_impacts(
            *edges, ev_impacts, ee_impacts);
    }

    int TimeBarrierConstraint::number_of_constraints()
    {
        return int(ee_impacts.size() * 2);
    }

    void TimeBarrierConstraint::resetBarrierEpsilon()
    {
        // Assumes the collisions have already been detected
        if (this->ee_impacts.size() > 0) {
            switch (initial_epsilon) {
            case InitialBarrierEpsilon::MIN_TOI:
                this->barrier_epsilon = std::min_element(
                    this->ee_impacts.begin(), this->ee_impacts.end(),
                    compare_impacts_by_time<EdgeEdgeImpact>)
                                            ->time;
                break;
            case InitialBarrierEpsilon::MAX_TOI:
                this->barrier_epsilon = std::max_element(
                    this->ee_impacts.begin(), this->ee_impacts.end(),
                    compare_impacts_by_time<EdgeEdgeImpact>)
                                            ->time;
                break;
            case InitialBarrierEpsilon::CUSTOM:
                this->barrier_epsilon = custom_inital_epsilon;
                break;
            }
        } else {
            this->barrier_epsilon = 0;
        }
    }

    void TimeBarrierConstraint::compute_constraints(
        const Eigen::MatrixXd& Uk, Eigen::VectorXd& barriers)
    {
        std::vector<double> v_barriers;
        compute_constraints_per_impact(Uk, v_barriers);
        barriers = Eigen::Map<Eigen::VectorXd>(
            v_barriers.data(), int(v_barriers.size()));
    }

    void TimeBarrierConstraint::compute_constraints_jacobian(

        const Eigen::MatrixXd& Uk, Eigen::MatrixXd& barriers_jacobian)
    {
        std::vector<DScalar> barriers;
        compute_constraints_per_impact(Uk, barriers);
        assemble_jacobian(barriers, barriers_jacobian);
    }

    void TimeBarrierConstraint::compute_constraints_hessian(
        const Eigen::MatrixXd& Uk,
        std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)
    {
        std::vector<DScalar> barriers;
        compute_constraints_per_impact(Uk, barriers);
        assemble_hessian(barriers, barriers_hessian);
    }

    void TimeBarrierConstraint::compute_constraints_and_derivatives(
        const Eigen::MatrixXd& Uk,
        Eigen::VectorXd& barriers,
        Eigen::MatrixXd& barriers_jacobian,
        std::vector<Eigen::SparseMatrix<double>>& barriers_hessian)

    {
        std::vector<DScalar> v_barriers;
        compute_constraints_per_impact(Uk, v_barriers);
        PROFILE(assemble_constraints(v_barriers, barriers),
            ProfiledPoint::COMPUTING_CONSTRAINTS)
        PROFILE(assemble_jacobian(v_barriers, barriers_jacobian),
            ProfiledPoint::COMPUTING_GRADIENT)
        PROFILE(assemble_hessian(v_barriers, barriers_hessian),
            ProfiledPoint::COMPUTING_HESSIAN)
    }

    template <typename T>
    void TimeBarrierConstraint::compute_constraints_per_impact(
        const Eigen::MatrixXd& displacements, std::vector<T>& constraints)
    {
        constraints.clear();
        constraints.reserve(2 * ee_impacts.size());

        if (differentiable<T>()) {
            // !!! All definitions using DScalar must be done after this !!!
            DiffScalarBase::setVariableCount(8);
        }

        // -----------------------------------------------------------
        // TODO: fix 2x trick
        T t = T(1);                             // Final time for the time step
        T time_scale = t + 2 * barrier_epsilon; // 2x for safety
        // -----------------------------------------------------------

        for (size_t i = 0; i < ee_impacts.size(); ++i) {
            ImpactTData<T> data
                = get_impact_data<T>(displacements, ee_impacts[i]);

            // -----------------------------------------------------------
            // TODO: fix this trick
            for (uint ii = 0; ii < 4; ii++) {
                data.u[ii] *= time_scale;
            }
            // -----------------------------------------------------------

            T toi, alpha_ij, alpha_kl;
            T vol_ij(0), vol_kl(0);

            if (compute_toi_alpha(data, toi, alpha_ij, alpha_kl)) {
                T barrier = opt::spline_barrier<T>(
                    time_scale * toi - t, barrier_epsilon);

                vol_ij = (data.v[0] - data.v[1]).norm() * barrier;
                vol_kl = (data.v[2] - data.v[3]).norm() * barrier;
            }
            constraints.push_back(vol_ij);
            constraints.push_back(vol_kl);
        }
    }

    template void TimeBarrierConstraint::compute_constraints_per_impact<double>(
        const Eigen::MatrixXd& displacements, std::vector<double>& constraints);

    template void TimeBarrierConstraint::compute_constraints_per_impact<DScalar>(
        const Eigen::MatrixXd& displacements,
        std::vector<DScalar>& constraints);
} // namespace opt
} // namespace ccd
