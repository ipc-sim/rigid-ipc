#include "constraint_factory.hpp"

#include <opt/barrier_constraint.hpp>
#include <opt/volume_constraint.hpp>

namespace ccd {
namespace opt {

    const ConstraintFactory& ConstraintFactory::factory()
    {
        static ConstraintFactory instance;

        return instance;
    }

    ConstraintFactory::ConstraintFactory()
    {
        constraints_.emplace("barrier_constraint", std::make_shared<BarrierConstraint>());
        constraints_.emplace("volume_constraint", std::make_shared<VolumeConstraint>());

        for (auto it = constraints_.begin(); it != constraints_.end(); ++it) {
            constraint_names_.push_back(it->first);
        }
    }

    std::shared_ptr<CollisionConstraint> ConstraintFactory::get_constraint(
        const std::string& constraint) const
    {
        auto it = constraints_.find(constraint);
        assert(it != constraints_.end());
        return it->second;
    }

} // namespace opt
} // namespace ccd
