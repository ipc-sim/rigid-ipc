#pragma once

#include <memory> // shared_ptr

#include <opt/collision_constraint.hpp>

namespace ccd {
namespace opt {
    class ConstraintFactory {
    public:
        static const ConstraintFactory& factory();

        std::shared_ptr<CollisionConstraint> get_constraint(const std::string& constraint) const;
        inline const std::vector<std::string>& get_constraint_names() const
        {
            return constraint_names_;
        }

    private:
        ConstraintFactory();
        std::map<std::string, std::shared_ptr<CollisionConstraint>> constraints_;
        std::vector<std::string> constraint_names_;
    };
} // namespace opt

} // namespace ccd
