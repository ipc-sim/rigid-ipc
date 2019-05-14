#pragma once
#include <opt/volume_constraint.hpp>
#include <opt/barrier_constraint.hpp>

namespace ccd {
void volume_constraint_menu(ccd::opt::VolumeConstraint& cstr);
void barrier_constraint_menu(ccd::opt::BarrierConstraint& cstr);
}
