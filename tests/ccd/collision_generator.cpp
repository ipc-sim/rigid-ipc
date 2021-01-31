#include "collision_generator.hpp"

#include <iomanip>
#include <stdlib.h> /* rand */

namespace ipc::rigid {
namespace unittests {

    std::ostream& operator<<(std::ostream& o, const TestImpact& p)
    {
        // clang-format off
        return o << std::setprecision(10)
                 << " Vi << " << p.Vi(0) << ", " << p.Vi(1) << ";\n"
                 << " Vj << " << p.Vj(0) << ", " << p.Vj(1) << ";\n"
                 << " Vk << " << p.Vk(0) << ", " << p.Vk(1) << ";\n"
                 << " Ui << " << p.Ui(0) << ", " << p.Ui(1) << ";\n"
                 << " Uj << " << p.Uj(0) << ", " << p.Uj(1) << ";\n"
                 << " Uk << " << p.Uj(0) << ", " << p.Uk(1) << ";\n"
                 << " toi = " << p.toi << ";\n"
                 << " alpha = " << p.alpha  << ";\n"
                 << std::endl;
        // clang-format on
    }

    TestImpactGenerator::TestImpactGenerator(size_t value, bool rigid)
        : rigid(rigid)
        , current_i(0)
        , max_i(value)
    {
        current = generate_random_impact(rigid);
    }

    bool TestImpactGenerator::next()
    {
        current = generate_random_impact(rigid);
        return ++current_i < max_i;
    }
    TestImpact const& TestImpactGenerator::get() const { return current; }

    Catch::Generators::GeneratorWrapper<TestImpact>
    random_impacts(size_t value, bool rigid)
    {
        return Catch::Generators::GeneratorWrapper<TestImpact>(
            std::make_unique<TestImpactGenerator>(value, rigid));
    }

    TestImpact generate_random_impact(const bool rigid)
    {
        static const double kEPSILON = 1E-8;
        static const double kSCALE = 10.0;

        TestImpact impact;

        // set first edge random positions keeping in larger than kEPSILON
        impact.Vi.setRandom(2);
        impact.Vj.setRandom(2);
        impact.Vi *= kSCALE;
        impact.Vk *= kSCALE;

        double len_ij = (impact.Vj - impact.Vi).norm();
        if (len_ij < kEPSILON) {
            impact.Vj += (impact.Vj - impact.Vi) / len_ij * kEPSILON;
            assert((impact.Vj - impact.Vi).norm() >= kEPSILON);
        }

        // set first edge random i.Velocities keeping in larger than kEPSILON
        // note that this doest NOT ensure its length is > 0 oi.Ver the full
        // time.
        impact.Ui.setRandom(2);
        impact.Ui *= kSCALE;
        if (rigid) {
            impact.Uj = impact.Ui;
        } else {
            impact.Uj.setRandom();
            len_ij = (impact.Vj + impact.Uj - impact.Vi - impact.Ui).norm();
            if (len_ij < kEPSILON) {
                impact.Uj += (impact.Vj + impact.Uj - impact.Vi - impact.Ui)
                    / len_ij * kEPSILON;
                assert(
                    (impact.Vj + impact.Uj - impact.Vi - impact.Ui).norm()
                    >= kEPSILON);
            }
        }

        // Generate a time and position for the impact.
        // How do we ensure it is the first time of impact?
        double toi = double(rand()) / RAND_MAX;
        double alpha = double(rand()) / RAND_MAX;
        Eigen::Vector2d Ve, Vi_toi, Vj_toi;
        Vi_toi = impact.Vi + impact.Ui * toi;
        Vj_toi = impact.Vj + impact.Uj * toi;
        Ve = Vi_toi + alpha * (Vj_toi - Vi_toi);

        // Ve = Vk + Uk * toi --> (Ve - Vk)/toi = Uk
        impact.Vk.setRandom(2);
        impact.Vk *= kSCALE;
        impact.Uk = (Ve - impact.Vk) / toi;

        impact.Vl.setRandom(2);
        impact.Vl *= kSCALE;
        double len_kl = (impact.Vl - impact.Vk).norm();
        if (len_kl < kEPSILON) {
            impact.Vl += (impact.Vl - impact.Vk) / len_kl * kEPSILON;
            assert((impact.Vl - impact.Vk).norm() >= kEPSILON);
        }

        if (rigid) {
            impact.Ul = impact.Uk;
        } else {
            impact.Ul.setRandom();
            len_kl = (impact.Vl + impact.Ul - impact.Vk - impact.Uk).norm();
            if (len_kl < kEPSILON) {
                impact.Ul += (impact.Vl + impact.Ul - impact.Vk - impact.Uk)
                    / len_kl * kEPSILON;
                assert(
                    (impact.Vl + impact.Ul - impact.Vk - impact.Uk).norm()
                    >= kEPSILON);
            }
        }

        impact.toi = toi;
        impact.alpha = alpha;

        return impact;
    }
} // namespace unittests
} // namespace ipc::rigid
