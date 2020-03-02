#include <fenv.h>
#include <iomanip>
#include <iostream>

#include <boost/numeric/interval.hpp>

// clang-format off
typedef boost::numeric::interval<
    double,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_std<double> >,
        boost::numeric::interval_lib::checking_base<double> > >
    Interval;
// clang-format on

void show_fe_current_rounding_method(void)
{
    switch (fegetround()) {
    case FE_TONEAREST:
        std::cout << "FE_TONEAREST" << std::endl;
        break;
    case FE_DOWNWARD:
        std::cout << "FE_DOWNWARD" << std::endl;
        break;
    case FE_UPWARD:
        std::cout << "FE_UPWARD" << std::endl;
        break;
    case FE_TOWARDZERO:
        std::cout << "FE_TOWARDZERO" << std::endl;
        break;
    default:
        std::cout << "unknown" << std::endl;
    };
}

int main(int argc, char* argv[])
{
    {
        Interval theta = Interval(0.79358805865013693);
        Interval output = cos(theta);
        std::cout << std::setprecision(
                         std::numeric_limits<double>::digits10 + 1)
                  << "[" << output.lower() << ", " << output.upper() << "]"
                  << std::endl;
        std::cout << "is empty: "
                  << (output.lower() > output.upper() ? "true" : "false")
                  << std::endl;
    }

    {
        double theta = 0.79358805865013693;
        show_fe_current_rounding_method();
        volatile double output1 = cos(theta);
        std::cout << std::setprecision(
                         std::numeric_limits<double>::digits10 + 1)
                  << output1 << std::endl;
        fesetround(FE_DOWNWARD);
        show_fe_current_rounding_method();
        volatile double output2 = cos(theta);
        std::cout << std::setprecision(
                         std::numeric_limits<double>::digits10 + 1)
                  << output2 << std::endl;
        fesetround(FE_UPWARD);
        show_fe_current_rounding_method();
        volatile double output3 = cos(theta);
        std::cout << std::setprecision(
                         std::numeric_limits<double>::digits10 + 1)
                  << output3 << std::endl;
        fesetround(FE_TOWARDZERO);
        show_fe_current_rounding_method();
        volatile double output4 = cos(theta);
        std::cout << std::setprecision(
                         std::numeric_limits<double>::digits10 + 1)
                  << output4 << std::endl;
    }
}
