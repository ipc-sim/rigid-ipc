
/******************************************************************************
*                    Code generated with jinja and sympy                     *
*                                                                            *
*                      See src/autogen/.py for originals                     *
*                                                                            * 
*                            DO NOT MODIFY THIS FILE                         *
******************************************************************************/
#ifndef CCD_AUTO_COLLISION_VOLUME_HPP
#define CCD_AUTO_COLLISION_VOLUME_HPP

#include <Eigen/Core>

namespace ccd {
namespace autogen {

    double collision_volume(const Eigen::Vector2d& Vi,
        const Eigen::Vector2d& Vj,
        const Eigen::Vector2d& /*Vk*/,
        const Eigen::Vector2d& /*Vl*/,
        const Eigen::Vector2d& Ui,
        const Eigen::Vector2d& Uj,
        const Eigen::Vector2d& /*Uk*/,
        const Eigen::Vector2d& /*Ul*/,
        const double epsilon);

}
}
#endif
