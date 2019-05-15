#include "InterferenceVolumeV2.h"

#include <algorithm> // std::swap
#include <iostream>
#include <vector>

#include "rpoly.h" // root finder

namespace IAGM {

bool time_of_collision(
    const Eigen::Vector2d x00, const Eigen::Vector2d x01,
    const Eigen::Vector2d x10, const Eigen::Vector2d x11,
    const Eigen::Vector2d x20, const Eigen::Vector2d x21,
    const double min_sep, const double tolerance,
    double& toi, double& alpha)
{

    using namespace Eigen;
    Vector2d x1(x00[0], x00[1]);
    Vector2d x2(x10[0], x10[1]);
    Vector2d x3(x20[0], x20[1]);

    Vector2d v1(x01[0] - x00[0], x01[1] - x00[1]);
    Vector2d v2(x11[0] - x10[0], x11[1] - x10[1]);
    Vector2d v3(x21[0] - x20[0], x21[1] - x20[1]);

    Vector2d x1x2 = x1 - x2;
    Vector2d v1v2 = v1 - v2;
    Vector2d x3x2 = x3 - x2;
    Vector2d v3v2 = v3 - v2;

    std::swap(x3x2[0], x3x2[1]);
    x3x2[1] *= -1.0;
    std::swap(v3v2[0], v3v2[1]);
    v3v2[1] *= -1.0;

    double a0 = x1x2.dot(x3x2);
    double a1 = x1x2.dot(v3v2) + x3x2.dot(v1v2);
    double a2 = v1v2.dot(v3v2);

    std::vector<double> coeffs;
    if (min_sep == 0.0) {
        coeffs.resize(3);
        coeffs[0] = a2;
        coeffs[1] = a1;
        coeffs[2] = a0;
    } else {
        double h = min_sep;
        coeffs.resize(5);

        coeffs[0] = a2 * a2;
        coeffs[1] = 2.0 * a1 * a2;
        coeffs[2] = a1 * a1 + 2.0 * a0 * a2;
        coeffs[3] = 2.0 * a0 * a1;
        coeffs[4] = a0 * a0;

        coeffs[2] -= h * h * v3v2.dot(v3v2);
        coeffs[3] -= h * h * 2.0 * x3x2.dot(v3v2);
        coeffs[4] -= h * h * x3x2.dot(x3x2);
    }

    // TODO: Coefficient pruning tests could speed this up
    //
    int degree = coeffs.size() - 1;
    while (coeffs[coeffs.size() - 1 - degree] == 0.0 && degree > 0)
        degree--;

    int cnt = 0;
    double rl[4], im[4];
    RootFinder rf;
    if (degree == 0) {
        cnt = 0;
    } else {
        // Find roots of polynomial
        //
        cnt = rf.rpoly(&coeffs[coeffs.size() - 1 - degree], degree, rl, im);
        /*
          std::cout<<"roots count: "<<cnt<<std::endl;
          std::cout<<"coeffs 0: "<<coeffs[0]<<std::endl;
          std::cout<<"coeffs 1: "<<coeffs[1]<<std::endl;
          std::cout<<"coeffs 2: "<<coeffs[2]<<std::endl;
          std::cout<<"coeffs 3: "<<coeffs[3]<<std::endl;
          std::cout<<"coeffs 4: "<<coeffs[4]<<std::endl;
          */
    }

    // Try vertex-segment two half circles collision
    Vector2d x1x3 = x1 - x3;
    Vector2d v1v3 = v1 - v3;
    double minSep = min_sep;

    // Test x1 collides wtih half circle of x3
    double coeffs_circle1[3];
    coeffs_circle1[0] = v1v3.dot(v1v3);
    coeffs_circle1[1] = 2 * x1x3.dot(v1v3);
    coeffs_circle1[2] = x1x3.dot(x1x3) - minSep * minSep;
    degree = 2;
    if (coeffs_circle1[0] == 0.0)
        degree = 1;
    double rl_circle1[2], im_circle1[2];
    int cnt_circle1 = rf.rpoly(&coeffs_circle1[2 - degree], degree, rl_circle1, im_circle1);

    // Test x1 collides with half circle of x2
    double coeffs_circle2[3];
    coeffs_circle2[0] = v1v2.dot(v1v2);
    coeffs_circle2[1] = 2 * x1x2.dot(v1v2);
    coeffs_circle2[2] = x1x2.dot(x1x2) - minSep * minSep;
    degree = 2;
    if (coeffs_circle2[0] == 0.0)
        degree = 1;
    double rl_circle2[2], im_circle2[2];
    int cnt_circle2 = rf.rpoly(&coeffs_circle2[2 - degree], degree, rl_circle2, im_circle2);

    //cnt_circle1 = 0;
    //cnt_circle2 = 0;
    int cnt_total = cnt + cnt_circle1 + cnt_circle2;
    double rl_total[8], im_total[8];
    int ctype[8];
    int root_count = 0;
    for (int i = 0; i < cnt; ++i) {
        rl_total[root_count] = rl[i];
        im_total[root_count] = im[i];
        ctype[root_count] = 0;
        root_count = root_count + 1;
    }

    for (int i = 0; i < cnt_circle1; ++i) {
        rl_total[root_count] = rl_circle1[i];
        im_total[root_count] = im_circle1[i];
        ctype[root_count] = 2;
        root_count = root_count + 1;
    }

    for (int i = 0; i < cnt_circle2; ++i) {
        rl_total[root_count] = rl_circle2[i];
        im_total[root_count] = im_circle2[i];
        ctype[root_count] = 1;
        root_count = root_count + 1;
    }

    double min_t = 10000;
    double min_s = 10000;
    double min_d = 10000;
    double min_ct = 3;

    Vector2d pe = x1 + v1;
    Vector2d ae = x2 + v2;
    Vector2d be = x3 + v3;

    Vector2d abe = be - ae;
    Vector2d ape = pe - ae;

    double se = abe.dot(ape) / abe.dot(abe);
    if (se < 0.00)
        se = 0.0;
    if (se > 1.00)
        se = 1.0;
    double de = (pe - (ae + abe * se)).squaredNorm() - minSep * minSep;

    Vector2d ps = x1;
    Vector2d as = x2;
    Vector2d bs = x3;

    Vector2d abs = bs - as;
    Vector2d aps = ps - as;

    double ss = abs.dot(aps) / abs.dot(abs);
    if (ss < 0.00)
        ss = 0.0;
    if (ss > 1.00)
        ss = 1.0;
    double ds = (ps - (as + abs * ss)).squaredNorm() - minSep * minSep;
    //assert(ds > -1e-10);
    if (ds < -1e-10)
        std::cout << "ds:" << ds << std::endl;
    //std::cout<<"begin of one p, a, b intersection."<<min_d<<std::endl;
    for (int i = 0; i < cnt_total; ++i) {
        //if (im_total[i] == 0.0 && rl_total[i] >= 0.0 && rl_total[i] <= 1.0)
        if ((std::abs(im_total[i]) < 1e-16) && (rl_total[i] >= -1e-16) && (rl_total[i] <= (1.0 + 1e-16))) {
            toi = rl_total[i];
            if (toi < 0)
                toi = 0;
            if (toi > 1)
                toi = 1;
            //std::cout<<"total roots: "<<cnt_total<<", i: "<<i<<std::endl;
            //std::cout<<"p: ["<<x01[0]<<","<<x01[1]<<"]"<<std::endl;
            //std::cout<<"a: ["<<x11[0]<<","<<x11[1]<<"]"<<std::endl;
            //std::cout<<"b: ["<<x21[0]<<","<<x21[1]<<"]"<<std::endl;

            // Check rl[i] for intersection
            //
            Vector2d p = x1 + v1 * toi;
            Vector2d a = x2 + v2 * toi;
            Vector2d b = x3 + v3 * toi;

            Vector2d ab = b - a;
            Vector2d ap = p - a;

            alpha = ab.dot(ap) / ab.dot(ab);

            //std::cout<<"s: "<<s<<std::endl;
            // Clamp to edge endpoints
            //
            if (alpha < 0.00)
                alpha = 0.0;
            if (alpha > 1.00)
                alpha = 1.0;

            // NOTE: If the code misses intersections, it's probably because of
            // this tolerance. There is a Eurographics 2012 paper which attempts to
            // address this, but for now it can be tricky to choose a reliable epsilon.
            //
            double d = (p - (a + ab * alpha)).squaredNorm() - minSep * minSep;
            //std::cout<<"d: "<<d<<std::endl;
            //std::cout<<"t: "<<t<<std::endl;

            // safeguard test if t is close to 1, check the ending position
            if ((d < tolerance) || (de <= 0))
            //if ((d < C->m_tol))
            {
                // safeguard test if t is close to 0, check the start position

                // To ensure a conservative interference volume, slightly enlarge
                // it by returning an earlier time
                /*
                std::cout<<"d smaller than tol "<<std::endl;
                std::cout<<"get intersection of p, a, b start config: "<<std::endl;
                std::cout<<"p: ["<<x00[0]<<","<<x00[1]<<"]"<<std::endl;
                std::cout<<"a: ["<<x10[0]<<","<<x10[1]<<"]"<<std::endl;
                std::cout<<"b: ["<<x20[0]<<","<<x20[1]<<"]"<<std::endl;
                std::cout<<"get intersection of p, a, b end config: "<<std::endl;
                std::cout<<"p: ["<<x01[0]<<","<<x01[1]<<"]"<<std::endl;
                std::cout<<"a: ["<<x11[0]<<","<<x11[1]<<"]"<<std::endl;
                std::cout<<"b: ["<<x21[0]<<","<<x21[1]<<"]"<<std::endl;
                std::cout<<"d: "<<d<<std::endl;
                std::cout<<"t: "<<t<<std::endl;
                std::cout<<"s: "<<s<<std::endl;
                std::cout<<"p to a: "<<(p-a).norm()<<std::endl;
                std::cout<<"p to b: "<<(p-b).norm()<<std::endl;
                */
                //t = max(0.0, t - 0);
                if (toi < min_t) {
                    //std::cout<<"smaller time accepted, i: "<<i<<std::endl;
                    min_t = toi;
                    min_s = alpha;
                    min_d = d;
                    min_ct = ctype[i];
                }
                //t = max(0.0, t-0.01);
                //return true;
            }
        }
    }

    if (min_t >= 0.0 && min_t <= 1.0) {
        //std::cout<<"End of one a,b,p intersection, end of d: "<<min_d<<std::endl;
        //std::cout<<"Smallest collision time:"<<min_t<<std::endl;
        //t = min_t;
        toi = std::max(0.0, min_t - 0);
        alpha = min_s;
        //cout<<"ctype: "<<min_ct<<endl;
        //cout<<"ctime: "<<min_t<<endl;
        return true;
    }

    return false;
}

double interference_volume(
    const Eigen::Vector2d x0[3],
    const Eigen::Vector2d x1[3],
    const double min_sep, const double tolerance,
    double toi, double alpha
    )
{

    if (toi < 0) {
        bool col_flag = time_of_collision(
                    x0[0], x1[0], x0[1], x1[1], x0[2], x1[2], min_sep, tolerance, toi, alpha);
        if (!col_flag){
            return 0.0;
        }
    }

    // compute the relative velocity between
    Eigen::Vector2d vel[3] = { x0[0] - x1[0], x0[1] - x1[1], x0[2] - x1[2] };

    Eigen::Vector2d edge_t0(x0[2] - x0[1]);
    Eigen::Vector2d edge_toi((x0[2] + vel[2] * toi) - (x0[1] + vel[1] * toi));



}

}
