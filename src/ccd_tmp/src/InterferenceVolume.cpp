// InterferenceVolume.cpp
//
//#define DEBUGCOL

#include "InterferenceVolume.h"
#include "Curve.h"
#include "rpoly.h"
#include <Eigen/Geometry>
#include <iostream>

namespace IAGM {

using namespace std;
using namespace __gnu_cxx;
using namespace Eigen;

InterferenceVolume::InterferenceVolume(Curve* C, Intersections& is, hash_map<size_t, size_t>& vertIndices)
{
    computeInterferenceVolume(C, is, vertIndices);
    computeInterferenceVolumeGradientDev(C, is, vertIndices);
    computeInterferenceVolumeGradientFDV2(C, is, vertIndices);
}

double InterferenceVolume::computeInterferenceVolume(
        const double *x0, const double *x1, Curve* C,
        Intersections& is, hash_map<size_t, size_t>& vertIndices, const bool recompute_time){

    double volume = 0.0;
    for (hash_map<size_t, size_t>::iterator iItr = vertIndices.begin(); iItr != vertIndices.end(); ++iItr) {
        Intersection& i = is[iItr->second];

        Vector2d q[3] = { Vector2d(&x0[2 * i.getVertex()]),
            Vector2d(&x0[2 * i.getEdgeVertex1()]),
            Vector2d(&x0[2 * i.getEdgeVertex2()]) };

        Vector2d p[3] = { Vector2d(&x1[2 * i.getVertex()]),
            Vector2d(&x1[2 * i.getEdgeVertex1()]),
            Vector2d(&x1[2 * i.getEdgeVertex2()]) };

        Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };

        // Which vertex does this summand belong to
        //
        size_t idx = 0;
        if (iItr->first == i.getEdgeVertex1())
            idx = 1;
        else if (iItr->first == i.getEdgeVertex2())
            idx = 2;

        // Intersection normal is the edge perpendicular
        //

        double t = i.getTime();
        double s = i.getAlpha();
        if (recompute_time){
            bool col_flag = getIntersectionTime(q[0], p[0], q[1], p[1], q[2], p[2], s, t, C);
            if (!col_flag){
                volume+=0.0;
                continue;
            }
        }

        Vector2d e((q[2] + v[2] * t) - (q[1] + v[1] * t));
        Vector2d n(e[1], -e[0]);

        // TODO volume for vertex-vertex collision, x0-x1, x0-x2, v[idx]*n_vv
        // where n_vv is in the direction of x0-x1 or x0-x2.
        // The corresponding Vgrad should be adapted too.
        t = max(0.0, t - 1e-6);
        assert((1.0 - t) != 0);

        double tmp = 0;

        if (idx == 0)
            tmp = (1.0 - t) * v[idx].dot(n);
        else
            tmp = (1.0 - t) * ((1 - s) * v[1] + s * v[2]).dot(n);

        // Always orient so that the summand is negative
        //
        if (tmp > 0.0)
            volume -= tmp;
        else
            volume += tmp;
    }
    return volume;
}

void InterferenceVolume::computeInterferenceVolume(Curve* C, Intersections& is, hash_map<size_t, size_t>& vertIndices)
{

    double* x0;
    double* x1;
    C->getStartConfiguration(x0);
    C->getEndConfiguration(x1);

    m_volume = computeInterferenceVolume(x0, x1, C, is, vertIndices);
}

bool InterferenceVolume::getIntersectionTime(Vector2d x00, Vector2d x01,
    Vector2d x10, Vector2d x11, Vector2d x20, Vector2d x21, double& s, double& t, Curve* C)
{

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

    swap(x3x2[0], x3x2[1]);
    x3x2[1] *= -1.0;
    swap(v3v2[0], v3v2[1]);
    v3v2[1] *= -1.0;

    double a0 = x1x2.dot(x3x2);
    double a1 = x1x2.dot(v3v2) + x3x2.dot(v1v2);
    double a2 = v1v2.dot(v3v2);

    vector<double> coeffs;
    if (C->getMinimumSeparation() == 0.0) {
        coeffs.resize(3);
        coeffs[0] = a2;
        coeffs[1] = a1;
        coeffs[2] = a0;
    } else {
        double h = C->getMinimumSeparation();
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
    double minSep = C->getMinimumSeparation();

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
            t = rl_total[i];
            if (t < 0)
                t = 0;
            if (t > 1)
                t = 1;
            //std::cout<<"total roots: "<<cnt_total<<", i: "<<i<<std::endl;
            //std::cout<<"p: ["<<x01[0]<<","<<x01[1]<<"]"<<std::endl;
            //std::cout<<"a: ["<<x11[0]<<","<<x11[1]<<"]"<<std::endl;
            //std::cout<<"b: ["<<x21[0]<<","<<x21[1]<<"]"<<std::endl;

            // Check rl[i] for intersection
            //
            Vector2d p = x1 + v1 * t;
            Vector2d a = x2 + v2 * t;
            Vector2d b = x3 + v3 * t;

            Vector2d ab = b - a;
            Vector2d ap = p - a;

            s = ab.dot(ap) / ab.dot(ab);

            //std::cout<<"s: "<<s<<std::endl;
            // Clamp to edge endpoints
            //
            if (s < 0.00)
                s = 0.0;
            if (s > 1.00)
                s = 1.0;

            // NOTE: If the code misses intersections, it's probably because of
            // this tolerance. There is a Eurographics 2012 paper which attempts to
            // address this, but for now it can be tricky to choose a reliable epsilon.
            //
            double d = (p - (a + ab * s)).squaredNorm() - minSep * minSep;
            //std::cout<<"d: "<<d<<std::endl;
            //std::cout<<"t: "<<t<<std::endl;

            // safeguard test if t is close to 1, check the ending position
            if ((d < C->m_tol) || (de <= 0))
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
                if (t < min_t) {
                    //std::cout<<"smaller time accepted, i: "<<i<<std::endl;
                    min_t = t;
                    min_s = s;
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
        t = max(0.0, min_t - 0);
        s = min_s;
        //cout<<"ctype: "<<min_ct<<endl;
        //cout<<"ctime: "<<min_t<<endl;
        return true;
    }

    return false;
}

// finite difference for volume gradient
void InterferenceVolume::computeInterferenceVolumeGradientFDV2(Curve* C, Intersections& is, hash_map<size_t, size_t>& vertIndices){
    double* x0;
    double* x1;
    C->getStartConfiguration(x0);
    C->getEndConfiguration(x1);

    double volume = computeInterferenceVolume(x0, x1, C, is, vertIndices);

    m_dVdq_fd = VectorXd::Zero(2 * C->getNbrVertices());
    Eigen::VectorXd x1_orig = VectorXd::Zero(2 * C->getNbrVertices());
    for (uint i=0; i < x1_orig.rows(); ++i){
        x1_orig[i] = x1[i];
    }
    for (hash_map<size_t, size_t>::iterator iItr = vertIndices.begin(); iItr != vertIndices.end(); ++iItr) {
        Intersection& i = is[iItr->second];
        VectorXi nodes(3);
        nodes << i.getVertex(), i.getEdgeVertex1(), i.getEdgeVertex2();

        const double eps = 1e-9;
        for (long j = 0; j < 3; ++j) {
            Eigen::Vector2d dv;
            for (long k = 0; k < 2; ++k) {
                x1[2 * nodes[j] + k] = x1_orig[2 * nodes[j] + k] + eps;
                double volume_eps = computeInterferenceVolume(x0, x1, C, is, vertIndices, true);
                dv[k] = (volume_eps - volume) / eps;

                x1[2 * nodes[j] + k] = x1_orig[2 * nodes[j] + k];
            }
             m_dVdq_fd.segment<2>(2 * nodes[j]) = dv;
        }
    }
}

void InterferenceVolume::computeInterferenceVolumeGradientFD(Curve* C, Intersections& is, hash_map<size_t, size_t>& vertIndices)
{
    double* x0;
    double* x1;
    C->getStartConfiguration(x0);
    C->getEndConfiguration(x1);

    m_dVdq_fd = VectorXd::Zero(2 * C->getNbrVertices());
    for (hash_map<size_t, size_t>::iterator iItr = vertIndices.begin(); iItr != vertIndices.end(); ++iItr) {
        Intersection& i = is[iItr->second];

        Vector2d q[3] = { Vector2d(&x0[2 * i.getVertex()]),
            Vector2d(&x0[2 * i.getEdgeVertex1()]),
            Vector2d(&x0[2 * i.getEdgeVertex2()]) };

        Vector2d p[3] = { Vector2d(&x1[2 * i.getVertex()]),
            Vector2d(&x1[2 * i.getEdgeVertex1()]),
            Vector2d(&x1[2 * i.getEdgeVertex2()]) };

        Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };

        // Which vertex does this summand belong to
        //
        size_t idx = 0;
        if (iItr->first == i.getEdgeVertex1())
            idx = 1;
        else if (iItr->first == i.getEdgeVertex2())
            idx = 2;

        // Intersection normal is the edge perpendicular
        //
        Vector2d e((q[2] + v[2] * i.getTime()) - (q[1] + v[1] * i.getTime()));
        Vector2d n(e[1], -e[0]);

        double sign = getSign(idx, i.getTime(), i.getAlpha(), q, v);
        //cout<<"alpha: "<<i.getAlpha()<<", time: "<<i.getTime()<<endl;

        double col_time = i.getTime();
        //col_time =  max(0.0, col_time - 0.0);
        //assert((1.0-col_time)!=0);
        double V = 0;
        if (idx == 0)
            V = (1.0 - col_time) * v[idx].dot(n);
        else
            V = (1.0 - col_time) * ((1 - i.getAlpha()) * v[1] + i.getAlpha() * v[2]).dot(n);
        //double V = v[idx].dot(n);

        // Always orient so that the summand is negative
        //
        //if (V > 0.0)
        //V = -V;
        V = V * sign;

        // finite diff for Vgrad
        const double eps = 1e-9;
        Vector2d dv[3];
        for (long j = 0; j < 3; ++j) {
            for (long k = 0; k < 2; ++k) {
                Vector2d p_eps[3] = { p[0], p[1], p[2] };
                p_eps[j][k] += eps;
                Vector2d v_eps[3] = { p_eps[0] - q[0], p_eps[1] - q[1], p_eps[2] - q[2] };

                double t = 1;
                double s = 0;

                bool col_flag = getIntersectionTime(q[0], p_eps[0], q[1], p_eps[1], q[2], p_eps[2], s, t, C);
                //cout<<"alpha_eps: "<<s<<", time_eps: "<<t<<endl;

                double V_eps = 0;
                if (col_flag)
                //if(true)
                {
                    Vector2d e_eps((q[2] + v_eps[2] * t) - (q[1] + v_eps[1] * t));
                    Vector2d n_eps(e_eps[1], -e_eps[0]);

                    t = max(0.0, t);
                    //assert((1.0-t)!=0);
                    if (idx == 0)
                        V_eps = (1.0 - t) * v_eps[idx].dot(n_eps);
                    else
                        V_eps = (1.0 - t) * ((1 - s) * v_eps[1] + s * v_eps[2]).dot(n_eps);

                    //if (V_eps > 0.0)
                    //V_eps = -V_eps;
                    double sign_eps = getSign(idx, t, s, q, v_eps);
                    V_eps = V_eps * sign_eps;
                } else {
                    V_eps = 0;
                }

                dv[j][k] = (V_eps - V) / eps;
            }
        }
        m_dVdq_fd.segment<2>(2 * i.getVertex()) += dv[0];
        m_dVdq_fd.segment<2>(2 * i.getEdgeVertex1()) += dv[1];
        m_dVdq_fd.segment<2>(2 * i.getEdgeVertex2()) += dv[2];
        // end of finite diff for Vgrad

    } // end of iter hash map
}

void InterferenceVolume::computeInterferenceVolumeGradient(Curve* C, Intersections& is, hash_map<size_t, size_t>& vertIndices)
{
    double* x0;
    double* x1;
    C->getStartConfiguration(x0);
    C->getEndConfiguration(x1);

    m_dVdq = VectorXd::Zero(2 * C->getNbrVertices());
    for (hash_map<size_t, size_t>::iterator iItr = vertIndices.begin(); iItr != vertIndices.end(); ++iItr) {
        Intersection& i = is[iItr->second];

        Vector2d q[3] = { Vector2d(&x0[2 * i.getVertex()]),
            Vector2d(&x0[2 * i.getEdgeVertex1()]),
            Vector2d(&x0[2 * i.getEdgeVertex2()]) };

        Vector2d p[3] = { Vector2d(&x1[2 * i.getVertex()]),
            Vector2d(&x1[2 * i.getEdgeVertex1()]),
            Vector2d(&x1[2 * i.getEdgeVertex2()]) };

        size_t idx = 0;
        if (iItr->first == i.getEdgeVertex1())
            idx = 1;
        else if (iItr->first == i.getEdgeVertex2())
            idx = 2;

        // s is the collision time alpha, not end configration alpha
        // double s = i.getAlpha();

        // get start config normal of edge and end config normal of edge
        Vector2d eStart(q[2] - q[1]);
        Vector2d eEnd(p[2] - p[1]);
        Vector2d nStart(eStart[1], -eStart[0]);
        Vector2d nEnd(eEnd[1], -eEnd[0]);

        Vector2d ap_start = q[0] - q[1];
        Vector2d ab_no = p[2] - p[1];
        Vector2d ap_no = p[0] - p[1];
        double s_no = ab_no.dot(ap_no) / ab_no.dot(ab_no);
        double s_start = eStart.dot(ap_start) / eStart.dot(eStart);
        //double d_no = (p[0]-(p[1] + ab_no * s_no)).norm() - m_minSep;
        //double d_no;
        //double startSign = ap_start.dot(nStart);
        double startSign = (q[0] - (q[1] + eStart * s_start)).dot(nStart);
        double endSign = (p[0] - (p[1] + ab_no * s_no)).dot(nEnd);

        Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };
        Vector2d e_test((q[2] + v[2] * i.getTime()) - (q[1] + v[1] * i.getTime()));
        Vector2d n_test(e_test[1], -e_test[0]);

        //s_no = s;
        //nEnd = n_test;
        if (s_no < 0)
            s_no = 0.0;
        if (s_no > 1)
            s_no = 1.0;
        Vector2d dvdq0, dvdq1, dvdq2;
        if (startSign == 0) {
            if (endSign > 0)
                nEnd = -nEnd;
            if (idx == 0) {
                dvdq0 = nEnd;
            } else {
                dvdq1 = -(1 - s_no) * nEnd;
                dvdq2 = -s_no * nEnd;
            }
            //dvdq0 = -dddq0*n.norm();
            //dvdq1 = -dddq1*n.norm();// - dist.norm()*dndq1 - m_minSep*dndq1;
            //dvdq2 = -dddq2*n.norm();// - dist.norm()*dndq2 - m_minSep*dndq2;
        }
        //	d_no = -m_minSep - (p[0]-(p[1] + ab_no * s_no)).norm();
        else if (endSign == 0) {

            if (startSign < 0)
                nEnd = -nEnd;
            if (idx == 0) {
                dvdq0 = nEnd;
            } else {
                dvdq1 = -(1 - s_no) * nEnd;
                dvdq2 = -s_no * nEnd;
            }
            //dvdq0 = dddq0*0;
            //dvdq1 = - m_minSep*dndq1*0;
            //dvdq2 = - m_minSep*dndq2*0;
        }
        //	d_no = -m_minSep;// + (p[0]-(p[1] + ab_no * s_no)).norm();

        if (startSign * endSign > 0) {

            // double d_tmp = -m_minSep + (p[0]-(p[1] + ab_no * s_no)).norm();
            if (endSign < 0)
                nEnd = -nEnd;
            //if(d_tmp < 0 && s_no < 1 && s_no > 0){
            //if(d_tmp < 0){
            if (idx == 0) {
                dvdq0 = nEnd;
            } else {
                dvdq1 = -(1 - s_no) * nEnd;
                dvdq2 = -s_no * nEnd;
            }
            //}
            //else
            //{

            //	dvdq0 = dddq0*0;
            //	dvdq1 = - m_minSep*dndq1*0;
            //	dvdq2 = - m_minSep*dndq2*0;

            //}
            //dvdq0 = dddq0*n.norm();
            //dvdq1 = dddq1*n.norm();// + dist.norm()*dndq1 - m_minSep*dndq1;
            //dvdq2 = dddq2*n.norm();// + dist.norm()*dndq2 - m_minSep*dndq2;

        }
        //	d_no = -m_minSep + (p[0]-(p[1] + ab_no * s_no)).norm();
        else if (startSign * endSign < 0) {

            if (endSign > 0)
                nEnd = -nEnd;
            if (idx == 0) {
                dvdq0 = nEnd;
            } else {
                dvdq1 = -(1 - s_no) * nEnd;
                dvdq2 = -s_no * nEnd;
            }

            //dvdq0 = -dddq0*n.norm();
            //dvdq1 = -dddq1*n.norm();// - dist.norm()*dndq1 - m_minSep*dndq1;
            //dvdq2 = -dddq2*n.norm();// - dist.norm()*dndq2 - m_minSep*dndq2;
        }
        //	d_no = -m_minSep - (p[0]-(p[1] + ab_no * s_no)).norm();
        if (idx == 0) {
            m_dVdq.segment<2>(2 * i.getVertex()) += dvdq0;
        } else {
            m_dVdq.segment<2>(2 * i.getEdgeVertex1()) += dvdq1;
            m_dVdq.segment<2>(2 * i.getEdgeVertex2()) += dvdq2;
        }
    }
}

void InterferenceVolume::computeInterferenceVolumeGradientDev(Curve* C, Intersections& is, hash_map<size_t, size_t>& vertIndices)
{
    double* x0;
    double* x1;
    C->getStartConfiguration(x0);
    C->getEndConfiguration(x1);

    // These gradients are with respect to low-level DOFs
    //
    m_dVdq = VectorXd::Zero(2 * C->getNbrVertices());
    for (hash_map<size_t, size_t>::iterator iItr = vertIndices.begin(); iItr != vertIndices.end(); ++iItr) {
        Intersection& i = is[iItr->second];

        Vector2d q[3] = { Vector2d(&x0[2 * i.getVertex()]),
            Vector2d(&x0[2 * i.getEdgeVertex1()]),
            Vector2d(&x0[2 * i.getEdgeVertex2()]) };

        Vector2d p[3] = { Vector2d(&x1[2 * i.getVertex()]),
            Vector2d(&x1[2 * i.getEdgeVertex1()]),
            Vector2d(&x1[2 * i.getEdgeVertex2()]) };

        Vector2d v[3] = { p[0] - q[0], p[1] - q[1], p[2] - q[2] };

        size_t idx = 0;
        if (iItr->first == i.getEdgeVertex1())
            idx = 1;
        else if (iItr->first == i.getEdgeVertex2())
            idx = 2;

        // Intersection normal is the edge perpendicular
        //
        Vector2d e((q[2] + v[2] * i.getTime()) - (q[1] + v[1] * i.getTime()));
        Vector2d n(e[1], -e[0]);

        // get dtdq
        Vector2d dtdq[3];
        getTimeGradient(i, 0, x0, x1, dtdq[0], C);
        getTimeGradient(i, 1, x0, x1, dtdq[1], C);
        getTimeGradient(i, 2, x0, x1, dtdq[2], C);

#ifdef DEBUGCOL
        // FD for dtdq
        Vector2d dtdq_eps[3];
        const double eps = 1e-9;
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 2; ++k) {
                Vector2d p_eps[3] = { p[0], p[1], p[2] };
                p_eps[j][k] += eps;
                Vector2d v_eps[3] = { p_eps[0] - q[0], p_eps[1] - q[1], p_eps[2] - q[2] };

                double t_eps = 1;
                double s_eps = 0;

                bool col_flag = getIntersectionTime(q[0], p_eps[0], q[1], p_eps[1], q[2], p_eps[2], s_eps, t_eps, C);
                dtdq_eps[j][k] = (t_eps - i.getTime()) / eps;
            }
        }

        cout << "Col Time:  " << i.getTime() << endl;
        cout << "Col Alpha: " << i.getAlpha() << endl;
        cout << "dtdq1    : " << dtdq[0].transpose() << endl;
        cout << "dtdq1_eps: " << dtdq_eps[0].transpose() << endl;
        cout << "dtdq2    : " << dtdq[1].transpose() << endl;
        cout << "dtdq2_eps: " << dtdq_eps[1].transpose() << endl;
        cout << "dtdq3    : " << dtdq[2].transpose() << endl;
        cout << "dtdq3_eps: " << dtdq_eps[2].transpose() << endl;
#endif

        // Constants
        double delta1, delta2, delta3;
        delta1 = (idx == 0) ? 1.0 : 0.0;
        delta2 = (idx == 1) ? 1.0 : 0.0;
        delta3 = (idx == 2) ? 1.0 : 0.0;

        MatrixXd Idm(2, 2);
        Idm(0, 0) = 1;
        Idm(1, 0) = 0;
        Idm(0, 1) = 0;
        Idm(1, 1) = 1;

        double sign = getSign(idx, i.getTime(), i.getAlpha(), q, v);

        MatrixXd TM(2, 2);
        TM(0, 0) = 0;
        TM(1, 1) = 0;
        TM(0, 1) = 1;
        TM(1, 0) = -1;

        // get dndq
        MatrixXd dndq1(2, 2);
        dndq1 = (v[2] - v[1]) * dtdq[0].transpose();
        dndq1 = TM * dndq1;

        MatrixXd dndq2(2, 2);
        dndq2 = (v[2] - v[1]) * dtdq[1].transpose() - Idm * i.getTime();
        dndq2 = TM * dndq2;

        MatrixXd dndq3(2, 2);
        dndq3 = (v[2] - v[1]) * dtdq[2].transpose() + Idm * i.getTime();
        dndq3 = TM * dndq3;

#ifdef DEBUGCOL
        // FD for dndq
        MatrixXd dndq1_eps(2, 2);
        MatrixXd dndq2_eps(2, 2);
        MatrixXd dndq3_eps(2, 2);

        // FD for dndq_1
        for (size_t k = 0; k < 2; ++k) {
            Vector2d p_eps[3] = { p[0], p[1], p[2] };
            p_eps[0][k] += eps;
            Vector2d v_eps[3] = { p_eps[0] - q[0], p_eps[1] - q[1], p_eps[2] - q[2] };

            double t_eps = 1;
            double s_eps = 0;

            bool col_flag = getIntersectionTime(q[0], p_eps[0], q[1], p_eps[1], q[2], p_eps[2], s_eps, t_eps, C);

            Vector2d e_eps((q[2] + v_eps[2] * t_eps) - (q[1] + v_eps[1] * t_eps));
            Vector2d n_eps(e_eps[1], -e_eps[0]);

            Vector2d tmpdndq = (n_eps - n) / eps;
            dndq1_eps(0, k) = tmpdndq(0);
            dndq1_eps(1, k) = tmpdndq(1);
        }
        // FD for dndq_2
        for (size_t k = 0; k < 2; ++k) {
            Vector2d p_eps[3] = { p[0], p[1], p[2] };
            p_eps[1][k] += eps;
            Vector2d v_eps[3] = { p_eps[0] - q[0], p_eps[1] - q[1], p_eps[2] - q[2] };

            double t_eps = 1;
            double s_eps = 0;

            bool col_flag = getIntersectionTime(q[0], p_eps[0], q[1], p_eps[1], q[2], p_eps[2], s_eps, t_eps, C);

            Vector2d e_eps((q[2] + v_eps[2] * t_eps) - (q[1] + v_eps[1] * t_eps));
            Vector2d n_eps(e_eps[1], -e_eps[0]);

            Vector2d tmpdndq = (n_eps - n) / eps;
            dndq2_eps(0, k) = tmpdndq(0);
            dndq2_eps(1, k) = tmpdndq(1);
        }
        // FD for dndq_3
        for (size_t k = 0; k < 2; ++k) {
            Vector2d p_eps[3] = { p[0], p[1], p[2] };
            p_eps[2][k] += eps;
            Vector2d v_eps[3] = { p_eps[0] - q[0], p_eps[1] - q[1], p_eps[2] - q[2] };

            double t_eps = 1;
            double s_eps = 0;

            bool col_flag = getIntersectionTime(q[0], p_eps[0], q[1], p_eps[1], q[2], p_eps[2], s_eps, t_eps, C);

            Vector2d e_eps((q[2] + v_eps[2] * t_eps) - (q[1] + v_eps[1] * t_eps));
            Vector2d n_eps(e_eps[1], -e_eps[0]);

            Vector2d tmpdndq = (n_eps - n) / eps;
            dndq3_eps(0, k) = tmpdndq(0);
            dndq3_eps(1, k) = tmpdndq(1);
        }

        cout << "dndq1:     \n"
             << dndq1 << endl;
        cout << "dndq1_eps: \n"
             << dndq1_eps << endl;
        cout << "dndq2:     \n"
             << dndq2 << endl;
        cout << "dndq2_eps: \n"
             << dndq2_eps << endl;
        cout << "dndq3:     \n"
             << dndq3 << endl;
        cout << "dndq3_eps: \n"
             << dndq3_eps << endl;
#endif

        // get dsdq
        Vector2d dsdq[3];
        if (i.getAlpha() == 0 || i.getAlpha() == 1) {
            dsdq[0] = dsdq[0] * 0;
            dsdq[1] = dsdq[1] * 0;
            dsdq[2] = dsdq[2] * 0;
        } else {
            getAlphaGradient(0, q[0] + i.getTime() * v[0], q[1] + i.getTime() * v[1], q[2] + i.getTime() * v[2], dsdq[0]);
            getAlphaGradient(1, q[0] + i.getTime() * v[0], q[1] + i.getTime() * v[1], q[2] + i.getTime() * v[2], dsdq[1]);
            getAlphaGradient(2, q[0] + i.getTime() * v[0], q[1] + i.getTime() * v[1], q[2] + i.getTime() * v[2], dsdq[2]);
        }

#ifdef DEBUGCOL
        // FD for dsdq
        Vector2d dsdq_eps[3];
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 2; ++k) {
                Vector2d r_eps[3] = { q[0] + i.getTime() * v[0], q[1] + i.getTime() * v[1], q[2] + i.getTime() * v[2] };
                r_eps[j][k] += eps;

                Vector2d ap_eps = r_eps[0] - r_eps[1];
                Vector2d ab_eps = r_eps[2] - r_eps[1];

                double s_eps = ap_eps.dot(ab_eps) / ab_eps.dot(ab_eps);
                if (s_eps > 1.0)
                    s_eps = 1;
                if (s_eps < 0.0)
                    s_eps = 0;

                dsdq_eps[j][k] = (s_eps - i.getAlpha()) / eps;
            }
        }
        cout << "dsdq1:     " << dsdq[0].transpose() << endl;
        cout << "dsdq1_eps: " << dsdq_eps[0].transpose() << endl;
        cout << "dsdq2:     " << dsdq[1].transpose() << endl;
        cout << "dsdq2_eps: " << dsdq_eps[1].transpose() << endl;
        cout << "dsdq3:     " << dsdq[2].transpose() << endl;
        cout << "dsdq3_eps: " << dsdq_eps[2].transpose() << endl;
#endif

        // get dwdq
        Vector2d dwdq[3];
        dwdq[0] = dwdq[0] * 0;
        dwdq[1] = dwdq[1] * 0;
        dwdq[2] = dwdq[2] * 0;
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 3; ++k) {
                if (j == k)
                    dwdq[j] += (v[k] * dtdq[j].transpose() + i.getTime() * Matrix2d::Identity()).transpose() * dsdq[k];
                else
                    dwdq[j] += (v[k] * dtdq[j].transpose()).transpose() * dsdq[k];
            }
        }

#ifdef DEBUGCOL
        // FD for dwdq
        Vector2d dwdq_eps[3];
        for (size_t j = 0; j < 3; ++j) {
            for (size_t k = 0; k < 2; ++k) {
                Vector2d p_eps[3] = { p[0], p[1], p[2] };
                p_eps[j][k] += eps;
                Vector2d v_eps[3] = { p_eps[0] - q[0], p_eps[1] - q[1], p_eps[2] - q[2] };

                double t_eps = 1;
                double s_eps = 0;

                bool col_flag = getIntersectionTime(q[0], p_eps[0], q[1], p_eps[1], q[2], p_eps[2], s_eps, t_eps, C);
                dwdq_eps[j][k] = (s_eps - i.getAlpha()) / eps;
            }
        }
        cout << "dwdq1:     " << dwdq[0].transpose() << endl;
        cout << "dwdq1_eps: " << dwdq_eps[0].transpose() << endl;
        cout << "dwdq2:     " << dwdq[1].transpose() << endl;
        cout << "dwdq2_eps: " << dwdq_eps[1].transpose() << endl;
        cout << "dwdq3:     " << dwdq[2].transpose() << endl;
        cout << "dwdq3_eps: " << dwdq_eps[2].transpose() << endl;
#endif

        // compute volume gradient
        /*
                m_dVdq.segment<2>(2*i.getVertex())      += sign * (-dtdq[0] * v[idx].dot(n) + 
                      delta1 * (1 - i.getTime()) * Idm * n + 
                      (1 - i.getTime()) * dndq1.transpose() * v[idx]);

                m_dVdq.segment<2>(2*i.getEdgeVertex1()) += sign * (-dtdq[1] * v[idx].dot(n) + 
                      delta2 * (1 - i.getTime()) * Idm * n +  
                      (1 - i.getTime()) * dndq2.transpose() * v[idx]);

                m_dVdq.segment<2>(2*i.getEdgeVertex2()) += sign * (-dtdq[2] * v[idx].dot(n) + 
                      delta3 * (1 - i.getTime()) * Idm * n +  
                      (1 - i.getTime()) * dndq3.transpose() * v[idx]);
                */

        if (idx == 0) {
            m_dVdq.segment<2>(2 * i.getVertex()) += sign * (-dtdq[0] * v[idx].dot(n) + delta1 * (1 - i.getTime()) * Idm * n + (1 - i.getTime()) * dndq1.transpose() * v[idx]);

            m_dVdq.segment<2>(2 * i.getEdgeVertex1()) += sign * (-dtdq[1] * v[idx].dot(n) + delta2 * (1 - i.getTime()) * Idm * n + (1 - i.getTime()) * dndq2.transpose() * v[idx]);

            m_dVdq.segment<2>(2 * i.getEdgeVertex2()) += sign * (-dtdq[2] * v[idx].dot(n) + delta3 * (1 - i.getTime()) * Idm * n + (1 - i.getTime()) * dndq3.transpose() * v[idx]);
            //m_dVdq.segment<2>(2*i.getVertex())      += sign * n;
        } else {
            m_dVdq.segment<2>(2 * i.getVertex()) += sign * (
                                                               // v[1] part
                                                               -dtdq[0] * (1 - i.getAlpha()) * v[1].dot(n) + 0 * n + (1 - i.getTime()) * (1 - i.getAlpha()) * dndq1.transpose() * v[1] + -dwdq[0] * (1 - i.getTime()) * v[1].dot(n) +
                                                               // v[2] part
                                                               -dtdq[0] * i.getAlpha() * v[2].dot(n) + 0 * n + (1 - i.getTime()) * i.getAlpha() * dndq1.transpose() * v[2] + dwdq[0] * (1 - i.getTime()) * v[2].dot(n));

            m_dVdq.segment<2>(2 * i.getEdgeVertex1()) += sign * (
                                                                    // v[1] part
                                                                    -dtdq[1] * (1 - i.getAlpha()) * v[1].dot(n) + (1 - i.getTime()) * (1 - i.getAlpha()) * n + (1 - i.getTime()) * (1 - i.getAlpha()) * dndq2.transpose() * v[1] + -dwdq[1] * (1 - i.getTime()) * v[1].dot(n) +
                                                                    // v[2] part
                                                                    -dtdq[1] * i.getAlpha() * v[2].dot(n) + 0 * n + (1 - i.getTime()) * i.getAlpha() * dndq2.transpose() * v[2] + dwdq[1] * (1 - i.getTime()) * v[2].dot(n));

            m_dVdq.segment<2>(2 * i.getEdgeVertex2()) += sign * (
                                                                    // v[1] part
                                                                    -dtdq[2] * (1 - i.getAlpha()) * v[1].dot(n) + 0 * n + (1 - i.getTime()) * (1 - i.getAlpha()) * dndq3.transpose() * v[1] + -dwdq[2] * (1 - i.getTime()) * v[1].dot(n) +
                                                                    // v[2] part
                                                                    -dtdq[2] * i.getAlpha() * v[2].dot(n) + (1 - i.getTime()) * i.getAlpha() * n + (1 - i.getTime()) * i.getAlpha() * dndq3.transpose() * v[2] + dwdq[2] * (1 - i.getTime()) * v[2].dot(n));
            //m_dVdq.segment<2>(2*i.getEdgeVertex1()) += sign * (1 - i.getAlpha()) * n;
            //m_dVdq.segment<2>(2*i.getEdgeVertex2()) += sign * i.getAlpha() * n;
        }
    }
}

double InterferenceVolume::getSign(size_t idx, double t, double a, Vector2d* q, Vector2d* v)
{
    Vector2d z1 = q[1] + v[1] * t;
    Vector2d z2 = q[2] + v[2] * t;
    Vector2d e(z2 - z1);
    Vector2d n(e[1], -e[0]);

    double sign = (idx) ? -1.0 : 1.0;
    double rv = (v[0] - (1.0 - a) * v[1] - a * v[2]).dot(n);
    if (rv > 0.0)
        return -sign;

    // Should never have an intersection between two elements
    // with zero relative velocity
    //
    // TODO: the relative speed maybe zero when vertex-vertex
    // collision is detected, for v-v collision, should use
    // direction of x0-x1 or x0-x2 instead of n to do
    // relative speed test
    if (rv == 0.0)
        cout << "relative speed: " << rv << endl;
    assert(rv != 0.0);

    return sign;
}

void InterferenceVolume::computeHandleGradient(Curve* C, std::set<size_t>& movable)
{
    // Project our gradient into the modeling subspace
    //
    m_dVdp = VectorXd::Zero(2 * C->getNbrCtrlVertices());
    C->getSubspaceVector(m_dVdq.data(), m_dVdp.data());

    m_dVdp_fd = VectorXd::Zero(2 * C->getNbrCtrlVertices());
    C->getSubspaceVector(m_dVdq_fd.data(), m_dVdp_fd.data());

    set<size_t> all;
    for (size_t i = 0; i < C->getNbrCtrlVertices(); ++i)
        all.insert(i);

    // Get the immovable vertices
    //
    set<size_t> immovable;
    set_difference(all.begin(), all.end(), movable.begin(), movable.end(),
        inserter(immovable, immovable.end()));

    // Zero out gradient for unselected vertices (they don't move)
    //
    for (set<size_t>::iterator sItr = immovable.begin(); sItr != immovable.end(); ++sItr)
        for (size_t i = 0; i < 2; ++i){
            m_dVdp[2 * (*sItr) + i] = 0.0;
            m_dVdp_fd[2 * (*sItr) + i] = 0.0;
        }
}

double InterferenceVolume::getBarycentricArea(size_t idx, size_t nbrVertices, double* x)
{
    // TODO: This assumes one spline, essentially. Ideally it should check if the vertex
    // is at the beginning or end of a new spline, not just at 0 or nbrVertices
    //
    double A = 0.0;

    if (idx > 0)
        A += 0.5 * Vector2d(x[2 * (idx - 1)] - x[2 * idx], x[2 * (idx - 1) + 1] - x[2 * idx + 1]).norm();

    if (idx < (nbrVertices - 1))
        A += 0.5 * Vector2d(x[2 * idx] - x[2 * (idx + 1)], x[2 * idx + 1] - x[2 * (idx + 1) + 1]).norm();

    return A;
}

double InterferenceVolume::getVolume(Vector2d& r,
    Vector2d x0, Vector2d x1, Vector2d x2,
    Vector2d v0, Vector2d v1, Vector2d v2,
    double A, double t)
{
    Vector2d z1 = x1 + v1 * t;
    Vector2d z2 = x2 + v2 * t;
    Vector2d e(z2 - z1);
    Vector2d n(e[1], -e[0]);

    return (1.0 - t) * r.dot(n) * A;
}

void InterferenceVolume::getDnDt(Vector2d* v, Vector2d& dndt)
{
    Vector2d e = v[2] - v[1];
    dndt = Vector2d(e[1], -e[0]);
}

void InterferenceVolume::getAlphaGradient(Intersection& i, size_t idx, double* x0, double* x1, Vector2d& dsdq)
{
    Vector2d a(x1[2 * i.getEdgeVertex1()], x1[2 * i.getEdgeVertex1() + 1]);
    Vector2d b(x1[2 * i.getEdgeVertex2()], x1[2 * i.getEdgeVertex2() + 1]);
    Vector2d p(x1[2 * i.getVertex()], x1[2 * i.getVertex() + 1]);

    // squared norm of edge ab
    double lenE = (b - a).dot(b - a);
    // squared norm of ap
    double lenV = (b - a).dot(p - a);

    // TODO seems to be wrong, should be the following function
    if (idx == 0) {
        dsdq = (b - a) / lenE;
    } else if (idx == 1) {
        dsdq = (p - a) / lenE - 2 * (b - a) * lenV / (lenE * lenE);
    } else {
        dsdq = (2 * a - p - b) / lenE + 2 * (b - a) * lenV / (lenE * lenE);
    }
}

void InterferenceVolume::getAlphaGradient(size_t idx, Vector2d p, Vector2d a, Vector2d b, Vector2d& dsdq)
{
    // squared norm of edge ab
    double lenE = (b - a).dot(b - a);
    // squared norm of ap
    double lenV = (b - a).dot(p - a);

    if (idx == 0) {
        dsdq = (b - a) / lenE;
    } else if (idx == 1) {
        dsdq = (2 * a - p - b) / lenE + 2 * (b - a) * lenV / (lenE * lenE);
    } else {
        dsdq = (p - a) / lenE - 2 * (b - a) * lenV / (lenE * lenE);
    }
}

void InterferenceVolume::getTimeGradient(Intersection& i, size_t idx, double* x0, double* x1, Vector2d& dtdq, Curve* C)
{
    Vector2d q1(x0[2 * i.getVertex()], x0[2 * i.getVertex() + 1]);
    Vector2d q2(x0[2 * i.getEdgeVertex1()], x0[2 * i.getEdgeVertex1() + 1]);
    Vector2d q3(x0[2 * i.getEdgeVertex2()], x0[2 * i.getEdgeVertex2() + 1]);
    Vector2d p1(x1[2 * i.getVertex()], x1[2 * i.getVertex() + 1]);
    Vector2d p2(x1[2 * i.getEdgeVertex1()], x1[2 * i.getEdgeVertex1() + 1]);
    Vector2d p3(x1[2 * i.getEdgeVertex2()], x1[2 * i.getEdgeVertex2() + 1]);

    Vector2d v1(p1 - q1);
    Vector2d v2(p2 - q2);
    Vector2d v3(p3 - q3);

    Vector2d x1x2 = q1 - q2;
    Vector2d v1v2 = v1 - v2;
    Vector2d x3x2 = q3 - q2;
    Vector2d v3v2 = v3 - v2;

    Vector2d x1x3 = q1 - q3;
    Vector2d v1v3 = v1 - v3;

    swap(x3x2[0], x3x2[1]);
    x3x2[1] *= -1.0;
    swap(v3v2[0], v3v2[1]);
    v3v2[1] *= -1.0;

    double a2 = v1v2.dot(v3v2);
    double a1 = x1x2.dot(v3v2) + v1v2.dot(x3x2);
    double a0 = x1x2.dot(x3x2);

    if (C->getMinimumSeparation() == 0.0) {
        double scale = -1.0 / (2.0 * a2 * i.getTime() + a1);

        Vector2d z1 = q2 + v2 * i.getTime();
        Vector2d z2 = q3 + v3 * i.getTime();
        Vector2d e(z2 - z1);
        Vector2d n(e[1], -e[0]);

        Vector2d num;

        if (idx == 0)
            num = n * i.getTime();
        else if (idx == 1)
            num = -(1.0 - i.getAlpha()) * n * i.getTime();
        else
            num = -i.getAlpha() * n * i.getTime();

        dtdq = scale * num;
    } else if (i.getAlpha() < 1 && i.getAlpha() > 0 && C->getMinimumSeparation() != 0.0) {
        double t = i.getTime();
        double h = C->getMinimumSeparation();
        MatrixXd TM(2, 2);
        TM(0, 0) = 0;
        TM(1, 1) = 0;
        TM(0, 1) = -1;
        TM(1, 0) = 1;

        double scale = -1.0 / (4.0 * (a2 * a2) * t * t * t + 3.0 * (2.0 * a1 * a2) * t * t + 2.0 * (a1 * a1 + 2.0 * a0 * a2 - h * h * v3v2.dot(v3v2)) * t + (2.0 * a0 * a1 - h * h * 2.0 * x3x2.dot(v3v2)));

        Vector2d num;
        if (idx == 0)
            num = (2.0 * a2 * v3v2) * t * t * t * t + (2.0 * a1 * v3v2 + 2.0 * a2 * x3x2) * t * t * t + (2.0 * a1 * x3x2 + 2.0 * a0 * v3v2) * t * t + (2.0 * a0 * x3x2) * t;
        else if (idx == 1)
            num = (2.0 * a2 * (-v3v2 - TM * v1v2)) * t * t * t * t + (2.0 * a1 * (-v3v2 - TM * v1v2) + 2.0 * a2 * (-TM * x1x2 - x3x2)) * t * t * t + (2.0 * a1 * (-TM * x1x2 - x3x2) + 2.0 * a0 * (-v3v2 - TM * v1v2) + h * h * 2.0 * (v3 - v2)) * t * t + (2.0 * a0 * (-TM * x1x2 - x3x2) + h * h * 2.0 * TM * x3x2) * t;
        else
            num = (2.0 * a2 * TM * v1v2) * t * t * t * t + (2.0 * a1 * (TM * v1v2) + 2.0 * a2 * (TM * x1x2)) * t * t * t + (2.0 * a1 * (TM * x1x2) + 2.0 * a0 * (TM * v1v2) - h * h * 2.0 * (v3 - v2)) * t * t + (2.0 * a0 * (TM * x1x2) - h * h * 2.0 * TM * x3x2) * t;

        dtdq = scale * num;
    } else if (i.getAlpha() == 1.0 && C->getMinimumSeparation() != 0.0) {
        double t = i.getTime();
        //   double h = C->getMinimumSeparation();

        a2 = v1v3.dot(v1v3);
        a1 = 2.0 * x1x3.dot(v1v3);

        double scale = -1.0 / (2.0 * a2 * t + a1);

        Vector2d num;
        if (idx == 0)
            num = (2.0 * v1v3) * t * t + (2.0 * x1x3) * t;
        else if (idx == 1)
            num = 0 * num;
        else
            num = (-2.0 * v1v3) * t * t + (-2.0 * x1x3) * t;

        dtdq = scale * num;
    } else if (i.getAlpha() == 0.0 && C->getMinimumSeparation() != 0.0) {
        double t = i.getTime();
        //   double h = C->getMinimumSeparation();

        a2 = v1v2.dot(v1v2);
        a1 = 2.0 * x1x2.dot(v1v2);

        double scale = -1.0 / (2.0 * a2 * t + a1);

        Vector2d num;
        if (idx == 0)
            num = (2.0 * v1v2) * t * t + (2.0 * x1x2) * t;
        else if (idx == 1)
            num = (-2.0 * v1v2) * t * t + (-2.0 * x1x2) * t;
        else
            num = 0 * num;

        dtdq = scale * num;
    }
}

} // end namespace IAGM
