#include "write_png.hpp"

#include <igl_stb_image.h>

namespace swr {

bool write_png(
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& R,
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& G,
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& B,
    const Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>& A,
    const std::string& filename)
{
    assert(R.rows() == G.rows() && R.cols() == G.cols());
    assert(G.rows() == B.rows() && G.cols() == B.cols());
    assert(B.rows() == A.rows() && B.cols() == A.cols());

    const int w = R.rows();               // Image width
    const int h = R.cols();               // Image height
    const int comp = 4;                   // 4 Channels Red, Green, Blue, Alpha
    const int stride_in_bytes = w * comp; // Length of one row in bytes
    std::vector<uint8_t> image(w * h * comp, 0); // The image itself

    for (unsigned wi = 0; wi < w; ++wi) {
        for (unsigned hi = 0; hi < h; ++hi) {
            // Flip the image too
            image[(hi * w * comp) + (wi * comp) + 0] = R(wi, h - 1 - hi);
            image[(hi * w * comp) + (wi * comp) + 1] = G(wi, h - 1 - hi);
            image[(hi * w * comp) + (wi * comp) + 2] = B(wi, h - 1 - hi);
            image[(hi * w * comp) + (wi * comp) + 3] = A(wi, h - 1 - hi);
        }
    }

    igl::stbi_write_png(
        filename.c_str(), w, h, comp, image.data(), stride_in_bytes);

    return true;
}

} // namespace swr
