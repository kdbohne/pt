#include "matrix.h"
#include "common.h"
#include <cstring>
#include <cmath>

Matrix4x4::Matrix4x4() : Matrix4x4(1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 1, 0,
                                   0, 0, 0, 1)
{
}

Matrix4x4::Matrix4x4(float mat[4][4])
{
    std::memcpy(m, mat, sizeof(float) * 16);
}

Matrix4x4::Matrix4x4(float m00, float m10, float m20, float m30,
                     float m01, float m11, float m21, float m31,
                     float m02, float m12, float m22, float m32,
                     float m03, float m13, float m23, float m33)
{
    m[0][0] = m00;
    m[0][1] = m01;
    m[0][2] = m02;
    m[0][3] = m03;

    m[1][0] = m10;
    m[1][1] = m11;
    m[1][2] = m12;
    m[1][3] = m13;

    m[2][0] = m20;
    m[2][1] = m21;
    m[2][2] = m22;
    m[2][3] = m23;

    m[3][0] = m30;
    m[3][1] = m31;
    m[3][2] = m32;
    m[3][3] = m33;
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4 &mat) const
{
    Matrix4x4 r;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            r.m[i][j] = m[0][j] * mat.m[i][0] +
                m[1][j] * mat.m[i][1] +
                m[2][j] * mat.m[i][2] +
                m[3][j] * mat.m[i][3];
        }
    }
    return r;
}

Matrix4x4 transpose(const Matrix4x4 &m)
{
    return Matrix4x4(m.m[0][0], m.m[0][1], m.m[0][2], m.m[0][3],
                     m.m[1][0], m.m[1][1], m.m[1][2], m.m[1][3],
                     m.m[2][0], m.m[2][1], m.m[2][2], m.m[2][3],
                     m.m[3][0], m.m[3][1], m.m[3][2], m.m[3][3]);
}

// NOTE: this is Inverse() from pbrt-v3.
Matrix4x4 inverse(const Matrix4x4 &m)
{
    int indxc[4], indxr[4];
    int ipiv[4] = {0, 0, 0, 0};
    float minv[4][4];
    std::memcpy(minv, m.m, 4 * 4 * sizeof(float));
    for (int i = 0; i < 4; i++) {
        int irow = 0, icol = 0;
        float big = 0.f;
        // Choose pivot
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (std::abs(minv[j][k]) >= big) {
                            big = float(std::abs(minv[j][k]));
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k] > 1)
                        error("Singular matrix in MatrixInvert%s", "");
                }
            }
        }
        ++ipiv[icol];
        // Swap rows _irow_ and _icol_ for pivot
        if (irow != icol) {
            for (int k = 0; k < 4; ++k) std::swap(minv[irow][k], minv[icol][k]);
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (minv[icol][icol] == 0.f) error("Singular matrix in MatrixInvert%s", "");

        // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        float pivinv = 1. / minv[icol][icol];
        minv[icol][icol] = 1.;
        for (int j = 0; j < 4; j++) minv[icol][j] *= pivinv;

        // Subtract this row from others to zero out their columns
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                float save = minv[j][icol];
                minv[j][icol] = 0;
                for (int k = 0; k < 4; k++) minv[j][k] -= minv[icol][k] * save;
            }
        }
    }
    // Swap columns to reflect permutation
    for (int j = 3; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
        }
    }
    return Matrix4x4(minv);
}
