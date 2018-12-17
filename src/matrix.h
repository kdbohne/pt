#pragma once

struct Matrix4x4
{
    float m[4][4];

    Matrix4x4();
    Matrix4x4(float mat[4][4]);
    Matrix4x4(float m00, float m10, float m20, float m30,
              float m01, float m11, float m21, float m31,
              float m02, float m12, float m22, float m32,
              float m03, float m13, float m23, float m33);

    Matrix4x4 operator*(const Matrix4x4 &mat) const;
};

Matrix4x4 transpose(const Matrix4x4 &m);
Matrix4x4 inverse(const Matrix4x4 &m);
