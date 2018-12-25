#pragma once

#include "vector.h"
#include "point.h"

float uniform_float();

Vector3f uniform_sample_sphere(const Point2f &u);
Point2f uniform_sample_triangle(const Point2f &u);
Point2f concentric_sample_disk(const Point2f &u);
Vector3f cosine_sample_hemisphere(const Point2f &u);

float uniform_cone_pdf(float cos_theta_max);
