#pragma once

#include <string>

struct Scene;
struct Camera;
struct Integrator;

bool parse_pbrt(const std::string &path, Scene *scene, Camera **camera, Integrator **integrator);
