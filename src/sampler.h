#pragma once

#include "camera.h"
#include "point.h"

struct Sampler
{
    int samples_per_pixel; // TODO: int64_t?
    int index;
    
    // TODO: int64_t?
    Sampler(int samples_per_pixel) : samples_per_pixel(samples_per_pixel), index(0) {}

    CameraSample get_camera_sample(const Point2i &pixel);

    virtual void start_pixel(const Point2i &pixel);
    virtual bool start_next_sample();

    virtual float get_1d() = 0;
    virtual Point2f get_2d() = 0;
};

struct RandomSampler : public Sampler
{
    RandomSampler(int samples_per_pixel) : Sampler(samples_per_pixel) {}

    float get_1d() override;
    Point2f get_2d() override;
};
