#include "sampler.h"
#include "sampling.h"

CameraSample Sampler::get_camera_sample(const Point2i &pixel)
{
    CameraSample cs;
    cs.p_film = (Point2f)pixel + get_2d();
    cs.time = get_1d();

    return cs;
}

void Sampler::start_pixel(const Point2i &pixel)
{
    index = 0;
}

bool Sampler::start_next_sample()
{
    return (index++ < samples_per_pixel);
}

float RandomSampler::get_1d()
{
    return uniform_float();
}

Point2f RandomSampler::get_2d()
{
    return Point2f(uniform_float(), uniform_float());
}
