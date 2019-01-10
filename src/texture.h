#pragma once

#include "point.h"
#include "spectrum.h"

#include <vector>
#include <string>
#include <cstring>

struct Intersection;

RgbSpectrum *read_image(const std::string &path, Point2i *resolution);

template<typename T>
struct Texture
{
    virtual ~Texture() {}
    virtual T evaluate(const Intersection &its) const = 0;
};

template<typename T>
struct ConstantTexture : public Texture<T>
{
    T value;

    ConstantTexture(const T &value) : value(value) {}

    T evaluate(const Intersection &its) const override
    {
        return value;
    }
};

template<typename T>
struct Mipmap
{
    struct Level
    {
        Point2i resolution;
        T *data; // TODO SPEED: BlockedArray
    };

    Point2i resolution;
    std::vector<Level> pyramid;

    Mipmap(const Point2i &resolution, const T *data);

    T lookup(const Point2f &st, float width = 0) const;

    const T &texel(int level, int s, int t) const;
    T triangle(int level, const Point2f &st) const;
};

struct ResampleWeight
{
    int first_texel;
    float weight[4];
};

static float lanczos(float x)
{
    x = std::abs(x);
    if (x < 1e-5)
        return 1;
    return std::sin(PI * x) / (PI * x);
}

static ResampleWeight *resample_weights(int old_res, int new_res)
{
    assert(new_res >= old_res);

    ResampleWeight *wt = new ResampleWeight[new_res];
    float filter_width = 2;

    for (int i = 0; i < new_res; ++i)
    {
        float center = (i + 0.5) * (float)old_res / (float)new_res;
        wt[i].first_texel = std::floor((center - filter_width) + 0.5);

        for (int j = 0; j < 4; ++j)
        {
            float pos = wt[i].first_texel + j + 0.5;
            wt[i].weight[j] = lanczos((pos - center) / filter_width);
        }

        float inv_sum_wts = 1 / (wt[i].weight[0] + wt[i].weight[1] +
                                 wt[i].weight[2] + wt[i].weight[3]);

        for (int j = 0; j < 4; ++j)
            wt[i].weight[j] *= inv_sum_wts;
    }

    return wt;
}

template<typename T>
Mipmap<T>::Mipmap(const Point2i &res, const T *data)
    : resolution(res)
{
    T *resampled_image = nullptr;
    if (!is_power_of_2(resolution[0]) || !is_power_of_2(resolution[1]))
    {
        Point2i res_pow2(round_up_power_of_2(resolution[0]), round_up_power_of_2(resolution[1]));

        // Resample in s-direction.
        ResampleWeight *s_weights = resample_weights(resolution[0], res_pow2[0]);
        resampled_image = new T[res_pow2[0] * res_pow2[1]];
        for (int t = 0; t < resolution[1]; ++t)
        {
            for (int s = 0; s < res_pow2[0]; ++s)
            {
                resampled_image[t * res_pow2[0] + s] = 0;
                for (int j = 0; j < 4; ++j)
                {
                    int s_orig = s_weights[s].first_texel + j;
                    // TODO: wrap mode option; always clamping for now
                    s_orig = clamp(s_orig, 0, resolution[0] - 1);

                    if ((s_orig >= 0) && (s_orig < resolution[0]))
                        resampled_image[t * res_pow2[0] + s] += s_weights[s].weight[j] * data[t * resolution[0] + s_orig];
                }
            }
        }

        // TODO FIXME
        // Resample in t-direction.
#if 0
        ResampleWeight *t_weights = resample_weights(resolution[1], res_pow2[1]);
        for (int s = 0; s < res_pow2[0]; ++s)
        {
            for (int t = 0; t < res_pow2[1]; ++t)
            {
            }

            for (int t = 0; t < res_pow2[1]; ++t)
                resampled_image[t * res_pow2[0] + s] = clamp(?????);
        }
#endif

        resolution = res_pow2;
    }

    int levels_count = 1 + log2int(std::max(resolution[0], resolution[1]));
    pyramid.resize(levels_count);

    pyramid[0].resolution = resolution;
    pyramid[0].data = new T[resolution[0] * resolution[1]];
    std::memcpy(pyramid[0].data, resampled_image ? resampled_image : data, sizeof(T) * resolution[0] * resolution[1]);

    for (int i = 1; i < levels_count; ++i)
    {
        Level &l = pyramid[i];

        int s_res = std::max(1, pyramid[i - 1].resolution[0] / 2);
        int t_res = std::max(1, pyramid[i - 1].resolution[1] / 2);

        l.resolution = Point2i(s_res, t_res);
        l.data = new T[s_res * t_res];

        for (int t = 0; t < t_res; ++t)
        {
            for (int s = 0; s < s_res; ++s)
            {
                l.data[t * s_res + s] = 0.25 * (texel(i - 1, 2 * s,     2 * t) +
                                                texel(i - 1, 2 * s + 1, 2 * t) +
                                                texel(i - 1, 2 * s,     2 * t + 1) +
                                                texel(i - 1, 2 * s + 1, 2 * t + 1));
            }
        }
    }
}

template<typename T>
T Mipmap<T>::lookup(const Point2f &st, float width) const
{
    int levels_count = (int)pyramid.size();
    float level = levels_count - 1 + std::log2(std::max(width, (float)1e-8));

    if (level < 0)
        return triangle(0, st);

    if (level >= levels_count - 1)
        return texel(levels_count - 1, 0, 0);

    int i_level = std::floor(level);
    float delta = level - i_level;
    return lerp(delta, triangle(i_level, st), triangle(i_level + 1, st));
}

template<typename T>
const T &Mipmap<T>::texel(int level, int s, int t) const
{
    const Level &l = pyramid[level];

    // TODO: wrap mode option; always clamping for now
    s = clamp(s, 0, l.resolution[0] - 1);
    t = clamp(t, 0, l.resolution[1] - 1);

    return l.data[t * l.resolution[0] + s];
}

template<typename T>
T Mipmap<T>::triangle(int level, const Point2f &st) const
{
    int levels_count = (int)pyramid.size();
    level = clamp(level, 0, levels_count - 1);

    float s = st[0] * pyramid[level].resolution[0] - 0.5;
    float t = st[1] * pyramid[level].resolution[1] - 0.5;

    int s0 = std::floor(s);
    int t0 = std::floor(t);
    float ds = s - s0;
    float dt = t - t0;

    return (1 - ds) * (1 - dt) * texel(level, s0,     t0) +
           (1 - ds) * dt       * texel(level, s0,     t0 + 1) +
           ds       * (1 - dt) * texel(level, s0 + 1, t0) +
           ds       * dt       * texel(level, s0 + 1, t0 + 1);
}
