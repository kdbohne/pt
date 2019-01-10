#pragma once

struct Intersection;

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
