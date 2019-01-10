#include "texture.h"

#include <OpenEXR/ImfRgba.h>
#include <OpenEXR/ImfRgbaFile.h>

static RgbSpectrum *read_image_exr(const std::string &path, int *width, int *height)
{
    Imf::RgbaInputFile file(path.c_str());
    Imath::Box2i dw = file.dataWindow();

    *width  = dw.max.x - dw.min.x + 1;
    *height = dw.max.y - dw.min.y + 1;

    std::vector<Imf::Rgba> pixels(*width * *height);
    file.setFrameBuffer(&pixels[0] - dw.min.x - dw.min.y * *width, 1, *width);
    file.readPixels(dw.min.y, dw.max.y);

    RgbSpectrum *result = new RgbSpectrum[*width * *height];
    for (int i = 0; i < *width * *height; ++i)
        result[i] = RgbSpectrum(pixels[i].r, pixels[i].g, pixels[i].b);

    return result;
}

RgbSpectrum *read_image(const std::string &path, Point2i *resolution)
{
    // TODO FIXME: assuming input file is .exr
    return read_image_exr(path, &resolution->x, &resolution->y);
}
