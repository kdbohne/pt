#include "film.h"
#include <fstream>

void Film::write_ppm(const std::string &path) const
{
    std::ofstream file(path);

    file << "P3\n";
    file << resolution.x << " " << resolution.y << "\n";
    file << "255" << "\n";

    for (int y = 0; y < resolution.y; ++y)
    {
        for (int x = 0; x < resolution.x; ++x)
        {
            // TODO: tone mapping
            const Pixel *p = &pixels[y * resolution.x + x];
            int r = (int)(p->rgb[0] * 255);
            int g = (int)(p->rgb[1] * 255);
            int b = (int)(p->rgb[2] * 255);

            file << r << " " << g << " " << b << " ";
        }

        file << "\n";
    }
}
