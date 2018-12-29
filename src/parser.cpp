#include "parser.h"
#include "common.h"
#include "transform.h"
#include "math.h"
#include "film.h"
#include "camera.h"
#include "integrator.h"
#include "geometry.h"
#include "scene.h"
#include "light.h"

#include <cstdio>

struct Token
{
    const char *s;
    int length;

    bool operator==(const char *str) const
    {
        int i;
        for (i = 0; *str; ++i, ++str)
        {
            if (i >= length)
                return false;
            if (*str != s[i])
                return false;
        }

        return (i == length);
    }
};

template<typename T>
struct Parameter
{
    std::string name;
    std::vector<T> values;
};

struct ParameterList
{
    std::vector<Parameter<int>> ints;
    std::vector<Parameter<float>> floats;
    std::vector<Parameter<std::string>> strings;
    std::vector<Parameter<Spectrum>> spectrums;
    std::vector<Parameter<Point3f>> point3fs;
    std::vector<Parameter<std::string>> textures;

    // TODO: reduce duplication?
    float find_int(const std::string &name, int default_value) const
    {
        for (const Parameter<int> &param : ints)
        {
            if ((param.name == name) && (param.values.size() == 1))
                return param.values[0];
        }

        return default_value;
    }

    float find_float(const std::string &name, float default_value) const
    {
        for (const Parameter<float> &param : floats)
        {
            if ((param.name == name) && (param.values.size() == 1))
                return param.values[0];
        }

        return default_value;
    }

    Spectrum find_spectrum(const std::string &name, const Spectrum &default_value) const
    {
        for (const Parameter<Spectrum> &param : spectrums)
        {
            if ((param.name == name) && (param.values.size() == 1))
                return param.values[0];
        }

        return default_value;
    }

    Point3f find_point3f(const std::string &name, const Point3f &default_value) const
    {
        for (const Parameter<Point3f> &param : point3fs)
        {
            if ((param.name == name) && (param.values.size() == 1))
                return param.values[0];
        }

        return default_value;
    }
};

struct GraphicsState
{
    // TODO: store more state, see GraphicsState in pbrt's api.cpp, line 201
    std::string material;
};

static bool is_newline(char c)
{
    return (c == '\n') || (c == '\r');
}

static bool is_whitespace(char c)
{
    return (c == ' ') || (c == '\t') || is_newline(c);
}

struct Parser
{
    std::string source;
    char *c;

    Parser(const std::string &path)
    {
        FILE *file = std::fopen(path.c_str(), "r");
        if (!file)
        {
            error("Failed to open file: %s", path.c_str());
            c = nullptr;

            return;
        }

        std::fseek(file, 0, SEEK_END);
        int length = std::ftell(file);
        std::fseek(file, 0, SEEK_SET);

        source.resize(length);

        if (std::fread(&source[0], length, 1, file) != 1)
            fatal("Failed to read file: %s", path.c_str());

        std::fclose(file);

        c = &source[0];
    }

    void eat_whitespace()
    {
        while (true)
        {
            if (is_whitespace(*c))
            {
                // TODO: advance()
                ++c;
            }
            else if (*c == '#')
            {
                while (!is_newline(*c))
                    ++c; // TODO: advance()

                // TODO: advance()
                ++c;
            }
            else
            {
                break;
            }
        }
    }

    Token next()
    {
        eat_whitespace();

        const char *start = c;

        while (!is_whitespace(*c))
            ++c; // TODO: advance()

        Token token;
        token.s = start;
        token.length = c - token.s;

        return token;
    }

    double parse_number()
    {
        // TODO: robustness?
        eat_whitespace();
        return std::strtod(c, &c);
    }

    std::string parse_string()
    {
        eat_whitespace();

        if (*c != '"')
        {
            // TODO: line/column numbers
            fatal("Expected quoted string; received \"%c\" instead.", *c);
            return "";
        }

        // TODO: advance()
        ++c;

        const char *start = c;

        // TODO: handle escapes
        while (*c != '"')
            ++c; // TODO: advance()

        const char *end = c;

        // TODO: advance()
        ++c;

        return std::string(start, end - start);
    }

    void parse_parameter(ParameterList *params)
    {
        eat_whitespace();

        std::string decl = parse_string();

        size_t space = decl.find_first_of(' ');
        std::string type = decl.substr(0, space);
        std::string name = decl.substr(space + 1);

        eat_whitespace();

//        printf("Type: \"%s\", Name: \"%s\"\n", type.c_str(), name.c_str());

        // TODO CLEANUP: reduce duplication
        // TODO CLEANUP: reorder types
        if (type == "integer")
        {
            params->ints.emplace_back();

            Parameter<int> &param = params->ints.back();
            param.name = name;

            if (*c == '[')
            {
                // TODO: advance()
                ++c;
                while (*c != ']')
                    param.values.push_back(parse_number());
                // TODO: advance()
                ++c;
            }
            else
            {
                param.values.push_back(parse_number());
            }
        }
        else if (type == "float")
        {
            params->floats.emplace_back();

            Parameter<float> &param = params->floats.back();
            param.name = name;

            if (*c == '[')
            {
                // TODO: advance()
                ++c;
                while (*c != ']')
                    param.values.push_back(parse_number());
                // TODO: advance()
                ++c;
            }
            else
            {
                param.values.push_back(parse_number());
            }
        }
        else if (type == "string")
        {
            params->strings.emplace_back();

            Parameter<std::string> &param = params->strings.back();
            param.name = name;

            // TODO: does this need to handle arrays?
            param.values.push_back(parse_string());
        }
        else if (type == "rgb")
        {
            params->spectrums.emplace_back();

            Parameter<Spectrum> &param = params->spectrums.back();
            param.name = name;

            if (*c == '[')
            {
                // TODO: advance()
                ++c;
                while (*c != ']')
                {
                    float r = parse_number();
                    float g = parse_number();
                    float b = parse_number();
                    param.values.push_back(Spectrum(r, g, b));

                    eat_whitespace();
                }
                // TODO: advance()
                ++c;
            }
            else
            {
                float r = parse_number();
                float g = parse_number();
                float b = parse_number();
                param.values.push_back(Spectrum(r, g, b));
            }
        }
        // TODO: xyz
        else if ((type == "point") || (type == "point3"))
        {
            params->point3fs.emplace_back();

            Parameter<Point3f> &param = params->point3fs.back();
            param.name = name;

            if (*c == '[')
            {
                // TODO: advance()
                ++c;
                while (*c != ']')
                {
                    float x = parse_number();
                    float y = parse_number();
                    float z = parse_number();
                    param.values.push_back(Point3f(x, y, z));

                    eat_whitespace();
                }
                // TODO: advance()
                ++c;
            }
            else
            {
                float x = parse_number();
                float y = parse_number();
                float z = parse_number();
                param.values.push_back(Point3f(x, y, z));
            }
        }
        else if (type == "blackbody")
        {
            params->spectrums.emplace_back();

            Parameter<Spectrum> &param = params->spectrums.back();
            param.name = name;

            // TODO: does this need to handle arrays?
            ++c; // TODO: advance()
            float temperature = parse_number();
            float scale = parse_number();
            ++c; // TODO: advance()
            // TODO FIXME: convert blackbody emitter values to Spectrumb
//            param.values.push_back(Spectrum());
        }
        else if (type == "texture")
        {
            params->textures.emplace_back();

            Parameter<std::string> &param = params->textures.back();
            param.name = name;

            // TODO: does this need to handle arrays?
            param.values.push_back(parse_string());
        }
        else
        {
            // TODO: line/column numbers
            fatal("Unknown parameter type: \"%s\"", type.c_str());
        }
    }

    void parse_parameters(ParameterList *params)
    {
        // Parameters are always double-quoted, so continue parsing parameters
        // until a non-quoted token is found.
        while (true)
        {
            eat_whitespace();

            if (*c != '"')
                return;

            parse_parameter(params);
        }
    }
};

bool parse_pbrt(const std::string &path, Scene *scene, Camera **camera, Integrator **integrator)
{
    Parser parser(path);
    if (!parser.c)
        return false;

    Transform camera_to_world;
    float camera_fov = radians(90);
    Point2i film_resolution(640, 480);

    std::vector<GraphicsState> graphics_states;
    GraphicsState graphics_state;

    Transform transform;

    while (true)
    {
        Token token = parser.next();
        if (token.length == 0)
            break;

        // TODO SPEED: switch on first character of token
        // TODO CLEANUP: reorder/alphabetize
        if (token == "LookAt")
        {
            float v[9];
            for (int i = 0; i < 9; ++i)
                v[i] = parser.parse_number();

            Point3f eye(v[0], v[1], v[2]);
            Point3f target(v[3], v[4], v[5]);
            Vector3f up(v[6], v[7], v[8]);

            camera_to_world = look_at(eye, target, up);
        }
        else if (token == "Camera")
        {
            std::string type = parser.parse_string();

            ParameterList params;
            parser.parse_parameters(&params);

            camera_fov = radians(params.find_float("fov", 90));
        }
        else if (token == "Sampler")
        {
            // TODO UNUSED
            std::string type = parser.parse_string();

            ParameterList params;
            parser.parse_parameters(&params);
        }
        else if (token == "Integrator")
        {
            // TODO UNUSED
            std::string type = parser.parse_string();

            ParameterList params;
            parser.parse_parameters(&params);
        }
        else if (token == "Film")
        {
            std::string type = parser.parse_string();

            ParameterList params;
            parser.parse_parameters(&params);

            film_resolution.x = params.find_int("xresolution", 640);
            film_resolution.y = params.find_int("yresolution", 480);
        }
        else if (token == "WorldBegin")
        {
            // TODO: ?
        }
        else if (token == "WorldEnd")
        {
            break;
        }
        else if (token == "LightSource")
        {
            std::string type = parser.parse_string();

            ParameterList params;
            parser.parse_parameters(&params);

            if (type == "point")
            {
                // TODO FIXME
            }
            else if (type == "distant")
            {
                Spectrum L = params.find_spectrum("L", Spectrum(1, 1, 1));
                Point3f from = params.find_point3f("from", Point3f(0, 0, 0));
                Point3f to = params.find_point3f("to", Point3f(0, 0, 1));
                Vector3f w = normalize(from - to);

                DirectionalLight *light = new DirectionalLight(transform, L, w);
                scene->lights.push_back(light);
            }
            else
            {
                // TODO: line/column numbers
                error("Unknown light source type \"%s\"; ignoring.", type.c_str());
            }
        }
        else if (token == "AttributeBegin")
        {
            graphics_states.push_back(graphics_state);
        }
        else if (token == "AttributeEnd")
        {
            assert(graphics_states.size() > 0);

            graphics_state = graphics_states.back();
            graphics_states.pop_back();
        }
        else if (token == "Material")
        {
            std::string name = parser.parse_string();
            graphics_state.material = name;

            // TODO FIXME
            ParameterList params;
            parser.parse_parameters(&params);
        }
        else if (token == "Shape")
        {
            std::string type = parser.parse_string();

            ParameterList params;
            parser.parse_parameters(&params);

            if (type == "sphere")
            {
                // NOTE: only using "radius" parameter for now.
                float radius = params.find_float("radius", 1.0);

                // TODO: reduce duplication here; creating a new BSDF per object
                Bsdf *bsdf = nullptr;
                if (graphics_state.material == "matte")
                    bsdf = make_matte_material(params.find_spectrum("Kd", Spectrum(0.5, 0.5, 0.5)));

                Sphere *sphere = new Sphere(transform, inverse(transform), radius);
                scene->entities.push_back(Entity(sphere, nullptr, bsdf));
            }
            else
            {
                error("Unknown shape type \"%s\"; ignoring.\n", type.c_str());
            }
        }
        else if (token == "Texture")
        {
            // TODO FIXME UNUSED
            std::string name = parser.parse_string();
            std::string type = parser.parse_string();
            std::string class_ = parser.parse_string();

            ParameterList params;
            parser.parse_parameters(&params);
        }
        else if (token == "Translate")
        {
            float v[3];
            for (int i = 0; i < 3; ++i)
                v[i] = parser.parse_number();

            transform = translate(Vector3f(v[0], v[1], v[2])) * transform;
        }
        else
        {
            // TODO: line/column numbers
            fatal("Unknown identifier: \"%.*s\"\n", token.length, token.s);
        }
    }

    Film *film = new Film(film_resolution);
    *camera = new Camera(camera_to_world, camera_fov, film);

    *integrator = new WhittedIntegrator(5);

    return true;
}
