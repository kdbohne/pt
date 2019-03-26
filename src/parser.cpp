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
#include "sampler.h"
#include "material.h"
#include "texture.h"

#include <cstdio>

// The directory of the input path passed to parse_pbrt().
static std::string input_directory;

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
    std::vector<Parameter<bool>> bools;
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

    bool find_bool(const std::string &name, bool default_value) const
    {
        for (const Parameter<bool> &param : bools)
        {
            if ((param.name == name) && (param.values.size() == 1))
                return param.values[0];
        }

        return default_value;
    }

    std::string find_string(const std::string &name, const std::string &default_value) const
    {
        for (const Parameter<std::string> &param : strings)
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

    const std::vector<int> *find_ints(const std::string &name) const
    {
        for (const Parameter<int> &param : ints)
        {
            if (param.name == name)
                return &param.values;
        }

        return nullptr;
    }

    const std::vector<Point3f> *find_point3fs(const std::string &name) const
    {
        for (const Parameter<Point3f> &param : point3fs)
        {
            if (param.name == name)
                return &param.values;
        }

        return nullptr;
    }
};

struct GraphicsState
{
    // TODO: store more state, see GraphicsState in pbrt's api.cpp, line 201

    std::string material_type;
    ParameterList material_params;

    std::string area_light_type;
    ParameterList area_light_params;
};

static bool is_newline(char c)
{
    return (c == '\n') || (c == '\r');
}

static bool is_whitespace(char c)
{
    return (c == ' ') || (c == '\t') || is_newline(c);
}

static std::string read_file(const std::string &path)
{
    std::string source;

    FILE *file = std::fopen(path.c_str(), "r");
    if (!file)
    {
        error("Failed to open file: %s", path.c_str());
        return source;
    }

    std::fseek(file, 0, SEEK_END);
    int length = std::ftell(file);
    std::fseek(file, 0, SEEK_SET);

    source.resize(length);

    if (std::fread(&source[0], length, 1, file) != 1)
        fatal("Failed to read file: %s", path.c_str());

    std::fclose(file);

    return source;
}

struct File
{
    std::string source;
    char *c;
};

struct Parser
{
    std::vector<File> file_stack;

    Parser(const std::string &path)
    {
        include(path);
    }

    void include(const std::string &path)
    {
        file_stack.emplace_back();

        File &file = file_stack.back();
        file.source = read_file(path);
        file.c = &file.source[0];
    }

    char *c()
    {
        assert(!file_stack.empty());
        return file_stack.back().c;
    }

    void advance()
    {
        assert(!file_stack.empty());

        // TODO: column/line numbers
        ++file_stack.back().c;
    }

    void eat_whitespace()
    {
        while (true)
        {
            if (is_whitespace(*c()))
            {
                advance();
            }
            else if (*c() == '#')
            {
                while (!is_newline(*c()))
                    advance();

                advance();
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

        const char *start = c();

        while (*c() && !is_whitespace(*c()))
            advance();

        Token token;
        token.s = start;
        token.length = c() - token.s;

        return token;
    }

    double parse_number()
    {
        // TODO: robustness?
        eat_whitespace();
        return std::strtod(c(), &file_stack.back().c);
    }

    std::string parse_string()
    {
        eat_whitespace();

        if (*c() != '"')
        {
            // TODO: file/line/column info
            fatal("Expected quoted string; received \"%c\" instead.", *c());
            return "";
        }

        advance();

        const char *start = c();

        // TODO: handle escapes
        while (*c() != '"')
            advance();

        const char *end = c();

        advance();

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

        // TODO CLEANUP: reduce duplication
        // TODO CLEANUP: reorder types
        if (type == "integer")
        {
            params->ints.emplace_back();

            Parameter<int> &param = params->ints.back();
            param.name = name;

            if (*c() == '[')
            {
                advance();
                while (*c() != ']')
                {
                    param.values.push_back(parse_number());
                    eat_whitespace();
                }
                advance();
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

            if (*c() == '[')
            {
                advance();
                while (*c() != ']')
                {
                    param.values.push_back(parse_number());
                    eat_whitespace();
                }
                advance();
            }
            else
            {
                param.values.push_back(parse_number());
            }
        }
        else if (type == "bool")
        {
            params->bools.emplace_back();

            Parameter<bool> &param = params->bools.back();
            param.name = name;

            // TODO CLEANUP: reduce duplication
            if (*c() == '[')
            {
                advance();
                while (*c() != ']')
                {
                    std::string string = parse_string();
                    if (string == "true")
                        param.values.push_back(true);
                    else if (string == "false")
                        param.values.push_back(false);
                    else
                        error("Invalid boolean value: \"%s\"", string.c_str()); // TODO: file/line/column info

                    eat_whitespace();
                }
                advance();
            }
            else
            {
                std::string string = parse_string();
                if (string == "true")
                    param.values.push_back(true);
                else if (string == "false")
                    param.values.push_back(false);
                else
                    error("Invalid boolean value: \"%s\"", string.c_str()); // TODO: file/line/column info
            }
        }
        else if (type == "string")
        {
            params->strings.emplace_back();

            Parameter<std::string> &param = params->strings.back();
            param.name = name;

            if (*c() == '[')
            {
                advance();
                while (*c() != ']')
                {
                    param.values.push_back(parse_string());
                    eat_whitespace();
                }
                advance();
            }
            else
            {
                param.values.push_back(parse_string());
            }
        }
        else if ((type == "rgb") || (type == "color"))
        {
            params->spectrums.emplace_back();

            Parameter<Spectrum> &param = params->spectrums.back();
            param.name = name;

            if (*c() == '[')
            {
                advance();
                while (*c() != ']')
                {
                    float r = parse_number();
                    float g = parse_number();
                    float b = parse_number();
                    param.values.push_back(Spectrum(r, g, b));

                    eat_whitespace();
                }
                advance();
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

            if (*c() == '[')
            {
                advance();
                while (*c() != ']')
                {
                    float x = parse_number();
                    float y = parse_number();
                    float z = parse_number();
                    param.values.push_back(Point3f(x, y, z));

                    eat_whitespace();
                }
                advance();
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
            advance();
            float temperature = parse_number();
            float scale = parse_number();
            advance();
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
            // TODO: file/line/column info
            fatal("Unknown parameter type: \"%s\"", type.c_str());
        }
    }

    ParameterList parse_parameters()
    {
        ParameterList params;

        // Parameters are always double-quoted, so continue parsing parameters
        // until a non-quoted token is found.
        while (true)
        {
            eat_whitespace();

            if (*c() != '"')
                return params;

            parse_parameter(&params);
        }

        return params;
    }
};

static Material *make_material(const std::string &type, const ParameterList &params)
{
    // TODO MEMORY: creating a new material/textures per call
    // TODO: allow non-constant textures

    if (type == "matte")
    {
        Spectrum Kd = params.find_spectrum("Kd", Spectrum(0.5));
        return new MatteMaterial(new ConstantTexture<Spectrum>(Kd));
    }

    if (type == "plastic")
    {
        Spectrum Kd = params.find_spectrum("Kd", Spectrum(0.25));
        Spectrum Ks = params.find_spectrum("Ks", Spectrum(0.25));
        float roughness = params.find_float("roughness", 0.1);
        bool remap_roughness = params.find_bool("remaproughness", true);

        return new PlasticMaterial(new ConstantTexture<Spectrum>(Kd),
                                   new ConstantTexture<Spectrum>(Ks),
                                   new ConstantTexture<float>(roughness),
                                   remap_roughness);
    }

    if (type == "glass")
    {
        Spectrum Kr = params.find_spectrum("Kr", Spectrum(1));
        Spectrum Kt = params.find_spectrum("Kt", Spectrum(1));
        float eta = params.find_float("index", 1.5);
        float u_roughness = params.find_float("uroughness", 0);
        float v_roughness = params.find_float("vroughness", 0);
        bool remap_roughness = params.find_bool("remaproughness", true);

        return new GlassMaterial(new ConstantTexture<Spectrum>(Kr),
                                 new ConstantTexture<Spectrum>(Kt),
                                 new ConstantTexture<float>(u_roughness),
                                 new ConstantTexture<float>(v_roughness),
                                 new ConstantTexture<float>(eta),
                                 remap_roughness);
    }

    // TODO: file/line/column info
    error("Unknown material type: \"%s\"", type.c_str());
    return nullptr;
}

static std::vector<Geometry *> make_geometries(const std::string &type, const ParameterList &params, const Transform &transform)
{
    std::vector<Geometry *> geometries;

    if (type == "sphere")
    {
        // NOTE: only using "radius" parameter for now.
        float radius = params.find_float("radius", 1.0);

        Sphere *sphere = new Sphere(transform, inverse(transform), radius);
        geometries.push_back(sphere);
    }
    else if (type == "trianglemesh")
    {
        MeshData data;

        const std::vector<int> *indices = params.find_ints("indices");
        if (indices)
        {
            if (indices->size() % 3 != 0)
            {
                // TODO: file/line/column info
                error("Triangle mesh indices list contains incomplete triangle(s): %d indices", (int)indices->size());
            }
            else
            {
                for (int i = 0; i < (int)indices->size(); i += 3)
                {
                    data.tris.emplace_back();
                    TriangleData &tri = data.tris.back();
                    for (int j = 0; j < 3; ++j)
                    {
                        tri.pi[j]  = (*indices)[i + j];
                        tri.ni[j]  = (*indices)[i + j];
                        tri.uvi[j] = (*indices)[i + j];
                    }
                }
            }
        }
        else
        {
            // TODO: file/line/column info
            error("Triangle mesh does not have indices.%s", ""); // TODO FIXME HACK: 0-arg error()
        }

        const std::vector<Point3f> *positions = params.find_point3fs("P");
        if (positions)
        {
            // TODO: check positions->size() against maximum index in indices list
            for (size_t i = 0; i < positions->size(); ++i)
                data.p.push_back((*positions)[i]);
        }
        else
        {
            // TODO: file/line/column info
            error("Triangle mesh does not have positions.%s", ""); // TODO FIXME HACK: 0-arg error()
        }

        Mesh *mesh = new Mesh(transform, inverse(transform), data);

        for (const TriangleData &t : mesh->data.tris)
        {
            Triangle *triangle = new Triangle(mesh->object_to_world, mesh->world_to_object, mesh, t);
            geometries.push_back(triangle);
        }
    }
    else if (type == "plymesh")
    {
        // TODO: alpha/shadowalpha texture parameters
        std::string filename = params.find_string("filename", "");

        MeshData data = load_ply(input_directory + "/" + filename);

        Mesh *mesh = new Mesh(transform, inverse(transform), data);

        for (const TriangleData &t : mesh->data.tris)
        {
            Triangle *triangle = new Triangle(mesh->object_to_world, mesh->world_to_object, mesh, t);
            geometries.push_back(triangle);
        }
    }
    else if (type == "loopsubdiv")
    {
        MeshData data;

        // TODO FIXME UNUSED: currently just using the given position data
        // without doing any Loop subdivision
        int levels = params.find_int("nlevels", 3);

        // TODO CLEANUP: this is copy/pasted from the trianglemesh case
        const std::vector<int> *indices = params.find_ints("indices");
        if (indices)
        {
            if (indices->size() % 3 != 0)
            {
                // TODO: file/line/column info
                error("Loop subdivision surface indices list contains incomplete triangle(s): %d indices", (int)indices->size());
            }
            else
            {
                for (int i = 0; i < (int)indices->size(); i += 3)
                {
                    data.tris.emplace_back();
                    TriangleData &tri = data.tris.back();
                    for (int j = 0; j < 3; ++j)
                    {
                        tri.pi[j]  = (*indices)[i + j];
                        tri.ni[j]  = (*indices)[i + j];
                        tri.uvi[j] = (*indices)[i + j];
                    }
                }
            }
        }
        else
        {
            // TODO: file/line/column info
            error("Loop subdivision surface does not have indices.%s", ""); // TODO FIXME HACK: 0-arg error()
        }

        const std::vector<Point3f> *positions = params.find_point3fs("P");
        if (positions)
        {
            // TODO: check positions->size() against maximum index in indices list
            for (size_t i = 0; i < positions->size(); ++i)
                data.p.push_back((*positions)[i]);
        }
        else
        {
            // TODO: file/line/column info
            error("Loop subdivision surface does not have positions.%s", ""); // TODO FIXME HACK: 0-arg error()
        }

        Mesh *mesh = new Mesh(transform, inverse(transform), data);

        for (const TriangleData &t : mesh->data.tris)
        {
            Triangle *triangle = new Triangle(mesh->object_to_world, mesh->world_to_object, mesh, t);
            geometries.push_back(triangle);
        }
    }
    else
    {
        // TODO: file/line/column info
        error("Unknown shape type \"%s\"; ignoring.", type.c_str());
    }

    return geometries;
}

bool parse_pbrt(const std::string &path, Scene *scene, Camera **camera, Integrator **integrator)
{
    size_t last_slash = path.find_last_of('/');
    if (last_slash == std::string::npos)
        input_directory = path;
    else
        input_directory = path.substr(0, last_slash);

    Parser parser(path);
    if (!parser.c())
        return false;

    int samples_per_pixel;

    Transform camera_to_world;
    float camera_fov = radians(90);
    Point2i film_resolution(640, 480);

    std::vector<GraphicsState> graphics_states;
    GraphicsState graphics_state;

    std::vector<Transform> transforms;
    Transform transform;

    while (true)
    {
        Token token = parser.next();
        while (token.length == 0)
        {
            if (parser.file_stack.empty())
                break;

            parser.file_stack.pop_back();
            token = parser.next();
        }

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

            transform = look_at(eye, target, up) * transform;
        }
        else if (token == "Camera")
        {
            std::string type = parser.parse_string();

            ParameterList params = parser.parse_parameters();
            camera_fov = radians(params.find_float("fov", 90));
            camera_to_world = inverse(transform);
        }
        else if (token == "Sampler")
        {
            // TODO FIXME UNUSED
            std::string type = parser.parse_string();

            ParameterList params = parser.parse_parameters();
            samples_per_pixel = params.find_int("pixelsamples", 1);
        }
        else if (token == "Integrator")
        {
            // TODO FIXME UNUSED
            std::string type = parser.parse_string();

            ParameterList params = parser.parse_parameters();
        }
        else if (token == "Film")
        {
            std::string type = parser.parse_string();

            ParameterList params = parser.parse_parameters();
            film_resolution.x = params.find_int("xresolution", 640);
            film_resolution.y = params.find_int("yresolution", 480);
        }
        else if (token == "WorldBegin")
        {
            graphics_states.clear();
            graphics_state = {};

            transforms.clear();
            transform = Transform();
        }
        else if (token == "WorldEnd")
        {
            break;
        }
        else if (token == "LightSource")
        {
            std::string type = parser.parse_string();

            ParameterList params = parser.parse_parameters();
            Spectrum scale = params.find_spectrum("scale", Spectrum(1, 1, 1));

            if (type == "point")
            {
                // TODO FIXME
                assert(false);
            }
            else if (type == "distant")
            {
                Spectrum L = params.find_spectrum("L", Spectrum(1, 1, 1));
                Point3f from = params.find_point3f("from", Point3f(0, 0, 0));
                Point3f to = params.find_point3f("to", Point3f(0, 0, 1));
                Vector3f w = normalize(from - to);

                DirectionalLight *light = new DirectionalLight(transform, L * scale, w);
                scene->lights.push_back(light);
            }
            else if (type == "infinite")
            {
                Spectrum L = params.find_spectrum("L", Spectrum(1, 1, 1));
                int samples_count = params.find_int("samples", 1);
                std::string relative_texture_path = params.find_string("mapname", "");
                std::string absolute_texture_path = input_directory + "/" + relative_texture_path;

                InfiniteAreaLight *light = new InfiniteAreaLight(transform, L * scale, samples_count, absolute_texture_path);
                scene->lights.push_back(light);
            }
            else
            {
                // TODO: file/line/column info
                error("Unknown light source type \"%s\"; ignoring.", type.c_str());
            }
        }
        else if (token == "AttributeBegin")
        {
            graphics_states.push_back(graphics_state);
            transforms.push_back(transform);
        }
        else if (token == "AttributeEnd")
        {
            assert(!graphics_states.empty());
            graphics_state = graphics_states.back();
            graphics_states.pop_back();

            assert(!transforms.empty());
            transform = transforms.back();
            transforms.pop_back();
        }
        else if (token == "Material")
        {
            std::string type = parser.parse_string();
            graphics_state.material_type = type;
            graphics_state.material_params = parser.parse_parameters();
        }
        else if (token == "Shape")
        {
            std::string type = parser.parse_string();
            ParameterList params = parser.parse_parameters();

            std::vector<Geometry *> geometries = make_geometries(type, params, transform);

            bool has_area_light = (graphics_state.area_light_type.length() > 0);
            for (Geometry *geometry : geometries)
            {
                AreaLight *area_light = nullptr;
                if (has_area_light)
                {
                    Spectrum L = graphics_state.area_light_params.find_spectrum("L", Spectrum(1, 1, 1));
                    bool two_sided = graphics_state.area_light_params.find_bool("twosided", false); // TODO UNUSED
                    int samples_count = graphics_state.area_light_params.find_int("samples", 1);

                    area_light = new DiffuseAreaLight(transform, samples_count, L, geometry);
                    scene->lights.push_back(area_light);
                }

                // TODO: avoid duplicating material per geometry
                Material *material = make_material(graphics_state.material_type, graphics_state.material_params);

                scene->entities.push_back(Entity(geometry, area_light, material));
            }
        }
        else if (token == "Texture")
        {
            // TODO FIXME UNUSED
            std::string name = parser.parse_string();
            std::string type = parser.parse_string();
            std::string class_ = parser.parse_string();

            ParameterList params = parser.parse_parameters();
        }
        else if (token == "Translate")
        {
            float v[3];
            for (int i = 0; i < 3; ++i)
                v[i] = parser.parse_number();

            transform = translate(Vector3f(v[0], v[1], v[2])) * transform;
        }
        else if (token == "Rotate")
        {
            float v[4];
            for (int i = 0; i < 4; ++i)
                v[i] = parser.parse_number();

            transform = rotate(v[0], Vector3f(v[1], v[2], v[3])) * transform;
        }
        else if (token == "Scale")
        {
            float v[3];
            for (int i = 0; i < 3; ++i)
                v[i] = parser.parse_number();

            transform = scale(v[0], v[1], v[2]) * transform;
        }
        else if (token == "CoordSysTransform")
        {
            // TODO FIXME UNUSED
            std::string name = parser.parse_string();
        }
        else if (token == "AreaLightSource")
        {
            std::string type = parser.parse_string();
            ParameterList params = parser.parse_parameters();

            if (type != "diffuse")
            {
                type = "";
                params = ParameterList();

                // TODO: file/line/column info
                error("Unsupported area light source type \"%s\"; ignoring.", type.c_str());
            }

            graphics_state.area_light_type = type;
            graphics_state.area_light_params = params;
        }
        else if (token == "Include")
        {
            std::string relative_path = parser.parse_string();
            std::string absolute_path = input_directory + "/" + relative_path;

            parser.include(absolute_path);
        }
        else if (token == "PixelFilter")
        {
            // TODO FIXME UNUSED
            std::string type = parser.parse_string();

            ParameterList params = parser.parse_parameters();
        }
        else
        {
            // TODO: file/line/column info
            fatal("Unknown identifier: \"%.*s\"", token.length, token.s);
        }
    }

    Film *film = new Film(film_resolution);
    *camera = new Camera(camera_to_world, camera_fov, film);

    RandomSampler *sampler = new RandomSampler(samples_per_pixel);
    *integrator = new WhittedIntegrator(5, sampler);

    return true;
}
