#ifndef TEXTURE_H
#define TEXTURE_H

#include "raytracing.h"
#include "perlin.h"
namespace schwi {
    class texture {
    public:
        virtual color value(double u, double v, const Point3d& p) const = 0;
    };

    class solid_color : public texture {
    public:
        solid_color() {}
        solid_color(color c) : color_value(c) {}

        solid_color(double red, double green, double blue)
            : solid_color(color(red, green, blue)) {}

        virtual color value(double u, double v, const Point3d& p) const override {
            return color_value;
        }

    private:
        color color_value;
    };

    class checker_texture : public texture {
    public:
        checker_texture() {}

        checker_texture(shared_ptr<texture> _even, shared_ptr<texture> _odd)
            : even(_even), odd(_odd) {}

        checker_texture(color c1, color c2)
            : even(make_shared<solid_color>(c1)), odd(make_shared<solid_color>(c2)) {}

        virtual color value(double u, double v, const Point3d& p) const override {
            auto sines = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z);
            if (sines < 0)
                return odd->value(u, v, p);
            else
                return even->value(u, v, p);
        }

    public:
        shared_ptr<texture> odd;
        shared_ptr<texture> even;
    };

    class noise_texture : public texture {
    public:
        noise_texture() {}
        noise_texture(double sc) : scale(sc) {}

        virtual color value(double u, double v, const Point3d& p) const override {
            return color(1, 1, 1) * 0.5 * (1 + sin(scale * p.z + 10 * noise.turb(p)));
        }

    public:
        perlin noise;
        double scale;
    };

    class image_texture : public texture {
    public:
        const static int bytes_per_pixel = 3;
        Image img;

        image_texture() {}

        image_texture(const char* filename) {
            auto components_per_pixel = bytes_per_pixel;

            if (!img.load_file(filename)) {
                std::cerr << "Texture image load fail" << std::endl;
            }
        }

        virtual color value(double u, double v, const Point3d& p) const override {
            // If we have no texture data, then return solid cyan as a debugging aid.
            if (img.empty())
                return color(0, 1, 1);

            // Clamp input texture coordinates to [0,1] x [1,0]
            u = clamp(u, 0.0, 1.0);
            v = 1.0 - clamp(v, 0.0, 1.0);  // Flip V to image coordinates
            return C2color(img.getColor(u, v));
        }
    };
}
#endif