#ifndef PDF_H
#define PDF_H
#include "point3.h"
#include "onb.h"

namespace schwi {
    class pdf {
    public:
        virtual ~pdf() {}

        virtual double value(const Vector3d& direction) const = 0;
        virtual Vector3d generate() const = 0;
    };

    class cosine_pdf : public pdf {
    public:
        cosine_pdf(const Vector3d& w) { uvw.build_from_w(w); }

        virtual double value(const Vector3d& direction) const override {
            auto cosine = dot(normalize(direction), uvw.w);
            return (cosine <= 0) ? 0 : cosine*InvPi;
        }

        virtual Vector3d generate() const override {
            return uvw.local(random_cosine_direction());
        }

    public:
        onb uvw;
    };

    class hittable_pdf : public pdf {
    public:
        hittable_pdf(shared_ptr<hittable> p, const Point3d& origin) : ptr(p), o(origin) {}

        virtual double value(const Vector3d& direction) const override {
            return ptr->pdf_value(o, direction);
        }

        virtual Vector3d generate() const override {
            return ptr->random(o);
        }

    public:
        Point3d o;
        shared_ptr<hittable> ptr;
    };

    class mixture_pdf : public pdf {
    public:
        mixture_pdf(shared_ptr<pdf> p0, shared_ptr<pdf> p1) {
            p[0] = p0;
            p[1] = p1;
        }

        virtual double value(const Vector3d& direction) const override {
            return 0.5 * p[0]->value(direction) + 0.5 * p[1]->value(direction);
        }

        virtual Vector3d generate() const override {
            if (random_double() < 0.5)
                return p[0]->generate();
            else
                return p[1]->generate();
        }

    public:
        shared_ptr<pdf> p[2];
    };

    inline Vector3d random_to_sphere(double radius, double distance_squared) {
        auto r1 = random_double();
        auto r2 = random_double();
        auto z = 1 + r2 * (sqrt(1 - radius * radius / distance_squared) - 1);

        auto phi = 2 * Pi * r1;
        auto x = cos(phi) * sqrt(1 - z * z);
        auto y = sin(phi) * sqrt(1 - z * z);

        return Vector3d(x, y, z);
    }
}
#endif