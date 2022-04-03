#ifndef SCHWI_RAY_H
#define SCHWI_RAY_H

#include "point3.h"
#include "schwi.h"
#include "medium.h"
namespace schwi {
    class Ray {
    public:
        Point3d o;
        Vector3d d;
        double tMax;
        double t;
        const Medium* medium;

    public:
        Ray() : tMax(Infinity), t(0.f), medium(nullptr) { }
        Ray(const Point3d& o, const Vector3d& d, double tMax = Infinity,
            double t = 0.f, const Medium* medium = nullptr)
            : o(o), d(d), tMax(tMax), t(t), medium(medium) { }

        Point3d origin() const { return o; }
        Vector3d direction() const { return d; }
        double time() const { return t; }

        Point3d at(double t) const {
            return o + t * d;
        }
    };
}
#endif