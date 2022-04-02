#ifndef RAY_H
#define RAY_H

#include "vec3.h"
namespace schwi {
    class ray {
    public:
        point3 orig;
        Vector3d dir;
        double tm;

    public:
        ray() {}
        ray(const point3& origin, const Vector3d& direction, double time = 0.0)
            : orig(origin), dir(direction), tm(time)
        {}

        point3 origin() const { return orig; }
        Vector3d direction() const { return dir; }
        double time() const { return tm; }

        point3 at(double t) const {
            return orig + t * dir;
        }
    };
}
#endif