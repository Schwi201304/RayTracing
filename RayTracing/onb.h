#ifndef ONB_H
#define ONB_H

#include "vec3.h"
namespace schwi {
    class onb {
    public:
        onb() {}

        inline Vector3d operator[](int i) const { 
            if (i == 0)return u;
            if (i == 1)return v;
            return w;
        }

        Vector3d local(double a, double b, double c) const {
            return a * u + b * v + c * w;
        }

        Vector3d local(const Vector3d& a) const {
            return a.x * u + a.y * v + a.z * w;
        }

        void build_from_w(const Vector3d&);

    public:
        Vector3d u,v,w;
    };


    inline void onb::build_from_w(const Vector3d& n) {
        w = unit_vector(n);
        Vector3d a = (fabs(w.x) > 0.9) ? Vector3d(0, 1, 0) : Vector3d(1, 0, 0);
        v = unit_vector(cross(w, a));
        u = cross(w, v);
    }
}
#endif