#include "SRandom.h"
#include "schwi.h"
namespace schwi {
    double random_double() {
        // Returns a random real in [0,1).
        return rand() / (RAND_MAX + 1.0);
    }

    double random_double(double min, double max) {
        // Returns a random real in [min,max).
        return min + (max - min) * random_double();
    }

    Vector3d random() {
        return Vector3d(random_double(), random_double(), random_double());
    }

    Vector3d random(double min, double max) {
        return Vector3d(random_double(min, max), random_double(min, max), random_double(min, max));
    }

    Vector3d random_in_unit_sphere() {
        while (true) {
            auto p = random(-1, 1);
            if (p.LengthSquared() >= 1) continue;
            return p;
        }
    }

    Vector3d random_unit_vector() {
        return normalize(random_in_unit_sphere());
    }

    Vector3d random_in_hemisphere(const Vector3d& normal) {
        Vector3d in_unit_sphere = random_in_unit_sphere();
        if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
            return in_unit_sphere;
        else
            return -in_unit_sphere;
    }

    Vector3d reflect(const Vector3d& v, const Vector3d& n) {
        return v - 2 * dot(v, n) * n;
    }

    Vector3d refract(const Vector3d& uv, const Vector3d& n, double etai_over_etat) {
        auto cos_theta = fmin(dot(-uv, n), 1.0);
        Vector3d r_out_perp = etai_over_etat * (uv + cos_theta * n);
        Vector3d r_out_parallel = -std::sqrt(fabs(1.0 - r_out_perp.LengthSquared())) * n;
        return r_out_perp + r_out_parallel;
    }

    Vector3d random_in_unit_disk() {
        while (true) {
            auto p = Vector3d(random_double(-1, 1), random_double(-1, 1), 0);
            if (p.LengthSquared() >= 1) continue;
            return p;
        }
    }

    Vector3d random_cosine_direction() {
        auto r1 = random_double();
        auto r2 = random_double();
        auto z = std::sqrt(1 - r2);

        auto phi = 2 * Pi * r1;
        auto x = cos(phi) * std::sqrt(r2);
        auto y = sin(phi) * std::sqrt(r2);

        return Vector3d(x, y, z);
    }
}