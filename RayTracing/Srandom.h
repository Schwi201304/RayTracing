#ifndef SCHWI_RANDOM_H
#define SCHWI_RANDOM_H

#include "vec3.h"
namespace schwi {
	double random_double();
	double random_double(double min, double max);
	Vector3d random();
	Vector3d random(double min, double max);
	Vector3d random_in_unit_sphere();
	Vector3d random_unit_vector();
	Vector3d random_in_hemisphere(const Vector3d& normal);
	Vector3d reflect(const Vector3d& v, const Vector3d& n);
	Vector3d refract(const Vector3d& uv, const Vector3d& n, double etai_over_etat);
	Vector3d random_in_unit_disk();
	Vector3d random_cosine_direction();
}
#endif
