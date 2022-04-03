#ifndef SCHWI_RAYDIFFERENTIAL_H
#define SCHWI_RAYDIFFERENTIAL_H

#include "ray.h"
#include "point3.h"

namespace schwi {
	class RayDifferential : public Ray {
	public:
		bool hasDifferentials;
		Point3d rxOrigin, ryOrigin;
		Vector3d rxDirection, ryDirection;

	public:
		RayDifferential() { hasDifferentials = false; }
		RayDifferential(const Point3d& o, const Vector3d& d,
			double tMax = Infinity, double time = 0.f,
			const Medium* medium = nullptr)
			: Ray(o, d, tMax, time, medium) {
			hasDifferentials = false;
		}
		RayDifferential(const Ray& ray) : Ray(ray) {
			hasDifferentials = false;
		}

		void ScaleDifferentials(double s) {
			rxOrigin = o + (rxOrigin - o) * s;
			ryOrigin = o + (ryOrigin - o) * s;
			rxDirection = d + (rxDirection - d) * s;
			ryDirection = d + (ryDirection - d) * s;
		}
	};
}
#endif