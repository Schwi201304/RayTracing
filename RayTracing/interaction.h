#ifndef SCHWI_INTERACTION_H
#define SCHWI_INTERACTION_H

#include "point2.h"
#include "point3.h"
#include "normal3.h"
#include "medium.h"
#include "shape.h"
#include "ray.h"
#include "RayDifferential.h"

namespace schwi {
	class Interaction {
	public:
		Point3d p;
		double time;
		Vector3d pError;
		Vector3d wo;
		Normal3d n;
		MediumInterface* mediumInterface;
	public:
		Interaction(const Point3d& p, const Normal3d& n, const Vector3d& pError,
			const Vector3d& wo, double time,
			MediumInterface* mediumInterface)
			: p(p), time(time), pError(pError), wo(wo), n(n),
			mediumInterface(mediumInterface) { }

		bool IsSurfaceInteraction() const {
			return n != Normal3d();
		}
	};

	class SurfaceInteraction : public Interaction {
	public:
		Point2d uv;
		Vector3d dpdu, dpdv;
		Normal3d dndu, dndv;
		const Shape* shape = nullptr;
		struct {
			Normal3d n;
			Vector3d dpdu, dpdv;
			Normal3d dndu, dndv;
		} shading;

	public:
		SurfaceInteraction(const Point3d& p,
			const Vector3d& pError, const Point2d& uv, const Vector3d& wo,
			const Vector3d& dpdu, const Vector3d& dpdv,
			const Normal3d& dndu, const Normal3d& dndv,
			double time, const Shape* shape)
			: Interaction(p, Normal3d(normalize(cross(dpdu, dpdv))), pError, wo,
				time, nullptr),
			uv(uv), dpdu(dpdu), dpdv(dpdv), dndu(dndu), dndv(dndv),
			shape(shape) {
		}

		void SetShadingGeometry(const Vector3d& dpdus,
			const Vector3d& dpdvs, const Normal3d& dndus,
			const Normal3d& dndvs, bool orientationIsAuthoritative) {
			shading.n = Normal3d(cross(dpdus, dpdvs).normalized());
			if (shape && (shape->reverseOrientation ^
				shape->transformSwapsHandedness))
				shading.n = -shading.n;
			if (orientationIsAuthoritative)
				n = Faceforward(n, Vector3d(shading.n));
			else
				shading.n = Faceforward(shading.n, Vector3d(n));

			shading.dpdu = dpdus;
			shading.dpdv = dpdvs;
			shading.dndu = dndus;
			shading.dndv = dndvs;
		}
	};
}

#endif