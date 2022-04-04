#ifndef SCHWI_SHAPE_H
#define SCHWI_SHAPE_H

#include "transform.h"
namespace schwi {
	class Shape {
	public:
		const Transform* ObjectToWorld, * WorldToObject;
		const bool reverseOrientation;
		const bool transformSwapsHandedness;
	public:
		Shape(const Transform* ObjectToWorld,
			const Transform* WorldToObject, bool reverseOrientation)
			: ObjectToWorld(ObjectToWorld), WorldToObject(WorldToObject),
			reverseOrientation(reverseOrientation),
			transformSwapsHandedness(ObjectToWorld->SwapsHandedness()) {
		}

		~Shape(){}

		virtual Bounds3d ObjectBound() const = 0;

		Bounds3d WorldBound() const {
			return (*ObjectToWorld)(ObjectBound());
		}
	};
}

#endif
