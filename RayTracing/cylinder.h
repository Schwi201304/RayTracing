#pragma once

#include "hittable.h"
#include "point3.h"

namespace schwi {
	class Cylinder :public hittable {
	public:
		Point3d center;
		double radius, height;
		shared_ptr<material> mat_ptr;
		double zMin, zMax;

	public:
		Cylinder(Point3d c,double r, double h, shared_ptr<material> m)
			:center(c),radius(r),height(h), mat_ptr(m) {
			zMin = c.z - h / 2;
			zMax = c.z + h / 2;
		}

		virtual bool hit(
			const Ray& r, double t_min, double t_max, hit_record& rec)const override;
		virtual bool bounding_box(
			double time0, double time1, aabb& output_box)const override;

		void get_uv(const Point3d& p, double& u, double& v)const;

	};

	bool Cylinder::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
		auto d = r.d;
		auto o = r.o;

		auto a = d.x * d.x + d.y * d.y;
		auto b = 2 * (d.x * o.x + d.y * o.y);
		auto c = o.x * o.x + o.y * o.y - radius * radius;
		auto discriminant = b * b - 4 * a * c;
		auto sqrtd = sqrt(discriminant);

		auto t = (-b + sqrtd) / (2*a);
		if (t<t_min || t>t_max) {
			t = (-b + sqrtd) / (2 * a);
			if (t<t_min || t>t_max)\
				return false;
		}
		auto p = r.at(t);
		if (p.z<zMin || p.z>zMax)
			return false;

		rec.t = t;
		rec.p = p;
		Vector3d outward_normal(p.x-center.x,p.y-center.y,0);
		outward_normal = outward_normal.normalized();
		rec.set_face_normal(r, outward_normal);
		get_uv(p, rec.u, rec.v);
		rec.mat_ptr = mat_ptr;

		return true;
	}

	bool Cylinder::bounding_box(
		double time0, double time1, aabb& output_box)const {
		output_box = aabb(
			center+Point3d(-radius, -radius, -height/2),
			center+Point3d(radius, radius, height/2)
		);
		return true;
	}

	void Cylinder::get_uv(const Point3d& p, double& u, double& v)const {
		auto phi = atan2(p.x, p.y) + Pi;

		u = phi * Inv2Pi;
		v = (p.z - zMin) / height;
	}

}