#pragma once
#include "hittable.h"
#include "point3.h"
#include "pdf.h"

namespace schwi {
	class sphere : public hittable {
	public:
		Point3d center;
		double radius;
		shared_ptr<material> mat_ptr;

	public:
		sphere(Point3d cen, double r, shared_ptr<material> m)
			: center(cen), radius(r), mat_ptr(m) {};

		virtual bool hit(
			const Ray& r, double t_min, double t_max, hit_record& rec) const override;
		virtual bool bounding_box(
			double time0, double time1, aabb& output_box) const override;
		virtual double pdf_value(const Point3d& o, const Vector3d& v) const override;
		virtual Vector3d random(const Point3d& o) const override;
	private:
		static void get_sphere_uv(const Vector3d& p, double& u, double& v) {
			// p: a given point on the sphere of radius one, centered at the origin.
			// u: returned value [0,1] of angle around the Y axis from X=-1.
			// v: returned value [0,1] of angle from Y=-1 to Y=+1.
			//     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
			//     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
			//     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

			auto theta = acos(-p.y);
			auto phi = atan2(-p.z, p.x) + Pi;

			u = phi *Inv2Pi;
			v = theta *InvPi;
		}
	};

	inline bool sphere::hit(const Ray& r, double t_min, double t_max, hit_record& rec) const {
		Vector3d oc = r.origin() - center;
		auto a = r.direction().LengthSquared();
		auto half_b = dot(oc, r.direction());
		auto c = oc.LengthSquared() - radius * radius;

		auto discriminant = half_b * half_b - a * c;
		if (discriminant < 0) return false;
		auto sqrtd = sqrt(discriminant);

		// Find the nearest root that lies in the acceptable range.
		auto root = (-half_b - sqrtd) / a;
		if (root < t_min || t_max < root) {
			root = (-half_b + sqrtd) / a;
			if (root < t_min || t_max < root)
				return false;
		}

		rec.t = root;
		rec.p = r.at(rec.t);
		Vector3d outward_normal = (rec.p - center) / radius;
		rec.set_face_normal(r, outward_normal);
		get_sphere_uv(outward_normal, rec.u, rec.v);
		rec.mat_ptr = mat_ptr;

		return true;
	}

	inline bool sphere::bounding_box(double time0, double time1, aabb& output_box) const {
		output_box = aabb(
			center - Vector3d(radius, radius, radius),
			center + Vector3d(radius, radius, radius));
		return true;
	}

	inline double sphere::pdf_value(const Point3d& o, const Vector3d& v) const {
		hit_record rec;
		if (!this->hit(Ray(o, v), 0.001, Infinity, rec))
			return 0;

		auto cos_theta_max = sqrt(1 - radius * radius / (center - o).LengthSquared());
		auto solid_angle = 2 * Pi * (1 - cos_theta_max);

		return  1 / solid_angle;
	}

	inline Vector3d sphere::random(const Point3d& o) const {
		Vector3d direction = center - o;
		auto distance_squared = direction.LengthSquared();
		onb uvw;
		uvw.build_from_w(direction);
		return uvw.local(random_to_sphere(radius, distance_squared));
	}
}