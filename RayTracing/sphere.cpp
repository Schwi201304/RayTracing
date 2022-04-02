#include "sphere.h"
#include "pdf.h"
namespace schwi {
	bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
		Vector3d oc = r.origin() - center;
		auto a = r.direction().length_squared();
		auto half_b = dot(oc, r.direction());
		auto c = oc.length_squared() - radius * radius;

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

	bool sphere::bounding_box(double time0, double time1, aabb& output_box) const {
		output_box = aabb(
			center - Vector3d(radius, radius, radius),
			center + Vector3d(radius, radius, radius));
		return true;
	}

	double sphere::pdf_value(const point3& o, const Vector3d& v) const {
		hit_record rec;
		if (!this->hit(ray(o, v), 0.001, infinity, rec))
			return 0;

		auto cos_theta_max = sqrt(1 - radius * radius / (center - o).length_squared());
		auto solid_angle = 2 * Pi * (1 - cos_theta_max);

		return  1 / solid_angle;
	}

	Vector3d sphere::random(const point3& o) const {
		Vector3d direction = center - o;
		auto distance_squared = direction.length_squared();
		onb uvw;
		uvw.build_from_w(direction);
		return uvw.local(random_to_sphere(radius, distance_squared));
	}
}