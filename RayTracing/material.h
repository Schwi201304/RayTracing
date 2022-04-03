#ifndef MATERIAL_H
#define MATERIAL_H

#include "raytracing.h"
#include "hittable.h"
#include "texture.h"
#include "onb.h"
#include "pdf.h"
namespace schwi {
	struct hit_record;

	struct scatter_record {
		Ray specular_ray;
		bool is_specular;
		color attenuation;
		shared_ptr<pdf> pdf_ptr;
	};

	class material {
	public:
		virtual color emitted(
			const Ray& r_in, const hit_record& rec, double u, double v, const Point3d& p
		) const {
			return color(0, 0, 0);
		}

		virtual bool scatter(
			const Ray& r_in, const hit_record& rec, scatter_record& srec
		) const {
			return false;
		}

		virtual double scattering_pdf(
			const Ray& r_in, const hit_record& rec, const Ray& scattered
		) const {
			return 0;
		}
	};

	class lambertian : public material {
	public:
		shared_ptr<texture> albedo;

	public:
		lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {}
		lambertian(shared_ptr<texture> a) : albedo(a) {}

		virtual bool scatter(
			const Ray& r_in, const hit_record& rec, scatter_record& srec
		) const override {
			srec.is_specular = false;
			srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
			srec.pdf_ptr = make_shared<cosine_pdf>(rec.normal);
			return true;
		}

		double scattering_pdf(
			const Ray& r_in, const hit_record& rec, const Ray& scattered
		) const {
			auto cosine = dot(rec.normal, normalize(scattered.direction()));
			return cosine < 0 ? 0 : cosine *InvPi;
		}
	};

	class metal : public material {
	public:
		color albedo;
		double fuzz;

	public:
		metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}

		virtual bool scatter(
			const Ray& r_in, const hit_record& rec, scatter_record& srec
		) const override {
			Vector3d reflected = reflect(normalize(r_in.direction()), rec.normal);
			srec.specular_ray = Ray(rec.p, reflected + fuzz * random_in_unit_sphere());
			srec.attenuation = albedo;
			srec.is_specular = true;
			srec.pdf_ptr = 0;
			return true;
		}
	};

	class dielectric : public material {
	public:
		double ir; // Index of Refraction

	public:
		dielectric(double index_of_refraction) : ir(index_of_refraction) {}

		virtual bool scatter(
			const Ray& r_in, const hit_record& rec, scatter_record& srec
		) const override {
			srec.is_specular = true;
			srec.pdf_ptr = nullptr;
			srec.attenuation = color(1.0, 1.0, 1.0);
			double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

			Vector3d unit_direction = normalize(r_in.direction());

			double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
			double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

			bool cannot_refract = refraction_ratio * sin_theta > 1.0;
			Vector3d direction;

			if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
				direction = reflect(unit_direction, rec.normal);
			else
				direction = refract(unit_direction, rec.normal, refraction_ratio);

			srec.specular_ray = Ray(rec.p, direction, r_in.time());
			return true;
		}
	private:
		static double reflectance(double cosine, double ref_idx) {
			// Use Schlick's approximation for reflectance.
			auto r0 = (1 - ref_idx) / (1 + ref_idx);
			r0 = r0 * r0;
			return r0 + (1 - r0) * pow((1 - cosine), 5);
		}
	};

	class diffuse_light : public material {
	public:
		diffuse_light(shared_ptr<texture> a) : emit(a) {}
		diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}

		virtual bool scatter(
			const Ray& r_in, const hit_record& rec, scatter_record& srecs
		) const override {
			return false;
		}

		virtual color emitted(const Ray& r_in, const hit_record& rec, double u, double v, const Point3d& p) const override {
			if (rec.front_face)
				return emit->value(u, v, p);
			else
				return color(0, 0, 0);
		}

	public:
		shared_ptr<texture> emit;
	};

	class isotropic : public material {
	public:
		isotropic(color c) : albedo(make_shared<solid_color>(c)) {}
		isotropic(shared_ptr<texture> a) : albedo(a) {}

		virtual bool scatter(
			const Ray& r_in, const hit_record& rec, scatter_record& srec
		) const override {
			srec.is_specular = true;
			srec.specular_ray = Ray(rec.p, random_in_unit_sphere(), r_in.time());
			srec.attenuation = albedo->value(rec.u, rec.v, rec.p);
			return true;
		}

	public:
		shared_ptr<texture> albedo;
	};
}
#endif
