#include "scene.h"
#include "pdf.h"
#include <iostream>

namespace schwi {
	color ray_color(
		const Ray& r, const color& background, const hittable& world,
		shared_ptr<hittable> lights, int depth
	) {
		hit_record rec;
		// If we've exceeded the ray bounce limit, no more light is gathered.
		if (depth <= 0)
			return color(0, 0, 0);

		// If the ray hits nothing, return the background color.
		if (!world.hit(r, 0.001, Infinity, rec))
			return background;

		scatter_record srec;
		color emitted = rec.mat_ptr->emitted(r, rec, rec.u, rec.v, rec.p);
		if (!rec.mat_ptr->scatter(r, rec, srec))
			return emitted;

		if (srec.is_specular) {
			return srec.attenuation
				* ray_color(srec.specular_ray, background, world, lights, depth - 1);
		}

		auto light_ptr = make_shared<hittable_pdf>(lights, rec.p);
		mixture_pdf p(light_ptr, srec.pdf_ptr);

		Ray scattered = Ray(rec.p, p.generate(), r.time());
		auto pdf_val = p.value(scattered.direction());

		return emitted
			+ srec.attenuation * rec.mat_ptr->scattering_pdf(r, rec, scattered)
			* ray_color(scattered, background, world, lights, depth - 1) / pdf_val;
	}
}

using namespace schwi;
int main() {

	// Image

	auto aspect_ratio = 16.0 / 9.0;
	int image_width = 400;
	int samples_per_pixel = 100;
	int max_depth = 50;
	color background(0, 0, 0);

	// World

	hittable_list world;

	Point3d lookfrom;
	Point3d lookat;
	auto vfov = 40.0;
	auto aperture = 0.0;

	switch (6) {
	case 1:
		world = random_scene();
		background = color(0.70, 0.80, 1.00);
		lookfrom = Point3d(13, 2, 3);
		lookat = Point3d(0, 0, 0);
		vfov = 20.0;
		aperture = 0.1;
		break;
	case 2:
		world = two_spheres();
		background = color(0.70, 0.80, 1.00);
		lookfrom = Point3d(13, 2, 3);
		lookat = Point3d(0, 0, 0);
		vfov = 20.0;
		break;
	case 3:
		world = two_perlin_spheres();
		background = color(0.70, 0.80, 1.00);
		lookfrom = Point3d(13, 2, 3);
		lookat = Point3d(0, 0, 0);
		vfov = 20.0;
		break;
	case 4:
		world = earth();
		background = color(0.70, 0.80, 1.00);
		lookfrom = Point3d(13, 2, 3);
		lookat = Point3d(0, 0, 0);
		vfov = 20.0;
		break;
	case 5:
		world = simple_light();
		background = color(0, 0, 0);
		lookfrom = Point3d(26, 3, 6);
		lookat = Point3d(0, 2, 0);
		vfov = 20.0;
		break;
	case 6:
		world = cornell_box();
		aspect_ratio = 1.0;
		image_width = 400;
		samples_per_pixel = 50;
		background = color(0, 0, 0);
		lookfrom = Point3d(278, 278, -800);
		lookat = Point3d(278, 278, 0);
		vfov = 40.0;
		break;
	case 7:
		world = cornell_smoke();
		aspect_ratio = 1.0;
		image_width = 600;
		samples_per_pixel = 200;
		lookfrom = Point3d(278, 278, -800);
		lookat = Point3d(278, 278, 0);
		vfov = 40.0;
		break;
	default:
	case 8:
		world = final_scene();
		aspect_ratio = 1.0;
		image_width = 800;
		samples_per_pixel = 1000;
		background = color(0, 0, 0);
		lookfrom = Point3d(478, 278, -600);
		lookat = Point3d(278, 278, 0);
		vfov = 40.0;
		break;
	}

	shared_ptr<hittable_list> lights=make_shared<hittable_list>();
	lights->add(make_shared<xz_rect>(213, 343, 227, 332, 554, make_shared<diffuse_light>(color(15, 15, 15))));
	lights->add(make_shared<sphere>(Point3d(190, 90, 190), 90, make_shared<dielectric>(1.5)));
	// Camera

	Vector3d vup(0, 1, 0);
	auto dist_to_focus = 10.0;
	int image_height = static_cast<int>(image_width / aspect_ratio);

	camera cam(lookfrom, lookat, vup, vfov, aspect_ratio, aperture, dist_to_focus, 0.0, 1.0);

	// Render
	Image img(image_width, image_height);
	std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

#pragma omp parallel for schedule(dynamic, 1) private(r)
	for (int j = image_height - 1; j >= 0; --j) {
		std::cerr << "\rScanlines remaining: " << j << ' ' << std::flush;
		for (int i = 0; i < image_width; ++i) {
			color pixel_color(0, 0, 0);
			for (int s = 0; s < samples_per_pixel; ++s) {
				auto u = (i + random_double()) / ((double)image_width - 1);
				auto v = (j + random_double()) / ((double)image_height - 1);
				Ray r = cam.get_ray(u, v);
				pixel_color += ray_color(r, background, world, lights, max_depth);
			}
			img.setColor(i, j, color_Correct(pixel_color , samples_per_pixel));
		}
	}
	img.write_file("out.png", true);

	std::cerr << "\nDone.\n";
}