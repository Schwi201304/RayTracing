#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <cmath>
#include <limits>
#include <memory>
#include "image.h"
#include <cstdlib>
// Common Headers

#include "Smath.h"
#include "ray.h"
#include "Srandom.h"

// Usings
namespace schwi {
	using std::shared_ptr;
	using std::make_shared;
	using std::sqrt;
	using BYTE = unsigned char;
	// Constants

	const double infinity = std::numeric_limits<double>::infinity();

	// Utility Functions

	inline double degrees_to_radians(double degrees) {
		return degrees * Pi / 180.0;
	}

	inline Color c2Color(color c) {
		return Color((BYTE)(c[0] * 256), (BYTE)(c[1] * 256), (BYTE)(c[2] * 256));
	}

	inline color C2color(Color c) {
		return color(c[0] / 256.0, c[1] / 256.0, c[2] / 256.0);
	}

	inline double clamp(double x, double min, double max) {
		if (x < min) return min;
		if (x > max) return max;
		return x;
	}

	inline Color color_Correct(color c, int samples_per_pixel) {
		// Replace NaN components with zero. See explanation in Ray Tracing: The Rest of Your Life.
		if (c[0] != c[0]) c[0] = 0.0;
		if (c[1] != c[1]) c[1] = 0.0;
		if (c[2] != c[2]) c[2] = 0.0;

		// Divide the color by the number of samples and gamma-correct for gamma=2.0.
		c[0] = clamp(sqrt(c[0] / samples_per_pixel), 0, 0.999);
		c[1] = clamp(sqrt(c[1] / samples_per_pixel), 0, 0.999);
		c[2] = clamp(sqrt(c[2] / samples_per_pixel), 0, 0.999);
		return c2Color(c);
	}

	inline int random_int(int min, int max) {
		// Returns a random integer in [min,max].
		return static_cast<int>(random_double(min, max + 1));
	}
}
#endif