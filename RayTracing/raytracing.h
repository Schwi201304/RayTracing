#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <cmath>
#include <limits>
#include <memory>
#include "image.h"
#include <cstdlib>
// Common Headers

#include "ray.h"

// Usings

using std::shared_ptr;
using std::make_shared;
using std::sqrt;
using schwi::Color;
using schwi::Image;

// Constants

const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// Utility Functions

inline double degrees_to_radians(double degrees) {
	return degrees * pi / 180.0;
}

inline Color c2Color(color c) {
	return Color((schwi::BYTE)c[0] * 256, (schwi::BYTE)c[1] * 256, (schwi::BYTE)c[2] * 256);
}

inline double clamp(double x, double min, double max) {
	if (x < min) return min;
	if (x > max) return max;
	return x;
}

inline Color color_Correct(color c) {
	c[0] = clamp(sqrt(c[0]), 0, 0.999);
	c[1] = clamp(sqrt(c[1]), 0, 0.999);
	c[2] = clamp(sqrt(c[2]), 0, 0.999);
	return c2Color(c);
}

inline int random_int(int min, int max) {
	// Returns a random integer in [min,max].
	return static_cast<int>(random_double(min, max + 1));
}
#endif