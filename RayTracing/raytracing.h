#ifndef RAYTRACING_H
#define RAYTRACING_H

#include <cmath>
#include <limits>
#include <memory>
#include "image.h"

// Common Headers

#include "ray.h"
#include "vec3.h"

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
	return Color(c[0] * 255, c[1] * 255, c[2] * 255);
}


#endif