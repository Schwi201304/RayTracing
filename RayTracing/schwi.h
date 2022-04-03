#ifndef SCHWI_GLOBAL_H
#define SCHWI_GLOBAL_H


// Global Include Files
#include <algorithm>
#include <assert.h>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

//Warning Disable
#pragma warning(disable : 4305)
#pragma warning(disable : 4244)
#pragma warning(disable : 4843)
#pragma warning(disable : 4267)
#pragma warning(disable : 4838)
#pragma warning(disable : 4996)

namespace schwi {
	// Global Forward Declarations
	template <typename T> class vec2;
	template <typename T> class vec3;
	template <typename T> class vec4;
	template <typename T> class Vector2;
	template <typename T> class Vector3;
	template <typename T> class Point2;
	template <typename T> class Point3;
	template <typename T> class Normal3;
	template <typename T> class Bounds3;
	class Shape;
	class Quaternion;
	class Transform;
	class Ray;
	class RayDifferential;
	class SurfaceInteraction;

	//Global Constants
	constexpr double Pi = 3.14159265358979323846;
	constexpr double InvPi = 0.31830988618379067154;
	constexpr double Inv2Pi = 0.15915494309189533577;
	constexpr double Inv4Pi = 0.07957747154594766788;
	constexpr double PiOver2 = 1.57079632679489661923;
	constexpr double PiOver4 = 0.78539816339744830961;
	constexpr double Sqrt2 = 1.41421356237309504880;
	constexpr double Log2 = 1.44269504088896322870;

	constexpr double Infinity = std::numeric_limits<double>::max();

	inline double Radians(double deg) {
		return (Pi / 180) * deg;
	}

	inline double Degrees(double rad) {
		return (180 * InvPi) * rad;
	}
}
#endif
