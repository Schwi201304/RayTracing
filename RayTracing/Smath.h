#ifndef SCHWI_MATH_H
#define SCHWI_MATH_H

namespace schwi {
	constexpr double Pi			= 3.14159265358979323846;
	constexpr double InvPi		= 0.31830988618379067154;
	constexpr double Inv2Pi		= 0.15915494309189533577;
	constexpr double Inv4Pi		= 0.07957747154594766788;
	constexpr double PiOver2	= 1.57079632679489661923;
	constexpr double PiOver4	= 0.78539816339744830961;
	constexpr double Sqrt2		= 1.41421356237309504880;
	constexpr double Log2		= 1.44269504088896322870;

	inline double Radians(double deg) {
		return (Pi / 180) * deg;
	}

	inline double Degrees(double rad) {
		return (180 * InvPi) * rad;
	}
}
#endif
