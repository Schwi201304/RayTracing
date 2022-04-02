#ifndef SCHWI_VEC3_H
#define SCHWI_VEC3_H

#ifndef NDEBUG
#define NDEBUG
#endif

#include<assert.h>
#include <cmath>
#include <iostream>

using std::isnan;
using std::sqrt;
namespace schwi {
	template<typename T>
	class vec3 {
	public:
		T x, y, z;
	public:
		vec3() = default;
		vec3(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {
			assert(isnan(_x) || isnan(_y)||isnan(_z));
		}

		vec3<T> operator-() const { return vec3<T>(-x, -y, -z); }
		T operator[](int i) const {
			assert(i >= 0 && i < 3);
			if (i == 0)return x;
			if (i == 1)return y;
			return z;
		}

		T& operator[](int i) {
			assert(i >= 0 && i < 3);
			if (i == 0)return x;
			if (i == 1)return y;
			return z;
		}

		vec3<T>& operator+=(const vec3<T>& v) {
			x += v.x;
			y += v.y;
			z += v.z;
			return *this;
		}

		template<typename T1>
		vec3<T>& operator*=(const T1 t) {
			x *= t;
			y *= t;
			z *= t;
			return *this;
		}

		template<typename T1>
		vec3<T>& operator/=(const T1 t) {
			return *this *= 1 / t;
		}

		double length() const {
			return std::sqrt(length_squared());
		}

		double length_squared() const {
			return x * x + y * y + z * z;
		}

		bool near_zero() const {
			// Return true if the vector is close to zero in all dimensions.
			const auto s = 1e-8;
			return (fabs(x) < s) && (fabs(y) < s) && (fabs(z) < s);
		}
	};

	template<typename T>
	inline std::ostream& operator<<(std::ostream& out, const vec3<T>& v) {
		return out << v.x << ' ' << v.y << ' ' << v.z;
	}

	template<typename T>
	inline vec3<T> operator+(const vec3<T>& u, const vec3<T>& v) {
		return vec3<T>(u.x + v.x, u.y + v.y, u.z + v.z);
	}

	template<typename T>
	inline vec3<T> operator-(const vec3<T>& u, const vec3<T>& v) {
		return vec3<T>(u.x - v.x, u.y - v.y, u.z - v.z);
	}

	template<typename T>
	inline vec3<T> operator*(const vec3<T>& u, const vec3<T>& v) {
		return vec3<T>(u.x * v.x, u.y * v.y, u.z * v.z);
	}

	template<typename T>
	inline vec3<T> operator*(T t, const vec3<T>& v) {
		return vec3<T>(t * v.x ,t * v.y, t * v.z);
	}

	template<typename T>
	inline vec3<T> operator*(const vec3<T>& v, T t) {
		return t * v;
	}

	template<typename T>
	inline vec3<T> operator/(const vec3<T>& v, T t) {
		return (1 / t) * v;
	}

	template<typename T>
	inline double dot(const vec3<T>& u, const vec3<T>& v) {
		return u.x * v.x
			+ u.y * v.y
			+ u.z * v.z;
	}

	template<typename T>
	vec3<double> cross(const vec3<T>& u, const vec3<T>& v) {
		return vec3<double>(u.y * v.z - u.z * v.y,
			u.z * v.x - u.x * v.z,
			u.x * v.y - u.y * v.x);
	}


	template<typename T>
	vec3<double> unit_vector(vec3<T> v) {
		return v / v.length();
	}

	using Vector3d = vec3<double>;
	using point3 = Vector3d;
	using color = Vector3d;
}
#endif