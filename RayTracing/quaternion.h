#ifndef SCHWI_QUATERNION_H
#define SCHWI_QUATERNION_H

#include "vector3.h"
#include "matrix.h"
#include "transform.h"

namespace schwi {
	class Quaternion {
	public:
		Vector3d v;
		double w;
	public:
		Quaternion() :v(0, 0, 0), w(0) {}

		Quaternion& operator+=(const Quaternion& q) {
			v += q.v;
			w += q.w;
			return *this;
		}
		friend Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
			Quaternion ret = q1;
			return ret += q2;
		}
		Quaternion& operator-=(const Quaternion& q) {
			v -= q.v;
			w -= q.w;
			return *this;
		}
		Quaternion operator-() const {
			Quaternion ret;
			ret.v = -v;
			ret.w = -w;
			return ret;
		}
		friend Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
			Quaternion ret = q1;
			return ret -= q2;
		}
		Quaternion& operator*=(double f) {
			v *= f;
			w *= f;
			return *this;
		}
		Quaternion operator*(double f) const {
			Quaternion ret = *this;
			ret.v *= f;
			ret.w *= f;
			return ret;
		}
		Quaternion& operator/=(double f) {
			v /= f;
			w /= f;
			return *this;
		}
		Quaternion operator/(double f) const {
			Quaternion ret = *this;
			ret.v /= f;
			ret.w /= f;
			return ret;
		}
		Transform ToTransform() const {
			double xx = v.x * v.x, yy = v.y * v.y, zz = v.z * v.z;
			double xy = v.x * v.y, xz = v.x * v.z, yz = v.y * v.z;
			double wx = v.x * w, wy = v.y * w, wz = v.z * w;

			Matrix4d m;
			m(0,0) = 1 - 2 * (yy + zz);
			m(0,1) = 2 * (xy + wz);
			m(0,2) = 2 * (xz - wy);
			m(1,0) = 2 * (xy - wz);
			m(1,1) = 1 - 2 * (xx + zz);
			m(1,2) = 2 * (yz + wx);
			m(2,0) = 2 * (xz + wy);
			m(2,1) = 2 * (yz - wx);
			m(2,2) = 1 - 2 * (xx + yy);

			// Transpose since we are left-handed.  Ugh.
			return Transform(m.transpose(), m);
		}
		Quaternion(const Transform& t){
			const Matrix4d& m = t.m;
			double trace = m(0,0) + m(1,1) + m(2,2);
			if (trace > 0.f) {
				// Compute w from matrix trace, then xyz
				// 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
				double s = std::sqrt(trace + 1.0f);
				w = s / 2.0f;
				s = 0.5f / s;
				v.x = (m(2,1) - m(1,2)) * s;
				v.y = (m(0,2) - m(2,0)) * s;
				v.z = (m(1,0) - m(0,1)) * s;
			}
			else {
				// Compute largest of $x$, $y$, or $z$, then remaining components
				const int nxt[3] = { 1, 2, 0 };
				double q[3];
				int i = 0;
				if (m(1,1) > m(0,0)) i = 1;
				if (m(2,2) > m(i,i)) i = 2;
				int j = nxt[i];
				int k = nxt[j];
				double s = std::sqrt((m(i,i) - (m(j,j) + m(k,k))) + 1.0);
				q[i] = s * 0.5f;
				if (s != 0.f) s = 0.5f / s;
				w = (m(k,j) - m(j,k)) * s;
				q[j] = (m(j,i) + m(i,j)) * s;
				q[k] = (m(k,i) + m(i,k)) * s;
				v.x = q[0];
				v.y = q[1];
				v.z = q[2];
			}
		}

		friend std::ostream& operator<<(std::ostream& os, const Quaternion& q) {
			os << q.v.x << ' ' << q.v.y << ' ' << q.v.z << ' ' << q.w;
			return os;
		}
	};

	Quaternion Slerp(double t, const Quaternion& q1, const Quaternion& q2);

	// Quaternion Inline Functions
	inline Quaternion operator*(double f, const Quaternion& q) { return q * f; }

	inline double dot(const Quaternion& q1, const Quaternion& q2) {
		return dot(q1.v, q2.v) + q1.w * q2.w;
	}

	inline Quaternion Normalize(const Quaternion& q) {
		return q / std::sqrt(dot(q, q));
	}

}

#endif
