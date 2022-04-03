#ifndef SCHWI_TRANSFORM_H
#define SCHWI_TRANSFORM_H

#include "schwi.h"
#include "matrix.h"
#include "vector3.h"
#include "bound.h"
#include "quaternion.h"
#include "interaction.h"

namespace schwi {
	class Transform {

	private:
		// Transform Private Data
		Matrix4d m, mInv;
		friend class AnimatedTransform;
		friend class Quaternion;
	public:
		// Transform Public Methods
		Transform() {}
		Transform(const double mat[4][4]) {
			m = Matrix4d({ mat[0][0], mat[0][1], mat[0][2], mat[0][3], mat[1][0],
				mat[1][1], mat[1][2], mat[1][3], mat[2][0], mat[2][1],
				mat[2][2], mat[2][3], mat[3][0], mat[3][1], mat[3][2],
				mat[3][3] });
			mInv = m.inverse();
		}
		Transform(const Matrix4d& m) : m(m), mInv(m.inverse()) {}
		Transform(const Matrix4d& m, const Matrix4d& mInv) : m(m), mInv(mInv) {}
		void Print(FILE* f) const;
		friend Transform Inverse(const Transform& t) {
			return Transform(t.mInv, t.m);
		}
		friend Transform Transpose(const Transform& t) {
			return Transform(t.m.transpose(), t.mInv.transpose());
		}
		bool operator==(const Transform& t) const {
			return t.m == m && t.mInv == mInv;
		}
		bool operator!=(const Transform& t) const {
			return t.m != m || t.mInv != mInv;
		}
		bool operator<(const Transform& t2) const {
			for (int i = 0; i < 4; ++i)
				for (int j = 0; j < 4; ++j) {
					if (m(i, j) < t2.m(i, j)) return true;
					if (m(i, j) > t2.m(i, j)) return false;
				}
			return false;
		}
		bool IsIdentity() const {
			return (m(0, 0) == 1.f && m(0, 1) == 0. && m(0, 2) == 0. &&
				m(0, 3) == 0. && m(1, 0) == 0. && m(1, 1) == 1. &&
				m(1, 2) == 0. && m(1, 3) == 0. && m(2, 0) == 0. &&
				m(2, 1) == 0. && m(2, 2) == 1. && m(2, 3) == 0. &&
				m(3, 0) == 0. && m(3, 1) == 0. && m(3, 2) == 0. &&
				m(3, 3) == 1.);
		}
		const Matrix4d& GetMatrix() const { return m; }
		const Matrix4d& GetInverseMatrix() const { return mInv; }
		bool HasScale() const {
			double la2 = (*this)(Vector3d(1, 0, 0)).LengthSquared();
			double lb2 = (*this)(Vector3d(0, 1, 0)).LengthSquared();
			double lc2 = (*this)(Vector3d(0, 0, 1)).LengthSquared();
#define NOT_ONE(x) ((x) < .999f || (x) > 1.001f)
			return (NOT_ONE(la2) || NOT_ONE(lb2) || NOT_ONE(lc2));
#undef NOT_ONE
		}
		template <typename T>
		inline Point3<T> operator()(const Point3<T>& p) const;
		template <typename T>
		inline Vector3<T> operator()(const Vector3<T>& v) const;
		template <typename T>
		inline Normal3<T> operator()(const Normal3<T>&) const;
		inline Ray operator()(const Ray& r) const;
		inline RayDifferential operator()(const RayDifferential& r) const;
		Bounds3d operator()(const Bounds3d& b) const;
		Transform operator*(const Transform& t2) const;
		bool SwapsHandedness() const;
		SurfaceInteraction operator()(const SurfaceInteraction& si) const;
		template <typename T>
		inline Point3<T> operator()(const Point3<T>& pt,
			Vector3<T>* absError) const;
		template <typename T>
		inline Point3<T> operator()(const Point3<T>& p, const Vector3<T>& pError,
			Vector3<T>* pTransError) const;
		template <typename T>
		inline Vector3<T> operator()(const Vector3<T>& v,
			Vector3<T>* vTransError) const;
		template <typename T>
		inline Vector3<T> operator()(const Vector3<T>& v, const Vector3<T>& vError,
			Vector3<T>* vTransError) const;
		inline Ray operator()(const Ray& r, Vector3d* oError,
			Vector3d* dError) const;
		inline Ray operator()(const Ray& r, const Vector3d& oErrorIn,
			const Vector3d& dErrorIn, Vector3d* oErrorOut,
			Vector3d* dErrorOut) const;

		friend std::ostream& operator<<(std::ostream& os, const Transform& t) {
			os << "t=" << t.m << ", inv=" << t.mInv;
			return os;
		}

	};

    Transform Translate(const Vector3d& delta);
    Transform Scale(double x, double y, double z);
    Transform RotateX(double theta);
    Transform RotateY(double theta);
    Transform RotateZ(double theta);
    Transform Rotate(double theta, const Vector3d& axis);
    Transform LookAt(const Point3d& pos, const Point3d& look, const Vector3d& up);
    Transform Orthographic(double znear, double zfar);
    Transform Perspective(double fov, double znear, double zfar);
    bool SolveLinearSystem2x2(const double A[2][2], const double B[2], double* x0,double* x1);

    // Transform Inline Functions
    template <typename T>
    inline Point3<T> Transform::operator()(const Point3<T>& p) const {
        T x = p.x, y = p.y, z = p.z;
        T xp = m(0,0) * x + m(0,1) * y + m(0,2) * z + m(0,3);
        T yp = m(1,0) * x + m(1,1) * y + m(1,2) * z + m(1,3);
        T zp = m(2,0) * x + m(2,1) * y + m(2,2) * z + m(2,3);
        T wp = m(3,0) * x + m(3,1) * y + m(3,2) * z + m(3,3);
        CHECK_NE(wp, 0);
        if (wp == 1)
            return Point3<T>(xp, yp, zp);
        else
            return Point3<T>(xp, yp, zp) / wp;
    }

    template <typename T>
    inline Vector3<T> Transform::operator()(const Vector3<T>& v) const {
        T x = v.x, y = v.y, z = v.z;
        return Vector3<T>(m(0,0) * x + m(0,1) * y + m(0,2) * z,
            m(1,0) * x + m(1,1) * y + m(1,2) * z,
            m(2,0) * x + m(2,1) * y + m(2,2) * z);
    }

    template <typename T>
    inline Normal3<T> Transform::operator()(const Normal3<T>& n) const {
        T x = n.x, y = n.y, z = n.z;
        return Normal3<T>(mInv(0,0) * x + mInv(1,0) * y + mInv(2,0) * z,
            mInv(0,1) * x + mInv(1,1) * y + mInv(2,1) * z,
            mInv(0,2) * x + mInv(1,2) * y + mInv(2,2) * z);
    }

    inline Ray Transform::operator()(const Ray& r) const {
        Vector3d oError;
        Point3d o = (*this)(r.o, &oError);
        Vector3d d = (*this)(r.d);
        double lengthSquared = d.LengthSquared();
        double tMax = r.tMax;
        if (lengthSquared > 0) {
            double dt = dot(abs(d), oError) / lengthSquared;
            o += d * dt;
            tMax -= dt;
        }
        return Ray(o, d, tMax, r.t, r.medium);
    }

    inline RayDifferential Transform::operator()(const RayDifferential& r) const {
        Ray tr = (*this)(Ray(r));
        RayDifferential ret(tr.o, tr.d, tr.tMax, tr.t, tr.medium);
        ret.hasDifferentials = r.hasDifferentials;
        ret.rxOrigin = (*this)(r.rxOrigin);
        ret.ryOrigin = (*this)(r.ryOrigin);
        ret.rxDirection = (*this)(r.rxDirection);
        ret.ryDirection = (*this)(r.ryDirection);
        return ret;
    }

    template <typename T>
    inline Point3<T> Transform::operator()(const Point3<T>& p,
        Vector3<T>* pError) const {
        T x = p.x, y = p.y, z = p.z;
        // Compute transformed coordinates from point _pt_
        T xp = (m(0,0) * x + m(0,1) * y) + (m(0,2) * z + m(0,3));
        T yp = (m(1,0) * x + m(1,1) * y) + (m(1,2) * z + m(1,3));
        T zp = (m(2,0) * x + m(2,1) * y) + (m(2,2) * z + m(2,3));
        T wp = (m(3,0) * x + m(3,1) * y) + (m(3,2) * z + m(3,3));

        // Compute absolute error for transformed point
        T xabsSum = (std::abs(m(0, 0) * x) + std::abs(m(0, 1) * y) +
            std::abs(m(0, 2) * z) + std::abs(m(0, 3));
        T yabsSum = (std::abs(m(1, 0) * x) + std::abs(m(1, 1) * y) +
            std::abs(m(1, 2) * z) + std::abs(m(1, 3)));
        T zabsSum = (std::abs(m(2, 0) * x) + std::abs(m(2, 1) * y) +
            std::abs(m(2, 2) * z) + std::abs(m(3, 3));
        *pError = gamma(3) * Vector3<T>(xabsSum, yabsSum, zabsSum);
        //CHECK_NE(wp, 0);
        if (wp == 1)
            return Point3<T>(xp, yp, zp);
        else
            return Point3<T>(xp, yp, zp) / wp;
    }

    template <typename T>
    inline Point3<T> Transform::operator()(const Point3<T>& pt,
        const Vector3<T>& ptError,
        Vector3<T>* absError) const {
        T x = pt.x, y = pt.y, z = pt.z;
        T xp = (m(0,0) * x + m(0,1) * y) + (m(0,2) * z + m(0,3));
        T yp = (m(1,0) * x + m(1,1) * y) + (m(1,2) * z + m(1,3));
        T zp = (m(2,0) * x + m(2,1) * y) + (m(2,2) * z + m(2,3));
        T wp = (m(3,0) * x + m(3,1) * y) + (m(3,2) * z + m(3,3));
        absError->x =
            (gamma(3) + (T)1) *
            (std::abs(m(0,0)) * ptError.x + std::abs(m(0,1)) * ptError.y +
                std::abs(m(0,2)) * ptError.z) +
            gamma(3) * (std::abs(m(0,0) * x) + std::abs(m(0,1) * y) +
                std::abs(m(0,2) * z) + std::abs(m(0,3)));
        absError->y =
            (gamma(3) + (T)1) *
            (std::abs(m(1,0)) * ptError.x + std::abs(m(1,1)) * ptError.y +
                std::abs(m(1,2)) * ptError.z) +
            gamma(3) * (std::abs(m(1,0) * x) + std::abs(m(1,1) * y) +
                std::abs(m(1,2) * z) + std::abs(m(1,3)));
        absError->z =
            (gamma(3) + (T)1) *
            (std::abs(m(2,0)) * ptError.x + std::abs(m(2,1)) * ptError.y +
                std::abs(m(2,2)) * ptError.z) +
            gamma(3) * (std::abs(m(2,0) * x) + std::abs(m(2,1) * y) +
                std::abs(m(2,2) * z) + std::abs(m(2,3)));
        //CHECK_NE(wp, 0);
        if (wp == 1.)
            return Point3<T>(xp, yp, zp);
        else
            return Point3<T>(xp, yp, zp) / wp;
    }

    template <typename T>
    inline Vector3<T> Transform::operator()(const Vector3<T>& v,
        Vector3<T>* absError) const {
        T x = v.x, y = v.y, z = v.z;
        absError->x =
            gamma(3) * (std::abs(m(0,0) * v.x) + std::abs(m(0,1) * v.y) +
                std::abs(m(0,2) * v.z));
        absError->y =
            gamma(3) * (std::abs(m(1,0) * v.x) + std::abs(m(1,1) * v.y) +
                std::abs(m(1,2) * v.z));
        absError->z =
            gamma(3) * (std::abs(m(2,0) * v.x) + std::abs(m(2,1) * v.y) +
                std::abs(m(2,2) * v.z));
        return Vector3<T>(m(0,0) * x + m(0,1) * y + m(0,2) * z,
            m(1,0) * x + m(1,1) * y + m(1,2) * z,
            m(2,0) * x + m(2,1) * y + m(2,2) * z);
    }

    template <typename T>
    inline Vector3<T> Transform::operator()(const Vector3<T>& v,
        const Vector3<T>& vError,
        Vector3<T>* absError) const {
        T x = v.x, y = v.y, z = v.z;
        absError->x =
            (gamma(3) + (T)1) *
            (std::abs(m(0,0)) * vError.x + std::abs(m(0,1)) * vError.y +
                std::abs(m(0,2)) * vError.z) +
            gamma(3) * (std::abs(m(0,0) * v.x) + std::abs(m(0,1) * v.y) +
                std::abs(m(0,2) * v.z));
        absError->y =
            (gamma(3) + (T)1) *
            (std::abs(m(1,0)) * vError.x + std::abs(m(1,1)) * vError.y +
                std::abs(m(1,2)) * vError.z) +
            gamma(3) * (std::abs(m(1,0) * v.x) + std::abs(m(1,1) * v.y) +
                std::abs(m(1,2) * v.z));
        absError->z =
            (gamma(3) + (T)1) *
            (std::abs(m(2,0)) * vError.x + std::abs(m(2,1)) * vError.y +
                std::abs(m(2,2)) * vError.z) +
            gamma(3) * (std::abs(m(2,0) * v.x) + std::abs(m(2,1) * v.y) +
                std::abs(m(2,2) * v.z));
        return Vector3<T>(m(0,0) * x + m(0,1) * y + m(0,2) * z,
            m(1,0) * x + m(1,1) * y + m(1,2) * z,
            m(2,0) * x + m(2,1) * y + m(2,2) * z);
    }

    inline Ray Transform::operator()(const Ray& r, Vector3d* oError,
        Vector3d* dError) const {
        Point3d o = (*this)(r.o, oError);
        Vector3d d = (*this)(r.d, dError);
        double tMax = r.tMax;
        double lengthSquared = d.LengthSquared();
        if (lengthSquared > 0) {
            double dt = dot(abs(d), *oError) / lengthSquared;
            o += d * dt;
            //        tMax -= dt;
        }
        return Ray(o, d, tMax, r.t, r.medium);
    }

    inline Ray Transform::operator()(const Ray& r, const Vector3d& oErrorIn,
        const Vector3d& dErrorIn, Vector3d* oErrorOut,
        Vector3d* dErrorOut) const {
        Point3d o = (*this)(r.o, oErrorIn, oErrorOut);
        Vector3d d = (*this)(r.d, dErrorIn, dErrorOut);
        double tMax = r.tMax;
        double lengthSquared = d.LengthSquared();
        if (lengthSquared > 0) {
            double dt = dot(abs(d), *oErrorOut) / lengthSquared;
            o += d * dt;
            //        tMax -= dt;
        }
        return Ray(o, d, tMax, r.t, r.medium);
    }
/*
    // AnimatedTransform Declarations
    class AnimatedTransform {
    public:
        // AnimatedTransform Public Methods
        AnimatedTransform(const Transform* startTransform, double startTime,
            const Transform* endTransform, double endTime);
        static void Decompose(const Matrix4d& m, Vector3d* T, Quaternion* R,
            Matrix4d* S);
        void Interpolate(double time, Transform* t) const;
        Ray operator()(const Ray& r) const;
        RayDifferential operator()(const RayDifferential& r) const;
        Point3d operator()(double time, const Point3d& p) const;
        Vector3d operator()(double time, const Vector3d& v) const;
        bool HasScale() const {
            return startTransform->HasScale() || endTransform->HasScale();
        }
        Bounds3d MotionBounds(const Bounds3d& b) const;
        Bounds3d BoundPointMotion(const Point3d& p) const;

    private:
        // AnimatedTransform Private Data
        const Transform* startTransform, * endTransform;
        const double startTime, endTime;
        const bool actuallyAnimated;
        Vector3d T[2];
        Quaternion R[2];
        Matrix4d S[2];
        bool hasRotation;
        struct DerivativeTerm {
            DerivativeTerm() {}
            DerivativeTerm(double c, double x, double y, double z)
                : kc(c), kx(x), ky(y), kz(z) {}
            double kc, kx, ky, kz;
            double Eval(const Point3d& p) const {
                return kc + kx * p.x + ky * p.y + kz * p.z;
            }
        };
        DerivativeTerm c1[3], c2[3], c3[3], c4[3], c5[3];
    };
*/
}

#endif
