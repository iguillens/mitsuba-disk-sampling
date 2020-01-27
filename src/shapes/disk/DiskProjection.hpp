#ifndef SRC_DISKPROJECTION_HPP_
#define SRC_DISKPROJECTION_HPP_

#include <cstdint>
#include <cmath>
#include <iomanip>

#include "vector3.hpp"

using v3l::vector3;

#define BOOST_MATH_MAX_SERIES_ITERATION_POLICY 100
#include <boost/math/special_functions/ellint_3.hpp>

using boost::math::ellint_3;

template<class T> inline
T sign(const T& x)
{
	return T(std::copysign(T(1), x));
}

template<class T> inline
T clamp(const T& x, const T& min, const T& max)
{
	return T(std::min(max, std::max(min, x)));
}

/*
 * Spherical ellipse constants
 */
template<class T>
struct SphericalEllipse
{
	// Spherical ellipse params
	T a, b;
	// Elliptical integral constants
	T c, n, k;
	// Precalculated constants
	T phi_coeff, t_coeff/*,asqrtb, bsqrta*/, mult;
	// Reprojection constants
	T yhp_norm, dhp_norm;
};

/*
 * Spherical wedge constants
 */
template<class T>
struct SphericalWedge
{
	// Spherical wedge aperture
	T beta;
	// Reprojection constants
	T yhp_norm, dhp_norm;
};

/*
 * Spherical circle constants
 */
template<class T>
struct SphericalCircle
{
	// Spherical circle angle's sine
	T a;
};

/*
 * A class for sampling oriented disks from any orientation
 */
template<class T>
class DiskProjection
{
public:
	// Type of projection
	enum { ELLIPSE, WEDGE, CIRCLE } type;

	// Data required to generate samples
	union {
		struct SphericalEllipse<T> ellipse;
		struct SphericalWedge<T>   wedge;
		struct SphericalCircle<T>  circle;
	} projection;

	// Disk coordinate system
	vector3<T> o, xd, yd, zd;
	// Distance to disk plane from shading point
	T d;

	// Surface area
	T S;
public:
	DiskProjection():
		type(DiskProjection::ELLIPSE),
		o(vector3<T>()), xd(vector3<T>()), yd(vector3<T>()), zd(vector3<T>()),
		d(T(0)), S(T(0))
	{}

	/*
	 * Construct a spherical ellipse by projecting a disk over the unit sphere
	 * 		p  : shading point
	 * 		c  : disk center
	 * 		nd : disk normal
	 * 		r  : disk radius
	 */
	DiskProjection(const vector3<T>& p, const vector3<T>& cd, const vector3<T>& nd, T r);

	~DiskProjection();
public:
	/*
	 * Analytic total area of the disk projection
	 */
	T S_a() const;

	/*
	 * A point over the unit sphere corresponding to the mapping from the [0..1]^2 square
	 * to the disk projected area
	 * 		u : x coordinate over the [0..1]^2 square
	 * 		v : y coordinate over the [0..1]^2 square
	 */
	vector3<T> fromUnitSquare(const T u, const T v) const;
private:
	//static constexpr T eps  = std::numeric_limits<T>::epsilon();
	static constexpr T eps  = T(1e-6);
	static constexpr T pi   = T(3.14159265358979323846);
	static constexpr T pi_h = T(3.14159265358979323846)/T(2);
	static constexpr T pi_2 = T(2)*T(3.14159265358979323846);

	/*
	 * Calculate total projection area by elliptic integral
	 */
	T S_PI_t() const;

	/*
	 * Calculate fractional projection area by elliptic integral
	 *		u : fraction of total area
	 */
	T S_PI(const T u) const;
};

template<class T>
DiskProjection<T>::DiskProjection(const vector3<T>& o, const vector3<T>& cd, const vector3<T>& nd, T r)
{
	this->S = 1.0f;
	// Radius nearly zero
	if (r < eps) {
		this->S = 0.0f;
		return;
	}

	vector3<T> l = cd - o;

	// Distance from shading point to disk plane
	T d = -v3l::dot(l, nd);

	// Check if the disk normal is pointing to the opposite side
	// and calculate disk's coordinates Z axis
	if (d < 0.0f) {
		this->zd = nd;
		d = -d;
	} else {
		this->zd = -nd;
	}
	this->d = d;

	// Point inside disk plane
	if (this->d < eps) {
		this->S = 0;
		return;
	}

	// Build reference system
	this->xd = cross(l, this->zd);
	if (length2(this->xd) < eps) {
		T x = this->zd.x, y = this->zd.y, z = this->zd.z;

		if (std::abs(x) > std::abs(y)) {
			T invLength = 1.0f / std::sqrt(x*x + z*z);
			this->xd = vector3<T>(-z*invLength, 0.0f, x*invLength);
		} else {
			T invLength = 1.0f / std::sqrt(y*y + z*z);
			this->xd = vector3<T>(0.0f, z*invLength, -y*invLength);
		}
	} else {
		this->xd = normalize(this->xd);
	}
	this->yd = cross(this->zd, this->xd);

	// Store local coordinate system origin
	this->o = o;

	// Shading point projected in disk plane
	vector3<T> o_p = o + d*(this->zd);

	// Distance from projection to disk center
	T h = v3l::length(cd - o_p);

	// Minor axis y coordinates in the disk
	T y0 = h - r;
	T y1 = h + r;

	T inv_sqrt0 = 1.0f / std::sqrt(d*d + y0*y0);
	T inv_sqrt1 = 1.0f / std::sqrt(d*d + y1*y1);

	// Minor axis y coordinates projected in the sphere
	T y0_p = y0 * inv_sqrt0;
	T y1_p = y1 * inv_sqrt1;

	// Minor axis z coordinates projected in the sphere
	T d0_p = d * inv_sqrt0;
	T d1_p = d * inv_sqrt1;

	// Spherical ellipse center coordinates
	T yh_p = (y0_p + y1_p) / 2.0f;
	T dh_p = (d0_p + d1_p) / 2.0f;

	// SE center projected in the disk
	T yh = (d * yh_p) / dh_p;

	// Disk cord length at SE center projection
	T r_p = std::sqrt(r*r - std::pow(h - yh, 2));

	// Mayor axis x coordinates
//	T x0 = -r_p;
	T x1 = r_p;

	// Mayor axis x coordinates projected in the sphere
//	T x0_p = x0 / std::sqrt(x0*x0 + yh*yh + d*d);
	T x1_p = x1 / std::sqrt(x1*x1 + yh*yh + d*d);

	// Spherical ellipse axis
	T a = x1_p;
	T b = std::sqrt(std::pow(y1_p - y0_p, 2) + std::pow(d1_p - d0_p, 2)) / T(2);

	// The projection is a spherical ellipse
	this->type = DiskProjection::ELLIPSE;

	// Compute spherical ellipse constants
	SphericalEllipse<T>& ellip = this->projection.ellipse;

	ellip.a = a;
	ellip.b = b;

	{
		// SA in double precision to avoid numerical instabilities
		// TODO: more stable calculations using atan2()
		double a = double(ellip.a);
		double b = double(ellip.b);
		double c = (b*(1.0 - a*a)) / (a * std::sqrt(1.0f - b*b));
		double k = std::sqrt((a*a - b*b) / (1.0f - b*b));
		double n = (a*a - b*b) / (a*a*(1.0f - b*b));

		ellip.c = T(c);
		ellip.n = T(n);
		ellip.k = T(k);

		// Spherical ellipse area by elliptical integrals
//		const double pi_half = (3.14159265358979323846/2.0);
//		this->S = T(pi_half - c * ellint_3(k, n, pi_half));
	}

	// Precalc values
	ellip.mult = a * std::sqrt(1.0f - b*b) * b * std::sqrt(1.0f - a*a);
	ellip.phi_coeff = (b*std::sqrt(1.0f - a*a)) / (a*std::sqrt(1.0f - b*b));
	ellip.t_coeff = 1.0f / ellip.phi_coeff;

	// Reprojection constants
	T axis_norm = 1.0f / std::sqrt(yh_p*yh_p + dh_p*dh_p);
	ellip.yhp_norm = yh_p * axis_norm;
	ellip.dhp_norm = dh_p * axis_norm;
}

template<class T>
DiskProjection<T>::~DiskProjection()
{

}

template<class T> inline
T DiskProjection<T>::S_a() const
{
	if (this->S == 0.0f) {
		return 0.0f;
	}
	const SphericalEllipse<T>& ellip = this->projection.ellipse;

	double a = double(ellip.a);
	double b = double(ellip.b);
	double c = (b*(1.0f - a*a)) / (a * std::sqrt(1.0f - b*b));
	double k = std::sqrt((a*a - b*b) / (1.0f - b*b));
	double n = (a*a - b*b) / (a*a*(1.0f - b*b));

	// Spherical ellipse area by elliptical integrals
	constexpr double pi_half = (3.14159265358979323846f / 2.0f);
	return 4.0f * T(pi_half - c * ellint_3(k, n, pi_half));
}

template<class T> inline
T DiskProjection<T>::S_PI_t() const
{
	const SphericalEllipse<T>& ellip = this->projection.ellipse;

	T integr = ellint_3(ellip.k, ellip.n, pi_h);
	return pi_h - (ellip.c * integr);
}

template<class T> inline
T DiskProjection<T>::S_PI(const T u) const
{
	const SphericalEllipse<T>& ellip = this->projection.ellipse;

	T phi_u = pi_h * u;
	T t_u = std::atan(ellip.t_coeff * std::tan(phi_u));
	if (t_u == 0.0f) t_u = pi_h; // Avoid numerical inestabilities

	T integr = ellint_3(ellip.k, ellip.n, t_u);
	return phi_u - ellip.c * integr;
}

template<class T> inline
vector3<T> DiskProjection<T>::fromUnitSquare(const T u, const T v) const
{
	// Coords over the disk
	T x, y, z;

	if (this->type == DiskProjection::CIRCLE) {
		// Sampling a spherical circle
		const SphericalCircle<T>& circle = this->projection.circle;

		T phi = pi_2 * u;
		T r   = std::sqrt(v * circle.a);

		T cosphi = std::cos(phi);
		T sinphi = std::sqrt(1.0f - cosphi*cosphi);

		T rsqrt = r * std::sqrt(2.0f - r*r);

		T xu = rsqrt * cosphi;
		T yu = rsqrt * sinphi;
		T zu = 1.0f - r*r;

		// Reproject into disk plane from sphere surface
		T prop = this->d / zu;

		x = xu * prop;
		y = yu * prop;
		z = this->d;
	} else if (this->type == DiskProjection::WEDGE) {
		// Sampling a spherical wedge
		const SphericalWedge<T>& wedge = this->projection.wedge;

		T e1 = 2.0f * u - 1.0f;
		T e2 = 2.0f * v - 1.0f;

		T phi   = wedge.beta * std::abs(e1);

		T cosphi = std::cos(phi);
		T sinphi = std::sqrt(1.0f - cosphi*cosphi);

		T sqrte2 = std::sqrt(1.0f - e2*e2);

		// Local coords
		T xu = sign(e1) * sqrte2 * sinphi;
		T yu = e2;
		T zu = sqrte2 * cosphi;

		// Coords over the sphere
		T yu_p =  yu * wedge.dhp_norm + zu * wedge.yhp_norm;
		T zu_p = -yu * wedge.yhp_norm + zu * wedge.dhp_norm;

		T prop = this->d / zu_p;

		// Coords projected on the disk
		x = xu   * prop;
		y = yu_p * prop;
		z = this->d;
	} else {
		// Sampling a spherical ellipse
		const SphericalEllipse<T>& ellip = this->projection.ellipse;

		T e1 = 2.0f * u - 1.0f;
		T e2 = 2.0f * v - 1.0f;

		T e1minus = 1.0f - std::abs(e1);

		// Initial guess, assume planar ellipse
		// TODO: Better initial approximation
		T phi = std::atan(ellip.phi_coeff * std::tan(e1minus * pi_h));
		// Avoid going outside of function domain
		phi = clamp(phi, 0.0f, e1minus * pi_h);

		// Good enough precision
		static constexpr T tol = 1e-5f;

		// It will be necessary after the Newton loop
		T sinphi = std::sin(phi);

		// Newton root finding
		for (size_t i = 0; i < 100; ++i) {
			T sinphi2 = sinphi*sinphi;
			T cosphi2 = 1.0f - sinphi2;

			// Parametric angle
			T t = std::atan(ellip.t_coeff * std::tan(phi));

			T sint = std::sin(t);
			T sint2 = sint*sint;

			// Surface area for phi
			T omega = phi - ellip.c * ellint_3(ellip.k, ellip.n, t);

			T val = omega - (e1minus * this->S);

			// Derivative at phi
			T der = 1.0f - ellip.c *
				(ellip.mult) / ((1.0f - ellip.n * sint2) * std::sqrt(1.0f - ellip.k*ellip.k * sint2) *
				((ellip.mult) * cosphi2 + (ellip.mult) * sinphi2));

			// Avoid nearly zero divisions
			if (std::abs(der) < eps) {
				break;
			}

			// Avoid going outside of function domain
			phi = clamp(phi - val / der, 0.0f, e1minus * pi_h);

			// Calculate for the next loop
			sinphi = std::sin(phi);

			// Check if required tolerance has been reached
			if (std::abs(val) < tol) {
				break;
			}
		}

		const T& a = ellip.a;
		const T& b = ellip.b;

		// Calculate spherical radius and height for the point over the spherical ellipse edge
		// delimited by angle phi
		T cosphi = std::sqrt(1.0f - sinphi*sinphi);

		T rphi = (a*b) / std::sqrt(a*a * sinphi*sinphi + b*b * cosphi*cosphi);
		T hphi = std::sqrt(1.0f - rphi*rphi);

		// Sample uniformly a point along the ellipse arc delimited by phi
		T hu = 1.0f - (1.0f - hphi) * std::abs(e2);
		T ru = std::sqrt(1.0f - hu*hu);

		T xu = sign(e1) * ru * cosphi;
		T yu = sign(e2) * ru * sinphi;
		T zu = hu;

		// Coords over the sphere
		T yu_p =  yu * ellip.dhp_norm + zu * ellip.yhp_norm;
		T zu_p = -yu * ellip.yhp_norm + zu * ellip.dhp_norm;

		T prop = this->d / zu_p;

		// Coords projected on the disk
		x = xu   * prop;
		y = yu_p * prop;
		z = this->d;
	}

	// Point over the disk
	return this->o + x*this->xd + y*this->yd + z*this->zd;
}

template<class T> inline
std::ostream & operator<<(std::ostream & os, const DiskProjection<T>& p)
{
	os << std::fixed << std::setprecision(10);
	if (p.type == DiskProjection<T>::CIRCLE) {
	   os << "CIRCLE" << "\n";
	   os << "a = " << p.projection.circle.a << "\n";
	} else if (p.type == DiskProjection<T>::WEDGE) {
	   os << "WEDGE" << "\n";
	   os << "beta = " << p.projection.wedge.beta << "\n";
	} else {
	   os << "ELLIPSE" << "\n";
	   os << "a = " << p.projection.ellipse.a << "\n";
	   os << "b = " << p.projection.ellipse.b << "\n";
	   os << "c = " << p.projection.ellipse.c << "\n";
	   os << "n = " << p.projection.ellipse.n << "\n";
	   os << "k = " << p.projection.ellipse.k << "\n";
	   os << "phi_coeff = " << p.projection.ellipse.phi_coeff << "\n";
	   os << "t_coeff = " << p.projection.ellipse.t_coeff << "\n";
	   os << "dhp_norm = " << p.projection.ellipse.dhp_norm << "\n";
	   os << "yhp_norm = " << p.projection.ellipse.yhp_norm << "\n";
	}
	os << "o = "  << p.o << "\n";
	os << "xd = " << p.xd << "\n";
	os << "yd = " << p.yd << "\n";
	os << "zd = " << p.zd << "\n";
	os << "d = "  << p.d << "\n";
	os << "S = "  << p.S << std::endl;

	return os ;
}

#endif /* SRC_DISKPROJECTION_HPP_ */
