/*
 * SphericalEllipseSampler.hpp
 *
 *      Author: ibon
 */

#ifndef __SPHERICALELLIPSESAMPLING_HPP__
#define __SPHERICALELLIPSESAMPLING_HPP__

#include <algorithm>
#include <numeric>
#include <iostream>

#define BOOST_MATH_MAX_SERIES_ITERATION_POLICY 500
#include <boost/math/special_functions/ellint_3.hpp>
namespace {
	using boost::math::ellint_3;
}
#include <boost/math/special_functions/ellint_1.hpp>
namespace {
	using boost::math::ellint_1;
}

#include "vector3.hpp"
namespace {
	using v3l::vector3;
}

#include "mathUtil.hpp"

// Select appropiate precomputed table size
//#include "CDFTable32x32.hpp"
//#include "CDFTable64x64.hpp"
//#include "CDFTable128x128.hpp"
//#include "CDFTable256x256.hpp"
#include "CDFTable512x512.hpp"
//#include "CDFTable1024x1024.hpp"


template<class T>
class SphericalEllipseSampler
{
public:
	enum class Type {
		ZERO, // Indicates unsampleable case
		POLAR,
		POLARTAB,
		CONCENTRIC,
		CONCENTRICTAB,
		CYLINDRICAL,
		CYLINDRICALTAB
	} type;

	// Parameters of each algorithm variant
	union {
		struct BoothData {
			// Elliptical integral constants
			T c, n, k;
			// Precalculated constants
			T phi_coeff, t_coeff, asqrtb, bsqrta;
			// Quadrant area
			T S;
		} boothData;
		struct {
			// Tabulated CDF
			double* CDF1D;
			// Quadrant area
			T S;
		} boothTabData;
		struct UrenaData {
			// Spherical ellipse arcs
			T alpha, beta;
			// Tangent ellipse params
			T a_p, b_p;
			// Elliptical integral constants
			T n, k, m;
			// Precalculated constants
			T p, c1;
			// Half area
			T S_2;
		} urenaData;
	} data;

	// Spherical ellipse params
	T a, b;
private:
	// Constants
	static constexpr T pi     = T(3.14159265358979323846);
	static constexpr T pi_2   = T(3.14159265358979323846/2.0);
	static constexpr T pi_4   = T(3.14159265358979323846/4.0);
	static constexpr T two_pi = T(2.0*3.14159265358979323846);
	static constexpr T eps    = T(1e-6);
	static constexpr T tol    = T(1e-5);
public:
	vector3<T> sample(const T e1, const T e2, T& pdf) const;

	T samplePdf(const vector3<T>& p) const;
private:
	/*
	 * Azimuth angle that covers the specified fractional area over the spherical ellipse,
	 * as found by inverting Booth's expression by Newton iterations
	 * 		u   : fraction of total area
	 * 		pdf : probability of the generated sample
	 * 		tol : tolerance
	 */
	T BoothInversionNewton(const T u, T &pdf, const T tol = tol) const;
	/*
	 * Azimuth angle that covers the specified fractional area over the spherical ellipse,
	 * as found by inverting a tabulated approximation by spherical triangles
	 * 		u   : fraction of total area
	 * 		pdf : probability of the generated sample
	 * 		h_i : height bound of the sampled spherical triangle
	 */
	T BoothInversionTable(const T u, T& pdf, T& h_i) const;
	/*
	 * Probability of sampling a given azimuth angle using table inversion
	 * 		phi : sampled azimuth angle
	 */
	T BoothPdfTable(const T phi) const;
	/*
	 * Fractional minor arc that covers the specified fractional area over the spherical ellipse,
	 * as found by inverting Urena's expression by Newton iterations
	 * 		u   : fraction of total area
	 * 		tol : tolerance
	 * 		pdf : probability of the generated sample
	 */
	T UrenaInversionNewton(T u, T &pdf, T tol = tol) const;
};

template<class T> inline
vector3<T> SphericalEllipseSampler<T>::sample(const T e1, const T e2, T& pdf) const
{
	vector3<T> p(T(0));

	switch(this->type) {
		case Type::POLAR :
		case Type::POLARTAB :
		{
			// This is a continuous variation of polar mapping that has a single discontinuity
			// and avoids the disk-to-square-and-back step. The branching could probably be
			// greatly reduced

			// Choose sampling quadrant
			T u;
			if (e1 >= 0.75f) {
				// Fourth quadrant
				u = 1.0f - 4.0f * (e1 - 0.75f);
			} else if (e1 >= 0.5f) {
				// Third quadrant
				u = 4.0f * (e1 - 0.5f);
			} else if (e1 >= 0.25f) {
				// Second quadrant
				u = 1.0f - 4.0f * (e1 - 0.25f);
			} else {
				// First quadrant
				u = 4.0f * e1;
			}

			T v = e2;

			// Sample azimuth angle
			T phi_u, h_i = 0.0f;
			if (this->type == Type::POLAR) {
				// Newton root finding
				phi_u = this->BoothInversionNewton(u, pdf);
			} else {
				// Tabulated CDF inversion
				phi_u = this->BoothInversionTable(u, pdf, h_i);
			}

			// Translate to sampling quadrant
			if (e1 >= 0.75f) {
				// Fourth quadrant
				phi_u = two_pi - phi_u;
			} else if (e1 >= 0.5f) {
				// Third quadrant
				phi_u = pi + phi_u;
			} else if (e1 >= 0.25f) {
				// Second quadrant
				phi_u = pi - phi_u;
			}

			const T a = this->a;
			const T b = this->b;

			// Calculate spherical radius and height for the point over the
			// spherical ellipse edge defined by phi
			T sinphi = std::sin(phi_u);
			T cosphi = std::cos(phi_u);

			T r_u = (a*b)/std::sqrt(a*a * sinphi*sinphi + b*b * cosphi*cosphi);
			T h_u = std::sqrt(1.0f - r_u*r_u);

			T h_vu;
			if (this->type == Type::POLAR) {
				// Sample uniformly a point along the spherical ellipse arc
				h_vu = 1.0f - (1.0f - h_u) * v;
			} else {
				// Sample uniformly a point along the spherical triangle arc
				h_vu = 1.0f - (1.0f - h_i) * v;

				// Check if the sampled point is inside the spherical ellipse
				// or reject otherwise
				if (h_vu < h_u) {
					pdf = 0.0f;
					return vector3<T>(0.0f);
				}
			}

			T r_vu = std::sqrt(1.0f - h_vu*h_vu);

			// Local coordinates
			T x_vu = r_vu * cosphi;
			T y_vu = r_vu * sinphi;
			T z_vu = h_vu;

			p = vector3<T>(x_vu, y_vu, z_vu);
		}
		break;
		case Type::CONCENTRIC :
		case Type::CONCENTRICTAB :
		{
			// A variant of Shirley's clever continuous mapping that unfortunately contains tons
			// of branching because we can only sample the first quadrant of the spherical
			// ellipse, so we need to go back and forth between quadrants

			T u = 2.0f * e1 - 1.0f;
			T v = 2.0f * e2 - 1.0f;

			// Choose sampling quadrant
			T r, theta;
			if (u == 0.0f && v == 0.0f) {
				r = 0.0f;
				theta = 0.0f;
			} else if (std::abs(u) > std::abs(v)) {
				r = u;
				if (u < 0.0f) {
					// Third quadrant
					theta = pi + pi_4 * v/u;
				} else if (v >= 0.0f) {
					// First quadrant positive
					theta = pi_4 * v/u;
				} else {
					// First quadrant negative
					theta = two_pi + pi_4 * v/u;
				}
			} else {
				r = v;
				if (v >= 0.0f) {
					// Second quadrant
					theta = pi_2 - pi_4 * u/v;
				} else {
					// Fourth quadrant
					theta = 3.0f*pi_2 - pi_4 * u/v;
				}
			}

			// Translate sampling coords
			T u_p;
			if (theta > 3.0f*pi_2) {
				u_p = 1.0f - (theta - 3.0f*pi_2) / pi_2;
			} else if (theta > pi) {
				u_p = (theta - pi) / pi_2;
			} else if (theta > pi_2) {
				u_p = 1.0f - (theta - pi_2) / pi_2;
			} else {
				u_p = theta / pi_2;
			}

			T v_p = r*r;

			// Sample azimuth angle
			T phi, h_i = 0.0f;
			if (this->type == Type::CONCENTRIC) {
				// Newton root finding
				phi = this->BoothInversionNewton(u_p, pdf);
			} else {
				// Tabulated CDF inversion
				phi = this->BoothInversionTable(u_p, pdf, h_i);
			}

			// Translate sampled angle to sampling quadrant
			T phi_u;
			if (theta > 3.0f*pi_2) {
				phi_u = two_pi - phi;
			} else if (theta > pi) {
				phi_u = pi + phi;
			} else if (theta > pi_2) {
				phi_u = pi - phi;
			} else {
				phi_u = phi;
			}

			const T a = this->a;
			const T b = this->b;

			// Calculate spherical radius and height for the point over the spherical ellipse edge
			// delimited by angle phi
			T sinphi = std::sin(phi_u);
			T cosphi = std::cos(phi_u);

			T r_u = (a*b)/std::sqrt(a*a * sinphi*sinphi + b*b * cosphi*cosphi);
			T h_u = std::sqrt(1.0f - r_u*r_u);

			T h_vu;
			if (this->type == Type::CONCENTRIC) {
				// Sample uniformly a point along the spherical ellipse arc
				h_vu = 1.0f - (1.0f - h_u) * v_p;
			} else {
				// Sample uniformly a point along the spherical triangle arc
				h_vu = 1.0f - (1.0f - h_i) * v_p;

				// Check if the sampled point is inside the spherical ellipse
				// or reject otherwise
				if (h_vu < h_u) {
					pdf = -1.0f;
					return vector3<T>(0.0f);
				}
			}

			T r_vu = std::sqrt(1.0f - h_vu*h_vu);

			// Local coordinates
			T x_vu = r_vu * cosphi;
			T y_vu = r_vu * sinphi;
			T z_vu = h_vu;

			p = vector3<T>(x_vu, y_vu, z_vu);
		}
		break;
		case Type::CYLINDRICAL :
		case Type::CYLINDRICALTAB :
		{
			T beta_t = this->UrenaInversionNewton(e1, pdf);

			const T a_p = this->data.urenaData.a_p;
			const T b_p = this->data.urenaData.b_p;

			// Aboslute coords. of the two points that delimit the tangent
			// ellipse's chord defined by beta
			T y_1 = std::tan(beta_t);
			T x_1 = a_p * std::sqrt(1.0f - (y_1 * y_1) / (b_p * b_p));

			// Limiting points reprojected into the sphere
//			vector3<T>s0_p = normalize(vector3<T>(-x_1, y_1 , 1.0f));
			vector3<T>s1_p = normalize(vector3<T>( x_1, y_1 , 1.0f));

			// Linear interpolation between limiting points
			T w = 2.0f * e2 - 1.0f;
			T f = std::sqrt((1.0f - std::pow(w * s1_p.x, 2))/
					        (1.0f - std::pow(s1_p.x, 2)));

			p = vector3<T>(s1_p.x * w, s1_p.y * f, s1_p.z * f);
		}
		break;
		case Type::ZERO :
		{
			pdf = 0.0f;
			p = vector3<T>(0.0f);
		}
	}

	return p;
}

//TODO: Find pdf by Table Inversion
template<class T> inline
T SphericalEllipseSampler<T>::samplePdf(const vector3<T>&) const
{
	switch(this->type) {
		case Type::POLAR :
		case Type::CONCENTRIC :
		{
			auto& data = this->data.boothData;
			return 1.0f / (4.0f * data.S);
		}
		break;
		case Type::POLARTAB :
		case Type::CONCENTRICTAB :
		{
			auto& data = this->data.boothTabData;
			return 1.0f / (4.0f * data.S);
		}
		break;
		case Type::CYLINDRICAL :
		case Type::CYLINDRICALTAB :
		{
			auto& data = this->data.urenaData;
			return 1.0f / (2.0f * data.S_2);
		}
		break;
		case Type::ZERO :
		{
			return 0.0f;
		}
	}

	return T(0);
}

template<class T> inline
T SphericalEllipseSampler<T>::BoothInversionNewton(const T u, T &pdf, const T tol) const
{
	auto& data = this->data.boothData;

	// Avoid going outside of function domain
	T phi_max = u * T(pi_2);

	// Initial guess, assume planar ellipse
	// TODO: Better initial approximation
	T phi = std::atan(data.phi_coeff * std::tan(phi_max));
	phi = sml::clamp(phi, T(0.0f), phi_max);

	// Target fractional area
	T S_u = u * data.S;

	// Newton iterations
	for (size_t i = 0; i < 100; ++i) {
		T sinphi = std::sin(phi);

		T sinphi2 = sinphi*sinphi;
		T cosphi2 = 1.0f - sinphi2;

		// Parametric angle
		T t = sml::clamp(std::atan(data.t_coeff * std::tan(phi)), T(0), T(pi_2));

		T sint = std::sin(t);
		T sint2 = sint*sint;

		// Surface area for phi
		T omega = phi - data.c * ellint_3(data.k, data.n, t);

		T diff = omega - S_u;

		// Check if required tolerance has been reached
		if (std::abs(diff) < tol) {
			break;
		}

		// Derivative at phi
		T der = 1.0f - ((data.c * data.asqrtb * data.bsqrta) /
			((1.0f - data.n * sint2) * std::sqrt(1.0f - data.k*data.k * sint2) *
			((data.bsqrta) * cosphi2 + (data.asqrtb) * sinphi2)));

		// Avoid nearly zero divisions
		if (std::abs(der) < T(eps)) {
			break;
		}

		// Avoid going outside of function domain
		phi = sml::clamp(phi - diff/der, T(0.0f), phi_max);
	}

	pdf = 1.0f / (4.0f * data.S);

	return phi;
}

template<class T> inline
T SphericalEllipseSampler<T>::BoothInversionTable(const T u, T& pdf, T& h_i) const
{
	auto& data = this->data.boothTabData;

	// Binary search
	double* ptr = std::upper_bound(data.CDF1D, data.CDF1D + CDF_SIZE, u);
	size_t index = size_t(std::max(ptrdiff_t(0), ptr - data.CDF1D));

	T CDF_0 = (index > 0) ? T(data.CDF1D[index - 1]) : 0.0f;
	T CDF_1 = T(data.CDF1D[index]);
	T CDF_D = CDF_1 - CDF_0;

	// Calculate height at the start of the table entry
	T phi_i = (T(index) * T(pi_2)) / T(CDF_SIZE);
	T sinphi_i = std::sin(phi_i);
	T cosphi_i = std::cos(phi_i);

	const T a = this->a;
	const T b = this->b;

	T rphi_i = (a*b)/std::sqrt(a*a * sinphi_i*sinphi_i + b*b * cosphi_i*cosphi_i);
	h_i = std::sqrt(1.0f - rphi_i*rphi_i);

	// Find interpolated azimuth angle for sampling
	T du = (u - CDF_0) / CDF_D;
	T phi = ((T(index) + du) * T(pi_2)) / T(CDF_SIZE);

	pdf = 1.0f / (4.0f * data.S);

	return phi;
}

template<class T> inline
T SphericalEllipseSampler<T>::UrenaInversionNewton(T u, T &pdf, T tol) const
{
	auto& data = this->data.urenaData;

	// Initial guess, assume linear relation
	T t = u;

	// Target fractional area
	T S_u = 2.0f * data.S_2 * u;

	// Newton iterations
	for (size_t i = 0; i < 100; ++i) {
		T beta_t = std::abs(data.beta * (2.0f * t - 1.0f));

		T Omega_t;
		if (t == 0.0f) {
			Omega_t = 0.0f;
		} else if (t == 0.5f) {
			Omega_t = data.S_2;
		} else if (t == 1.0f) {
			Omega_t = 2.0f * data.S_2;
		} else {
			T ampl = -std::asin(std::min(T(1.0f), std::tan(beta_t) / data.b_p));
			T I = data.c1 * data.b_p * (data.p * ellint_1(data.k, ampl) -
			      (data.p + 1.0f) * ellint_3(data.k, data.n, ampl));

			// Partial area
			Omega_t = data.S_2 + sml::sign(t - 0.5f) * I;
		}

		T diff = Omega_t - S_u;

		// Check if required tolerance is reached
		if (std::abs(diff) < tol) {
			break;
		}

		T tan_beta_t = std::tan(beta_t);
		T delta = data.c1 * std::sqrt((1.0f - data.p * tan_beta_t * tan_beta_t)/
				  (1.0f - data.m * data.p * tan_beta_t * tan_beta_t));

		T Omega_t_d = 2.0f * data.beta * delta;

		// Avoid nearly zero divisions
		if (std::abs(Omega_t_d) < T(eps)) {
			break;
		}

		t = t - (Omega_t - S_u) / Omega_t_d;
	}

	pdf = 1.0f / (2.0f * data.S_2);

	return (2.0f * t - 1.0f) * data.beta;
}

#endif // __SPHERICALELLIPSESAMPLING_HPP__
