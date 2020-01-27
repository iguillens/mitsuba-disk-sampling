/*
 * mathUtil.hpp
 *
 * Single header shitty miscellaneous math library
 *      Author:  Ibon Guillen
 *      Version: 2017-03-21
 */

#ifndef __MATHUTIL__
#define __MATHUTIL__

#include <algorithm>
#include <cmath>

namespace sml {
	template<class T> inline constexpr
	T sign(const T& x)
	{
		return std::copysign(T(1), x);
	}

	template<class T> inline constexpr
	T clamp(const T& x, const T& min, const T& max)
	{
		return std::min(max, std::max(min, x));
	}
}

#endif // __MATHUTIL__
