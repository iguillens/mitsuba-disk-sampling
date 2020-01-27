/*
 * Vector3.hpp
 *
 * Single header shitty vector library
 *      Author:  Ibon Guillen
 *      Version: 2017-03-21
 */

#ifndef __VECTOR3__
#define __VECTOR3__

#include <cmath>
#include <cstdint>
#include <iostream>

namespace v3l {

template<class T>
class vector3
{
public:
	constexpr vector3(T n = T(0)): x(n), y(n), z(n) {}
	constexpr vector3(T x, T y, T z): x(x), y(y), z(z) {}

	constexpr vector3(vector3<T>& v): x(v.x), y(v.y), z(v.z) {}
	constexpr vector3(const vector3<T>& v): x(v.x), y(v.y), z(v.z) {}

	~vector3() {}
public:
	T x, y, z;
};

// Defines

using vector3d = vector3<double>;
using vector3f = vector3<float>;

// Scalar operations

template<class T> constexpr inline
vector3<T> operator*(const vector3<T>& v, const T& n)
{
	return vector3<T>(v.x * n, v.y * n, v.z * n);
}

template<class T> constexpr inline
vector3<T> operator*(const T& n, const vector3<T>& v)
{
	return vector3<T>(n * v.x, n * v.y, n * v.z);
}

template<class T> constexpr inline
vector3<T> operator/(const vector3<T>& v, const T& n)
{
	return vector3<T>(v.x / n, v.y / n, v.z / n);
}

template<class T> constexpr inline
vector3<T> operator/(const T& n, const vector3<T>& v)
{
	return vector3<T>(n / v.x, n / v.y, n / v.z);
}

// Element-wise operations

template<class T> constexpr inline
vector3<T> operator-(const vector3<T>& v)
{
	return vector3<T>(-v.x, -v.y, -v.z);
}

template<class T> constexpr inline
vector3<T> operator+(const vector3<T>& v1, const vector3<T>& v2)
{
	return vector3<T>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

template<class T> constexpr inline
vector3<T> operator-(const vector3<T>& v1, const vector3<T>& v2)
{
	return vector3<T>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

template<class T> constexpr inline
vector3<T> operator*(const vector3<T>& v1, const vector3<T>& v2)
{
	return vector3<T>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
}

template<class T> constexpr inline
vector3<T> operator/(const vector3<T>& v1, const vector3<T>& v2)
{
	return vector3<T>(v1.x / v2.x, v1.y / v2.y, v1.z / v2.z);
}

// Vector operations

template<class T> constexpr inline
T dot(const vector3<T>& v1, const vector3<T>& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template<class T> constexpr inline
T operator|(const vector3<T>& v1, const vector3<T>& v2)
{
	return dot(v1, v2);
}

template<class T> constexpr inline
vector3<T> cross(const vector3<T>& v1, const vector3<T>& v2)
{
	return vector3<T>(v1.y * v2.z - v1.z * v2.y,
					  v1.z * v2.x - v1.x * v2.z,
	                  v1.x * v2.y - v1.y * v2.x);
}

template<class T> constexpr inline
vector3<T> operator^(const vector3<T>& v1, const vector3<T>& v2)
{
	return cross(v1, v2);
}

template<class T> constexpr inline
T length(const vector3<T>& v)
{
	return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

template<class T> constexpr inline
T length2(const vector3<T>& v)
{
	return v.x * v.x + v.y * v.y + v.z * v.z;
}

template<class T> constexpr inline
vector3<T> normalize(const vector3<T>& v)
{
	T norm = T(1)/length(v);
	return vector3<T>(v.x * norm, v.y * norm, v.z * norm);
}

template<class T> inline
std::ostream & operator<<(std::ostream & os, const vector3<T>& v)
{
   os << "(" << v.x << "," << v.y << "," << v.z << ")";
   return os ;
}

} // namespace v3l

#endif // __VECTOR3__
