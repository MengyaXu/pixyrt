//-----------------------------------------------------------------------------
/** @file vec2.h
*/
//-----------------------------------------------------------------------------
/*
* Copyright (c) 2009 mebiusbox software. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
*/
#pragma once

#define MAGE_MATH_PI_OVER_2 1.570796327f
// tol : ベクトルが同一線上にあるかどうかを判定するための許容値
#define MAGE_DEFAULT_SLERP_TOL 0.999f

// vectormath compatible version
class Vector2
{
public:
	union
	{
		struct { float x, y; };
		float v[2];
	};

	inline Vector2();
	inline Vector2(float inX, float inY);
	inline Vector2(const Vector2& other);
	explicit inline Vector2(float value);
	//-----------------------------------------------------------------------------
	/// Set the x element of a 2D vector.
	inline Vector2& setX(float value) { x = value; return *this; }
	/// Set the y element of a 2D vector.
	inline Vector2& setY(float value) { y = value; return *this; }
	/// Get the x element of a 2D vector.
	inline float getX() const { return x; }
	/// Get the y element of a 2D Vector.
	inline float getY() const { return y; }
	/// Set an x or y element of a 2D vector by index.
	inline Vector2& setElem(int idx, float value) { v[idx] = value; return *this; }
	/// Get an x or y element of a 2D vector by index.
	inline float getElem(int idx) const { return v[idx]; }
	//-----------------------------------------------------------------------------
	inline float& operator[](int idx) { return v[idx]; }
	inline float  operator[](int idx) const { return v[idx]; }
	//-----------------------------------------------------------------------------
	inline Vector2 operator+(const Vector2& rhs) const;
	inline Vector2 operator-(const Vector2& rhs) const;
	inline Vector2 operator*(float s) const;
	inline Vector2 operator/(float s) const;
	inline Vector2& operator+=(const Vector2& rhs);
	inline Vector2& operator-=(const Vector2& rhs);
	inline Vector2& operator*=(float s);
	inline Vector2& operator/=(float s);
	inline Vector2  operator-() const { return Vector2(-x, -y); }
	//-----------------------------------------------------------------------------
	inline operator float* () { return v; }
	inline operator const float* () const { return v; }
	//-------------------------------------------------------------------------
	/// Construct x axis.
	static inline const Vector2& xAxis();
	/// COnstruct y axix.
	static inline const Vector2& yAxis();
	/// Construct zero.
	static inline const Vector2& zero();
};
//-------------------------------------------------------------------------
/// Compute cross product of two 2D vectors.
inline float cross(const Vector2& v1, const Vector2& v2);
/// Compute the dot product of two 2D vectors.
inline float dot(const Vector2& v1, const Vector2& v2);
/// Compute the length of a 2D vector.
inline float length(const Vector2& v);
/// Compute the square of the length of a 2D vector.
inline float lengthSqr(const Vector2& v);
/// Compute the determinant of the matrix created from two 2D vectors columns.
inline float determinant(const Vector2& v1, const Vector2& v2);
/// Compute the absolute value of a 2D vector per element.
inline Vector2 absPerElem(const Vector2& v);
/// Maximum element of a 2D vector.
inline float maxElem(const Vector2& v);
/// Minimum element of a 2D vector.
inline float minElem(const Vector2& v);
/// Maximum of two 2D vectors per element.
inline Vector2 maxPerElem(const Vector2& v1, const Vector2& v2);
/// Minimum of two 2D vectors per element.
inline Vector2 minPerElem(const Vector2& v1, const Vector2& v2);
/// Multiply of two 2D vectors per element.
inline Vector2 mulPerElem(const Vector2& v1, const Vector2& v2);
/// Divice two 2D vectors per element.
inline Vector2 divPerElem(const Vector2& v1, const Vector2& v2);
/// ２次元ベクトルの各要素の逆数を計算します.
/// Compute the reciprocal of a 2D vector per element.
inline Vector2 recipPerElem(const Vector2& v);
/// Compute the square root of a 2D vector per element.
inline Vector2 sqrtPerElem(const Vector2& v);
/// Compute the reciprocal square root fo a 2D vecdtor per element.
inline Vector2 rsqrtPerElem(const Vector2& v);
/// Clamp each element of a 2D vector between corresponding elements specifying
/// minimum and maximum values.
inline Vector2 clampPerElem(const Vector2& v1, const Vector2& clampMin, const Vector2& clampMax);
/// Compute the angle of a 2D vector against the x-axis
inline float angle(const Vector2& v);
/// Compute the angle between two 2D vectors.
inline float angle(const Vector2& v1, const Vector2& v2);
/// ２次元ベクトルに対して垂直な２次元ベクトルを計算します.
/// Compute a 2D vector perpendicular to the 2D vector.
inline Vector2 perp(const Vector2& v);
/// Rotate a 2D vector.
inline Vector2 rotate(const Vector2& v, float angle);
/// Compute the sum of all elements of a 2D vector.
inline float sum(const Vector2& v);
/// Normalize a 2D vector.
inline Vector2 normalize(const Vector2& v);
/// Linear interpolation between two 2D vectors.
inline Vector2 lerp(float t, const Vector2& v1, const Vector2& v2);
/// Spherical linear interpolation between two 2D vectors.
inline Vector2 slerp(float t, const Vector2& v1, const Vector2& v2);
inline Vector2 operator*(float lhs, const Vector2& rhs);
inline void loadXY(Vector2& vec, const float* fptr);
inline void storeXY(const Vector2& vec, float* fptr);
inline Vector2 select(const Vector2& v1, const Vector2& v2, bool select2);