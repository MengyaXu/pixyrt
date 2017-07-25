//-----------------------------------------------------------------------------
/** @file vec2.inl
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
//-------------------------------------------------------------------------
inline Vector2::Vector2()
{
}

inline Vector2::Vector2(const Vector2& rhs)
	: x(rhs.x), y(rhs.y)
{

}

inline Vector2::Vector2(float value)
	: x(value), y(value)
{
}

inline Vector2::Vector2(float inX, float inY)
	: x(inX), y(inY)
{
}
//-----------------------------------------------------------------------------
inline Vector2 Vector2::operator+(const Vector2& rhs) const
{
	return Vector2(this->x + rhs.x, this->y + rhs.y);
}
//-------------------------------------------------------------------------
inline Vector2 Vector2::operator-(const Vector2& rhs) const
{
	return Vector2(this->x - rhs.x, this->y - rhs.y);
}
//-------------------------------------------------------------------------
inline Vector2 Vector2::operator*(float s) const
{
	return Vector2(this->x*s, this->y*s);
}
//-------------------------------------------------------------------------
inline Vector2 Vector2::operator/(float rhs) const
{
	return Vector2(this->x / rhs, this->y / rhs);
}

inline Vector2& Vector2::operator+=(const Vector2& rhs)
{
	x += rhs.x;
	y += rhs.y;
	return (*this);
}

inline Vector2& Vector2::operator-=(const Vector2& rhs)
{
	x -= rhs.x;
	y -= rhs.y;
	return (*this);
}

inline Vector2& Vector2::operator*=(float s)
{
	x *= s;
	y *= s;
	return (*this);
}

inline Vector2& Vector2::operator/=(float s)
{
	if (s != 0)
	{
		x /= s;
		y /= s;
	}
	return (*this);
}
//-------------------------------------------------------------------------
inline float cross(const Vector2& v1, const Vector2& v2)
{
	return v1.x*v2.y - v1.y*v2.x;
}
//-------------------------------------------------------------------------
inline float dot(const Vector2& v1, const Vector2& v2)
{
	return v1.x*v2.x + v1.y*v2.y;
}
//-------------------------------------------------------------------------
inline Vector2 normalize(const Vector2& v)
{
	float q = length(v);
	if (q > 0.0f) {
		return v / q;
	}
	return v;
}
//-------------------------------------------------------------------------
inline float length(const Vector2& v)
{
	return sqrtf(lengthSqr(v));
}
//-------------------------------------------------------------------------
inline float lengthSqr(const Vector2& v)
{
	return v.x*v.x + v.y*v.y;
}
//-------------------------------------------------------------------------
inline float determinant(const Vector2& v1, const Vector2& v2)
{
	return v1.getX() * v2.getY() - v2.getX() * v1.getY();
}
//-------------------------------------------------------------------------
inline Vector2 absPerElem(const Vector2& v)
{
	return Vector2(::fabsf(v.x), ::fabsf(v.y));
}
//-------------------------------------------------------------------------
inline float maxElem(const Vector2& v)
{
	return (v.x > v.y) ? v.x : v.y;
}
//-------------------------------------------------------------------------
inline float minElem(const Vector2& v)
{
	return (v.x < v.y) ? v.x : v.y;
}
//-------------------------------------------------------------------------
inline Vector2 maxPerElem(const Vector2& v1, const Vector2& v2)
{
	return Vector2(
		(v1.x > v2.x) ? v1.x : v2.x,
		(v1.y > v2.y) ? v1.y : v2.y);
}
//-------------------------------------------------------------------------
inline Vector2 minPerElem(const Vector2& v1, const Vector2& v2)
{
	return Vector2(
		(v1.x < v2.x) ? v1.x : v2.x,
		(v1.y < v2.y) ? v1.y : v2.y);
}
//-------------------------------------------------------------------------
inline Vector2 mulPerElem(const Vector2& v1, const Vector2& v2)
{
	return Vector2(v1.x * v2.x, v1.y * v2.y);
}
//-------------------------------------------------------------------------
inline Vector2 divPerElem(const Vector2& v1, const Vector2& v2)
{
	//	return Vector2(
	//		v2.x != 0 ? v1.x / v2.x : v1.x,
	//		v2.y != 0 ? v1.y / v2.y : v1.y);
	return Vector2(v1.x / v2.x, v1.y / v2.y);
}
//-------------------------------------------------------------------------
inline Vector2 recipPerElem(const Vector2& v)
{
	return Vector2(1.0f / v.getX(), 1.0f / v.getY());
}
//-------------------------------------------------------------------------
inline Vector2 sqrtPerElem(const Vector2& v)
{
	return Vector2(::sqrtf(v.getX()), ::sqrtf(v.getY()));
}
//-------------------------------------------------------------------------
inline Vector2 rsqrtPerElem(const Vector2& v)
{
	return Vector2(1.0f / ::sqrtf(v.getX()), 1.0f / ::sqrtf(v.getY()));
}
//-------------------------------------------------------------------------
inline Vector2 clampPerElem(const Vector2& v1, const Vector2& clampMin, const Vector2& clampMax)
{
	return maxPerElem(minPerElem(v1, clampMin), clampMax);
}
//-------------------------------------------------------------------------
inline float angle(const Vector2& v)
{
	return angle(Vector2::xAxis(), v);
	//	return ::atanf(v.getY(), v.getX());
}
//-------------------------------------------------------------------------
inline float angle(const Vector2& v1, const Vector2& v2)
{
	// stabler than acos
	return fabsf(atan2f(cross(v1, v2), dot(v1, v2)));
	//	float dotProduct = dot(v1, v2);
	//	float len1 = length(v1);
	//	float len2 = length(v2);
	//	float cosAngle = dotProduct / ::sqrtf(len1*len2);
	//	cosAngle = -1.0f > cosAngle ? -1.0f : cosAngle;
	//	cosAngle = +1.0f < cosAngle ? +1.0f : cosAngle;
	//	return ::acosf(cosAngle);
}
//-------------------------------------------------------------------------
inline Vector2 perp(const Vector2& v)
{
	return Vector2(-v.getY(), v.getX());
}
//-------------------------------------------------------------------------
/** ”½ŽžŒvŽü‚è‚É‰ñ“]
*/
inline Vector2 rotate(const Vector2& v, float angle)
{
	float c = cosf(angle);
	float s = sinf(angle);
	return Vector2(v.x*c - v.y*s, v.x*s + v.y*c);
}
//-------------------------------------------------------------------------
inline float sum(const Vector2& v)
{
	return v.x + v.y;
}
//-------------------------------------------------------------------------
inline Vector2 lerp(float t, const Vector2& v1, const Vector2& v2)
{
	return v1 + (v2 - v1) * t;
}
//-------------------------------------------------------------------------
inline Vector2 slerp(float t, const Vector2& v1, const Vector2& v2)
{
	float scale1, scale2;
	float cosAngle = dot(v1, v2);
	if (cosAngle < MAGE_DEFAULT_SLERP_TOL) {
		float angle = acosf(cosAngle);
		float recipSinAngle = 1.0f / sinf(angle);
		scale1 = sinf(((1.0f - t) * angle)) * recipSinAngle;
		scale2 = sinf((t*angle)) * recipSinAngle;
	}
	else {
		scale1 = 1.0f - t;
		scale2 = t;
	}

	return ((v1 * scale1) + (v2 * scale2));
}
//-------------------------------------------------------------------------
inline const Vector2& Vector2::xAxis()
{
	static Vector2 v(1.0f, 0.0f);
	return v;
}
//-------------------------------------------------------------------------
inline const Vector2& Vector2::yAxis()
{
	static Vector2 v(0.0f, 1.0f);
	return v;
}
//-------------------------------------------------------------------------
inline const Vector2& Vector2::zero()
{
	static Vector2 v(0.0f);
	return v;
}
//-------------------------------------------------------------------------
inline Vector2 operator*(float lhs, const Vector2& rhs)
{
	return rhs * lhs;
}
//-------------------------------------------------------------------------
inline void loadXY(Vector2& vec, const float* fptr)
{
	vec.setX(fptr[0]);
	vec.setY(fptr[1]);
}
//-------------------------------------------------------------------------
inline void storeXY(const Vector2& vec, float* fptr)
{
	fptr[0] = vec.getX();
	fptr[1] = vec.getY();
}
//-------------------------------------------------------------------------
inline Vector2 select(const Vector2& v1, const Vector2& v2, bool select2)
{
	return select2 ? v2 : v1;
}
