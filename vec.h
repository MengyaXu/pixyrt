#pragma once

#define _USE_MATH_DEFINES
#include <math.h> // M_PI, etc
#include <float.h> // FLT_MIN, FLT_MAX
#include <algorithm>
#include <iostream>

#define PI 3.14159265359f
#define PI2 6.28318530718f
#define RECIPROCAL_PI 0.31830988618f
#define RECIPROCAL_PI2 0.15915494f
#define LOG2 1.442695f
#define EPSILON 1e-6
#define GAMMA_FACTOR 2.2f

// https://github.com/kikikikina/drand48/blob/master/main.cpp
#include <random>
inline float drand48() {
	return float(((double)(rand()) / (RAND_MAX))); /* RAND_MAX = 32767 */
}
//inline double drand48() {
//	return ((double)(rand()) / (RAND_MAX)); /* RAND_MAX = 32767 */
//}

inline float pow2(float x) { return x*x; }
inline float pow3(float x) { return x*x*x; }
inline float pow4(float x) { return x*x*x*x; }
inline float pow5(float x) { return x*x*x*x*x; }
inline float clamp(float x, float a, float b) { return x < a ? a : x > b ? b : x; }
inline float saturate(float x) { return x < 0.0 ? 0.0 : x > 1.0 ? 1.0 : x; }
inline float recip(float x) { return 1.0/x; }
inline float mix(float a, float b, float t) { return a*(1.0f-t) + b*t; /* return a + (b-a) * t; */ }
inline float step(float edge, float x) { return (x < edge) ? 0.0f : 1.0f; }
inline float smoothstep(float a, float b, float t) { if (a>=b) return 0.0f; float x = saturate((t-a) / (b-a)); return x*x*(3-2*t); }
inline float radians(float deg) { return (deg/180)*M_PI; }
inline float degrees(float rad) { return (rad/M_PI)*180; }

#include <vectormath/scalar/cpp/vectormath_aos.h>
using namespace Vectormath::Aos;
#include "vec2.h"
#include "vec2.inl"
typedef Vector2 vec2;
typedef Vector3 vec3;
typedef Vector4 vec4;
typedef Vector3 col3;

inline Vector3 saturate(const Vector3& v)
{
	return minPerElem(maxPerElem(v, Vector3(0.0f)), Vector3(1.0f));
}

inline vec3 rndvec() 
{
	return vec3(drand48(), drand48(), drand48());
}

inline vec3 reflect(const vec3& v, const vec3& n)
{
	return v - 2*dot(v,n)*n; 
}

inline bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted)
{
	vec3 uv = normalize(v);
	float dt = dot(uv,n);
	float D = 1.0 - pow2(ni_over_nt) * (1- pow2(dt));
	if (D > 0) {
		refracted = ni_over_nt * (uv - n*dt) - n*sqrt(D);
		return true;
	}
	else {
		return false;
	}
}

inline vec3 linear_to_gamma(const vec3& v, float gammaFactor)
{
	float recipGammaFactor = recip(gammaFactor);
	return vec3(
		powf(v.getX(), recipGammaFactor),
		powf(v.getY(), recipGammaFactor),
		powf(v.getZ(), recipGammaFactor));
}

inline vec3 gamma_to_linear(const vec3& v, float gammaFactor)
{
	return vec3(
		powf(v.getX(), gammaFactor),
		powf(v.getY(), gammaFactor),
		powf(v.getZ(), gammaFactor));
}

inline vec3 de_nan(const vec3& c)
{
	return vec3(
		!(c[0] == c[0]) ? 0 : c[0],
		!(c[1] == c[1]) ? 0 : c[1],
		!(c[2] == c[2]) ? 0 : c[2]);
}

inline float schlick(float cosine, float roi)
{
	float r0 = pow2((1-roi) / (1+roi));
	return r0 + (1-r0) * pow5(1-cosine);
}

inline float schlick(float cosine, float f0, float f90)
{
	return f0 + (f90 - f0) * pow5(1 - cosine);
}

void get_sphere_uv(const vec3& p, float& u, float& v) {
	float phi = atan2(p.getZ(), p.getX());
	float theta = asin(p.getY());
	u = 1 - (phi + M_PI) / (2 * M_PI);
	v = (theta + M_PI / 2) / M_PI;
}