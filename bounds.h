#pragma once

#include "ray.h"

class aabb {
public:
	aabb() {}
	aabb(const vec3& a, const vec3& b) : m_min(a), m_max(b) {}

	const vec3& lower() const { return m_min; }
	const vec3& upper() const { return m_max; }

	bool hit(const ray& r, float tmin, float tmax) const;

	vec3 m_min;
	vec3 m_max;
};

//bool aabb::hit(const ray& r, float tmin, float tmax) const
//{
//	for (int a = 0; a<3; ++a)
//	{
//		float t0 = std::min(
//			(m_min[a] - r.origin()[a]) / r.direction()[a],
//			(m_max[a] - r.origin()[a]) / r.direction()[a]);
//		float t1 = std::max(
//			(m_min[a] - r.origin()[a]) / r.direction()[a],
//			(m_max[a] - r.origin()[a]) / r.direction()[a]);
//		tmin = std::min(tmin, t0);
//		tmax = std::max(tmax, t1);
//		if (tmax <= tmin) {
//			return false;
//		}
//	}
//
//	return true;
//}

bool aabb::hit(const ray& r, float tmin, float tmax) const
{
	for (int a=0; a<3; ++a)
	{
		float invD = 1.0f / r.direction()[a];
		float t0 = (lower()[a] - r.origin()[a]) * invD;
		float t1 = (upper()[a] - r.origin()[a]) * invD;
		if (invD < 0.0f) {
			std::swap(t0,t1);
		}
		tmin = t0 > tmin ? t0 : tmin;
		tmax = t1 < tmax ? t1 : tmax;
		if (tmax <= tmin) {
			return false;
		}
	}

	return true;
}

inline aabb surrounding_box(const aabb& a, const aabb& b)
{
	return aabb(minPerElem(a.lower(), b.lower()), maxPerElem(a.upper(), b.upper()));
}
