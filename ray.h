#pragma once
#include "vec.h"

class ray {
public:
	ray() {}
	ray(const vec3& origin, const vec3& direction, float ti = 0.0) : m_origin(origin), m_direction(direction), m_time(ti) {}
	const vec3& origin() const { return m_origin; }
	const vec3& direction() const { return m_direction; }
	float time() const { return m_time; }
	vec3 point_at_parameter(float t) const { return m_origin + m_direction * t; }

private:
	vec3 m_origin;
	vec3 m_direction;
	float m_time;
};