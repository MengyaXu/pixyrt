#pragma once
#include "vec.h"

inline vec3 random_cosine_direction()
{
	float r1 = drand48();
	float r2 = drand48();
	float z = sqrt(1 - r2 + FLT_EPSILON);
	float phi = 2 * M_PI*r1;
	float x = cos(phi) * 2 * sqrt(r2);
	float y = sin(phi) * 2 * sqrt(r2);
	return vec3(x, y, z);

}

// Ortho-normal Bases
class onb {
public:
	onb() {}
	inline vec3& operator[](int i) { return m_axis[i]; }
	inline const vec3& operator[](int i) const { return m_axis[i]; }
	const vec3& u() const { return m_axis[0]; }
	const vec3& v() const { return m_axis[1]; }
	const vec3& w() const { return m_axis[2]; }
	vec3 local(float a, float b, float c) const { return a*u() + b*v() + c*w(); }
	vec3 local(const vec3& a) const { return a.getX()*u() + a.getY()*v() + a.getZ()*w(); }
	void build_from_w(const vec3& n);
private:
	vec3 m_axis[3];
};

void onb::build_from_w(const vec3& n)
{
	m_axis[2] = normalize(n);
	vec3 a;
	if (fabs(w().getX()) > 0.9) {
		a = vec3(0,1,0);
	}
	else {
		a = vec3(1,0,0);
	}
	m_axis[1] = normalize(cross(w(), a));
	m_axis[0] = cross(w(), v());
}

