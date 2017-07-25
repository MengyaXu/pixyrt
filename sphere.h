#pragma once

#include "hitable.h"
#include "onb.h"

class sphere : public hitable {
public:
	sphere() {}
	sphere(const vec3& c, float r, material* matp) : m_center(c), m_radius(r), m_matp(matp) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override;
	virtual float pdf_value(const vec3& o, const vec3& v) const override;
	virtual vec3 random(const vec3& o) const override;

private:
	vec3 m_center;
	float m_radius;
	material* m_matp;
};

bool sphere::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	vec3 oc = r.origin() - m_center;
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc,oc) - m_radius*m_radius;
	float D = b*b - a*c;
	if (D > 0) {
		float sqrtD = sqrt(D);
		float t = (-b - sqrtD) / a;
		if (t < tmax && t > tmin) {
			hrec.t = t;
			hrec.p = r.point_at_parameter(t);
			hrec.n = (hrec.p - m_center) / m_radius;
			hrec.matp = m_matp;
			get_sphere_uv((hrec.p-m_center)/m_radius, hrec.u, hrec.v);
			return true;
		}

		t = (-b + sqrtD) / a;
		if (t < tmax && t > tmin) {
			hrec.t = t;
			hrec.p = r.point_at_parameter(t);
			hrec.n = (hrec.p - m_center) / m_radius;
			hrec.matp = m_matp;
			get_sphere_uv((hrec.p - m_center) / m_radius, hrec.u, hrec.v);
			return true;
		}
	}
	return false;
}

bool sphere::bounds(float t0, float t1, aabb& box) const
{
	box = aabb(m_center - vec3(m_radius), m_center + vec3(m_radius));
	return true;
}

float sphere::pdf_value(const vec3& o, const vec3& v) const
{
	hit_record hrec;
	if (this->hit(ray(o,v), 0.001, FLT_MAX, hrec)) {
		float cos_theta_max = sqrt(1-pow2(m_radius) * recip(lengthSqr(m_center - o)));
		float solid_angle = 2 * M_PI * (1-cos_theta_max);
		return recip(solid_angle);
	}
	else {
		return 0;
	}
}

vec3 sphere::random(const vec3& o) const
{
	vec3 direction = m_center - o;
	float distance_squared = lengthSqr(direction);
	onb uvw;
	uvw.build_from_w(direction);
	return uvw.local(random_to_sphere(m_radius, distance_squared));
}


class moving_sphere : public hitable {
public:
	moving_sphere() {}
	moving_sphere(const vec3& center0, const vec3& center1, float t0, float t1, float r, material* matp)
		: m_center0(center0), m_center1(center1), m_time0(t0), m_time1(t1), m_radius(r), m_matp(matp) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& aabb) const override;

	vec3 center(float time) const;

private:
	vec3 m_center0, m_center1;
	float m_time0, m_time1;
	float m_radius;
	material* m_matp;
};

vec3 moving_sphere::center(float time) const
{
	return m_center0 + ((time - m_time0) / (m_time1 - m_time0)) * (m_center1 - m_center0);
}

bool moving_sphere::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	vec3 oc = r.origin() - center(r.time());
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc, oc) - m_radius*m_radius;
	float D = b*b - a*c;
	if (D > 0) {
		float sqrtD = sqrt(D);
		float t = (-b - sqrtD) / a;
		if (t < tmax && t > tmin) {
			hrec.t = t;
			hrec.p = r.point_at_parameter(t);
			hrec.n = (hrec.p - center(r.time())) / m_radius;
			hrec.matp = m_matp;
			return true;
		}

		t = (-b + sqrtD) / a;
		if (t < tmax && t > tmin) {
			hrec.t = t;
			hrec.p = r.point_at_parameter(t);
			hrec.n = (hrec.p - center(r.time())) / m_radius;
			hrec.matp = m_matp;
			return true;
		}
	}
	return false;
}

bool moving_sphere::bounds(float t0, float t1, aabb& box) const
{
	box = surrounding_box(
		aabb(center(t0) - vec3(m_radius), center(t0) + vec3(m_radius)),
		aabb(center(t1) - vec3(m_radius), center(t1) + vec3(m_radius)));
	return true;
}