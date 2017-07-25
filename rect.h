#pragma once

#include "hitable.h"

class xy_rect : public hitable {
public:
	xy_rect() {}
	xy_rect(float x0, float x1, float y0, float y1, float k, material* matp) : m_x0(x0), m_x1(x1), m_y0(y0), m_y1(y1), m_k(k), m_matp(matp) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override;

private:
	material* m_matp;
	float m_x0, m_x1, m_y0, m_y1, m_k;
};

class xz_rect : public hitable {
public:
	xz_rect() {}
	xz_rect(float x0, float x1, float z0, float z1, float k, material* matp) : m_x0(x0), m_x1(x1), m_z0(z0), m_z1(z1), m_k(k), m_matp(matp) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override;

	virtual float pdf_value(const vec3& o, const vec3& v) const override {
		hit_record hrec;
		if (this->hit(ray(o,v), 0.001, FLT_MAX, hrec)) {
			float area = (m_x1 - m_x0) * (m_z1 - m_z0);
			float distance_squared = pow2(hrec.t) * lengthSqr(v);
			float cosine = fabs(dot(v,hrec.n)) / (length(v));
			return distance_squared / (cosine * area);
		}
		else {
			return 0;
		}
	}

	virtual vec3 random(const vec3& o) const override {
		vec3 random_point = vec3(
			m_x0 + drand48()*(m_x1 - m_x0),
			m_k,
			m_z0 + drand48()*(m_z1 - m_z0));
		return random_point - o;
	}

private:
	material* m_matp;
	float m_x0, m_x1, m_z0, m_z1, m_k;
};

class yz_rect : public hitable {
public:
	yz_rect() {}
	yz_rect(float y0, float y1, float z0, float z1, float k, material* matp) : m_y0(y0), m_y1(y1), m_z0(z0), m_z1(z1), m_k(k), m_matp(matp) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override;

private:
	material* m_matp;
	float m_y0, m_y1, m_z0, m_z1, m_k;
};

bool xy_rect::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	float t = (m_k - r.origin().getZ()) / r.direction().getZ();
	if (t < tmin || t > tmax) {
		return false;
	}
	float x = r.origin().getX() + t*r.direction().getX();
	float y = r.origin().getY() + t*r.direction().getY();
	if (x < m_x0 || x > m_x1 || y < m_y0 || y > m_y1) {
		return false;
	}

	hrec.u = (x - m_x0) / (m_x1 - m_x0);
	hrec.v = (y - m_y0) / (m_y1 - m_y0);
	hrec.t = t;
	hrec.matp = m_matp;
	hrec.p = r.point_at_parameter(t);
	hrec.n = vec3(0, 0, 1);
	return true;
}

bool xy_rect::bounds(float t0, float t1, aabb& box) const
{
	box = aabb(vec3(m_x0, m_y0, m_k - 0.0001), vec3(m_x1, m_y1, m_k + 0.0001));
	return true;
}

bool xz_rect::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	float t = (m_k - r.origin().getY()) / r.direction().getY();
	if (t < tmin || t > tmax) {
		return false;
	}
	float x = r.origin().getX() + t*r.direction().getX();
	float z = r.origin().getZ() + t*r.direction().getZ();
	if (x < m_x0 || x > m_x1 || z < m_z0 || z > m_z1) {
		return false;
	}

	hrec.u = (x - m_x0) / (m_x1 - m_x0);
	hrec.v = (z - m_z0) / (m_z1 - m_z0);
	hrec.t = t;
	hrec.matp = m_matp;
	hrec.p = r.point_at_parameter(t);
	hrec.n = vec3(0, 1, 0);
	return true;
}

bool xz_rect::bounds(float t0, float t1, aabb& box) const
{
	box = aabb(vec3(m_x0, m_k - 0.0001, m_z0), vec3(m_x1, m_k + 0.0001, m_z1));
	return true;
}

bool yz_rect::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	float t = (m_k - r.origin().getX()) / r.direction().getX();
	if (t < tmin || t > tmax) {
		return false;
	}
	float y = r.origin().getY() + t*r.direction().getY();
	float z = r.origin().getZ() + t*r.direction().getZ();
	if (y < m_y0 || y > m_y1 || z < m_z0 || z > m_z1) {
		return false;
	}

	hrec.u = (y - m_y0) / (m_y1 - m_y0);
	hrec.v = (z - m_z0) / (m_z1 - m_z0);
	hrec.t = t;
	hrec.matp = m_matp;
	hrec.p = r.point_at_parameter(t);
	hrec.n = vec3(1, 0, 0);
	return true;
}

bool yz_rect::bounds(float t0, float t1, aabb& box) const
{
	box = aabb(vec3(m_k - 0.0001, m_y0, m_z0), vec3(m_k + 0.0001, m_y1, m_z1));
	return true;
}