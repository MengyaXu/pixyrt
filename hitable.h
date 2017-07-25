#pragma once

#include "ray.h"
#include "bounds.h"

class material;

vec3 random_in_unit_sphere() {
	vec3 p;
	do {
		p = 2.0*rndvec() - vec3(1);
	} while (lengthSqr(p) >= 1.0);
	return p;
}

vec3 random_on_unit_sphere() {
	return normalize(random_in_unit_sphere());
}

inline vec3 random_to_sphere(float radius, float distance_squared) {
	float r1 = drand48();
	float r2 = drand48();
	float z = 1 + r2 * (sqrt(1-pow2(radius) * recip(distance_squared)) - 1);
	float phi = 2*M_PI*r1;
	float sqrt_z = sqrt(1-z*z);
	float x = cos(phi) * sqrt_z;
	float y = sin(phi) * sqrt_z;
	return vec3(x,y,z);
}

struct hit_record {
	float t;
	float u;
	float v;
	vec3 p;
	vec3 n;
	material* matp;
};

class hitable {
public:
	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const = 0;
	virtual bool bounds(float t0, float t1, aabb& box) const = 0;
	virtual float pdf_value(const vec3& o, const vec3& v) const { return 0; }
	virtual vec3 random(const vec3& o) const { return vec3(1,0,0); }
};

class flip_normals : public hitable {
public:
	flip_normals(hitable* p) : m_hitable(p) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override {
		if (m_hitable->hit(r, tmin, tmax, hrec)) {
			hrec.n = -hrec.n;
			return true;
		}
		else {
			return false;
		}
	}
	virtual bool bounds(float t0, float t1, aabb& box) const override {
		return m_hitable->bounds(t0, t1, box);
	}

	virtual float pdf_value(const vec3& o, const vec3& v) const {
		return m_hitable->pdf_value(o,v);
	}

	virtual vec3 random(const vec3& o) const {
		return m_hitable->random(o);
	}

private:
	hitable* m_hitable;
};

class translate : public hitable {
public:
	translate(hitable* p, const vec3& displacement) : m_hitable(p), m_offset(displacement) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override {
		ray moved_r(r.origin() - m_offset, r.direction(), r.time());
		if (m_hitable->hit(moved_r, tmin, tmax, hrec)) {
			hrec.p += m_offset;
			return true;
		}
		else {
			return false;
		}
	}
	virtual bool bounds(float t0, float t1, aabb& box) const override {
		if (m_hitable->bounds(t0, t1, box)) {
			box = aabb(box.lower() + m_offset, box.upper() + m_offset);
			return true;
		}
		else {
			return false;
		}
	}

private:
	hitable* m_hitable;
	vec3 m_offset;
};

class rotate_y : public hitable {
public:
	rotate_y(hitable* p, float angle);

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override {
		box = m_bounds;
		return m_has_box;
	}

	vec3 rotate_vec(const vec3& v) const
	{
		return vec3(
			m_cos_theta*v.getX() + m_sin_theta*v.getZ(),
			v.getY(),
			-m_sin_theta*v.getX() + m_cos_theta*v.getZ());
	}

	vec3 reverse_rotate_vec(const vec3& v) const
	{
		return vec3(
			m_cos_theta*v.getX() - m_sin_theta*v.getZ(),
			v.getY(),
			m_sin_theta*v.getX() + m_cos_theta*v.getZ());
	}

private:
	hitable* m_hitable;
	float m_sin_theta;
	float m_cos_theta;
	bool m_has_box;
	aabb m_bounds;
};

rotate_y::rotate_y(hitable* p, float angle)
	: m_hitable(p)
{
	float rad = radians(angle);
	m_sin_theta = sin(rad);
	m_cos_theta = cos(rad);
	m_has_box = p->bounds(0, 1, m_bounds);

	vec3 lower(FLT_MAX, FLT_MAX, FLT_MAX);
	vec3 upper(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	for (int i=0; i<2; ++i) {
		for (int j=0; j<2; ++j) {
			for (int k=0; k<2; ++k) {
				float x = i*m_bounds.upper().getX() + (1 - i) * m_bounds.lower().getX();
				float y = j*m_bounds.upper().getY() + (1 - j) * m_bounds.lower().getY();
				float z = k*m_bounds.upper().getZ() + (1 - k) * m_bounds.lower().getZ();
				vec3 tester = rotate_vec(vec3(x,y,z));
				for (int c=0; c<3; ++c) {
					if (tester[c] > upper[c]) {
						upper[c] = tester[c];
					}
					if (tester[c] < lower[c]) {
						lower[c] = tester[c];
					}
				}
			}
		}
	}

	m_bounds = aabb(lower, upper);
}

bool rotate_y::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	vec3 origin = reverse_rotate_vec(r.origin());
	vec3 direction = reverse_rotate_vec(r.direction());
	ray rotated_r(origin, direction, r.time());
	if (m_hitable->hit(rotated_r, tmin, tmax, hrec)) {
		hrec.p = rotate_vec(hrec.p);
		hrec.n = rotate_vec(hrec.n);
		return true;
	}
	else {
		return false;
	}
}