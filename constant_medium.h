#pragma once
#include "hitable.h"
#include "material.h"

class constant_medium : public hitable {
public:
	constant_medium(hitable* p, float d, texture* a)
		: m_boundary(p)
		, m_density(d) 
	{
		m_phase_function = new isotropic(a);
	}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override {
		return m_boundary->bounds(t0, t1, box);
	}

private:
	hitable* m_boundary;
	float m_density;
	material* m_phase_function;
};

bool constant_medium::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	bool db = (drand48() < 0.00001);
	db = false;
	hit_record hrec1, hrec2;
	if (m_boundary->hit(r, -FLT_MAX, FLT_MAX, hrec1)) {
		if (m_boundary->hit(r, hrec1.t + 0.0001, FLT_MAX, hrec2)) {
			if (db) std::cerr << "\nt0 t1 " << hrec1.t << " " << hrec2.t << "\n";
			hrec1.t = Max(hrec1.t, tmin);
			hrec2.t = Min(hrec2.t, tmax);
			if (hrec1.t >= hrec2.t) {
				return false;
			}

			hrec1.t = Max(hrec1.t, 0.0f);

			float len = length(r.direction());
			float distance_inside_boundary = (hrec2.t - hrec1.t) * len;
			float hit_distance = -(1/m_density) * log(drand48());
			if (hit_distance < distance_inside_boundary) {
				if (db) std::cerr << "hit distance = " << hit_distance << "\n";
				hrec.t = hrec1.t + hit_distance / len;
				if (db) std::cerr << "rec.t = " << hrec.t << "\n";
				hrec.p = r.point_at_parameter(hrec.t);
				if (db) std::cerr << "rec.p = " << hrec.p.getX() << "," << hrec.p.getY() << "," << hrec.p.getZ() << "\n";
				hrec.n = vec3(1,0,0); // arbitrary
				hrec.matp = m_phase_function;
				return true;
			}
		}
	}

	return false;
}
