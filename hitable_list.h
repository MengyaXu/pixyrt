#pragma once

#include "hitable.h"

class hitable_list : public hitable {
public:
	hitable_list() {}
	hitable_list(hitable** l, int n) : m_list(l), m_list_size(n) {}

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override;
	virtual float pdf_value(const vec3& o, const vec3& v) const override;
	virtual vec3 random(const vec3& o) const override;

private:
	hitable** m_list;
	int m_list_size;
};

bool hitable_list::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	hit_record rec;
	bool hit_anything = false;
	float closest_so_far = tmax;
	for (int i=0; i<m_list_size; ++i) {
		if (m_list[i]->hit(r, tmin, closest_so_far, rec)) {
			hit_anything = true;
			closest_so_far = rec.t;
			hrec = rec;
		}
	}
	return hit_anything;
}

bool hitable_list::bounds(float t0, float t1, aabb& box) const
{
	if (m_list_size < 1) return false;

	aabb tmp;
	if (m_list[0]->bounds(t0, t1, tmp) == false) {
		return false;
	}

	box = tmp;

	for (int i=1; i<m_list_size; ++i) {
		if (m_list[i]->bounds(t0, t1, tmp)) {
			box = surrounding_box(box, tmp);
		}
		//else {
		//	return false;
		//}
	}
	return true;
}

float hitable_list::pdf_value(const vec3& o, const vec3& v) const
{
	//float weight = 1.0 / 1;
	//float sum = 0;
	//for (int i=0; i<1; ++i) {
	//	sum += weight * m_list[i]->pdf_value(o,v);
	//}
	float weight = 1.0 / m_list_size;
	float sum = 0;
	for (int i=0; i<m_list_size; ++i) {
		sum += weight * m_list[i]->pdf_value(o,v);
	}
	return sum;
}

vec3 hitable_list::random(const vec3& o) const
{
	//return m_list[0]->random(o);
	int index = int(drand48() * m_list_size);
	if (index >= m_list_size) index = m_list_size-1;
	return m_list[index]->random(o);
}