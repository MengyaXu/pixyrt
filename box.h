#pragma once
#include "hitable_list.h"
#include "rect.h"

class box : public hitable {
public:
	box() {}
	box(const vec3& p0, const vec3& p1, material* matp);

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const { box = aabb(m_p0, m_p1); return true; }

private:
	vec3 m_p0, m_p1;
	material* m_matp;
	hitable* m_listp;
};

box::box(const vec3& p0, const vec3& p1, material* matp)
	: m_p0(p0)
	, m_p1(p1)
	, m_matp(matp)
	, m_listp(nullptr)
{
	hitable** list = new hitable*[6];
	list[0] = new xy_rect(m_p0.getX(), m_p1.getX(), m_p0.getY(), m_p1.getY(), m_p1.getZ(), m_matp);
	list[1] = new flip_normals(new xy_rect(m_p0.getX(), m_p1.getX(), m_p0.getY(), m_p1.getY(), m_p0.getZ(), m_matp));
	list[2] = new xz_rect(m_p0.getX(), m_p1.getX(), m_p0.getZ(), m_p1.getZ(), m_p1.getY(), m_matp);
	list[3] = new flip_normals(new xz_rect(m_p0.getX(), m_p1.getX(), m_p0.getZ(), m_p1.getZ(), m_p0.getY(), m_matp));
	list[4] = new yz_rect(m_p0.getY(), m_p1.getY(), m_p0.getZ(), m_p1.getZ(), m_p1.getX(), m_matp);
	list[5] = new flip_normals(new yz_rect(m_p0.getY(), m_p1.getY(), m_p0.getZ(), m_p1.getZ(), m_p0.getX(), m_matp));
	m_listp = new hitable_list(list, 6);
}

bool box::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	return m_listp->hit(r, tmin, tmax, hrec);
}