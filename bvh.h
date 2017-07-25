#pragma once

#include "hitable.h"
#include <iostream>

int box_x_compare(const void* a, const void* b)
{
	aabb lbox, rbox;
	hitable* ah = *(hitable**)a;
	hitable* bh = *(hitable**)b;
	if (!ah->bounds(0,0,lbox) || !bh->bounds(0,0,rbox)) {
		std::cerr << "no bounding box in bvh_node constructor" << std::endl;
	}
	if (lbox.lower().getX() - rbox.lower().getX() < 0.0) {
		return -1;
	}
	else {
		return 1;
	}
}

int box_y_compare(const void* a, const void* b)
{
	aabb lbox, rbox;
	hitable* ah = *(hitable**)a;
	hitable* bh = *(hitable**)b;
	if (!ah->bounds(0, 0, lbox) || !bh->bounds(0, 0, rbox)) {
		std::cerr << "no bounding box in bvh_node constructor" << std::endl;
	}
	if (lbox.lower().getY() - rbox.lower().getY() < 0.0) {
		return -1;
	}
	else {
		return 1;
	}
}

int box_z_compare(const void* a, const void* b)
{
	aabb lbox, rbox;
	hitable* ah = *(hitable**)a;
	hitable* bh = *(hitable**)b;
	if (!ah->bounds(0, 0, lbox) || !bh->bounds(0, 0, rbox)) {
		std::cerr << "no bounding box in bvh_node constructor" << std::endl;
	}
	if (lbox.lower().getZ() - rbox.lower().getZ() < 0.0) {
		return -1;
	}
	else {
		return 1;
	}
}

class bvh_node : public hitable {
public:
	bvh_node() {}
	bvh_node(hitable** l, int n, float t0, float t1);

	virtual bool hit(const ray& r, float tmin, float tmax, hit_record& hrec) const override;
	virtual bool bounds(float t0, float t1, aabb& box) const override;

private:
	hitable* m_left;
	hitable* m_right;
	aabb m_bounds;
};

bvh_node::bvh_node(hitable** l, int n, float t0, float t1)
{
	int axis = int (3*drand48());
	if (axis == 0) {
		qsort(l, n, sizeof(hitable*), box_x_compare);
	}
	else if (axis == 1) {
		qsort(l, n, sizeof(hitable*), box_y_compare);
	}
	else {
		qsort(l, n, sizeof(hitable*), box_z_compare);
	}

	if (n == 1) {
		m_left = m_right = l[0];
	}
	else if (n == 2) {
		m_left = l[0];
		m_right = l[1];
	}
	else {
		m_left = new bvh_node(l, n/2, t0, t1);
		m_right = new bvh_node(l+n/2, n-n/2, t0, t1);
	}

	aabb lbox, rbox;
	if (!m_left->bounds(t0, t1, lbox) || !m_right->bounds(t0, t1, rbox)) {
		std::cerr << "no bounding box in bvh_node constructor\n";
	}

	m_bounds = surrounding_box(lbox, rbox);
}

bool bvh_node::hit(const ray& r, float tmin, float tmax, hit_record& hrec) const
{
	if (m_bounds.hit(r, tmin, tmax)) {
		hit_record lrec, rrec;
		bool lhit = m_left->hit(r, tmin, tmax, lrec);
		bool rhit = m_right->hit(r, tmin, tmax, rrec);
		if (lhit && rhit) {
			if (lrec.t < rrec.t) {
				hrec = lrec;
			}
			else {
				hrec = rrec;
			}
			return true;
		}
		else if (lhit) {
			hrec = lrec;
			return true;
		}
		else if (rhit) {
			hrec = rrec;
			return true;
		}
		else {
			return false;
		}
	}
	else {
		return false;
	}
}

bool bvh_node::bounds(float t0, float t1, aabb& box) const
{
	box = m_bounds;
	return true;
}
