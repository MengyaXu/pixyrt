#pragma once

#include "ray.h"

inline vec3 random_in_unit_disk() 
{
	vec3 p;
	do {
		p = 2.0 * vec3(drand48(), drand48(), 0) - vec3(1,1,0);
	} while (dot(p,p) >= 1.0);
	return p;
}

class camera {
public:
	camera() {}
	camera(const vec3& lookfrom, const vec3& lookat, const vec3& vup, float vfov, float aspect, float aperture, float focus_dist, float t0, float t1) {
		m_time0 = t0;
		m_time1 = t1;
		m_lens_radius = aperture / 2;
		float theta = radians(vfov);
		float halfH = tan(theta/2);
		float halfW = aspect * halfH;
		m_origin = lookfrom;
		m_w = normalize(lookfrom - lookat);
		m_u = normalize(cross(vup, m_w));
		m_v = cross(m_w,m_u);
		m_lower_left_corner = m_origin - focus_dist*halfW*m_u - focus_dist*halfH*m_v - focus_dist*m_w;
		m_horizontal = 2*halfW*focus_dist*m_u;
		m_vertical = 2*halfH*focus_dist*m_v;
	}

	ray get_ray(float s, float t) const {
		vec3 rd = m_lens_radius * random_in_unit_disk();
		vec3 offset = m_u*rd.getX() + m_v*rd.getY();
		float time = m_time0 + drand48()*(m_time1-m_time0);
		return ray(m_origin + offset, m_lower_left_corner + s*m_horizontal + t*m_vertical - m_origin - offset, time);
	}

private:
	vec3 m_origin;
	vec3 m_lower_left_corner;
	vec3 m_horizontal;
	vec3 m_vertical;
	vec3 m_u, m_v, m_w;
	float m_time0, m_time1;
	float m_lens_radius;
};
