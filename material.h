#pragma once

#include "ray.h"
#include "hitable.h"
#include "texture.h"
#include "onb.h"
#include "pdf.h"

struct scatter_record
{
	vec3 albedo;
	ray  specular_ray;
	bool is_specular;
	pdf* pdfp;
};

class material {
public:
	virtual bool scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const = 0;
	virtual float scattering_pdf(const ray& r, const hit_record& hrec, scatter_record& srec) const { return 1; }
	virtual vec3 emitted(const ray& r, const hit_record& hrec, float u, float v, const vec3& p) const { return vec3(0); }
};

class diffuse_light : public material {
public:
	diffuse_light(texture* a) : m_emit(a) {}
	
	virtual bool scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const override;
	virtual vec3 emitted(const ray& r, const hit_record& hrec, float u, float v, const vec3& p) const override;

private:
	texture* m_emit;
};

bool diffuse_light::scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const
{
	return false;
}

vec3 diffuse_light::emitted(const ray& r, const hit_record& hrec, float u, float v, const vec3& p) const
{
	if (dot(hrec.n, r.direction()) < 0) {
		return m_emit->value(u, v, p);
	}
	else {
		return vec3(0);
	}
}

class sky_light : public material {
public:
	sky_light() {}

	virtual bool scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const override {
		return false;
	}

	virtual vec3 emitted(const ray& r, const hit_record& hrec, float u, float v, const vec3& p) const override {
		vec3 n = normalize(r.direction());
		float t = 0.5f * (n.getY() + 1.0f);
		return (1.0f - t) * vec3(1.0f) + t*vec3(0.5f, 0.7f, 1.0f);
	}
};

class isotropic : public material {
public:
	isotropic(texture* a) : m_albedo(a) {}
	
	virtual bool scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const override {
		srec.specular_ray = ray(hrec.p, random_in_unit_sphere());
		srec.albedo = m_albedo->value(hrec.u, hrec.v, hrec.p);
		srec.is_specular = false;
		return true;
	}

private:
	texture* m_albedo;
};

class lambertian : public material {
public:
	lambertian(texture* a) : m_albedo(a) {}

	virtual bool scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const override;
	virtual float scattering_pdf(const ray& r, const hit_record& hrec, scatter_record& srec) const override;

private:
	texture* m_albedo;
};

bool lambertian::scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const
{
	onb uvw;
	uvw.build_from_w(hrec.n);
	vec3 direction = uvw.local(random_cosine_direction());
	srec.is_specular = false;
	srec.specular_ray = ray(hrec.p, normalize(direction), r.time());
	srec.albedo = m_albedo->value(hrec.u, hrec.v, hrec.p);
	srec.pdfp = new cosine_pdf(hrec.n);
	return true;
}

float lambertian::scattering_pdf(const ray& r, const hit_record& hrec, scatter_record& srec) const
{
	float cosine = dot(hrec.n, normalize(srec.specular_ray.direction()));
	if (cosine < 0) cosine = 0;
	return cosine / M_PI;
}

class emission : public lambertian {
public:
	emission(texture* a, texture* e) : lambertian(a), m_emit(e) {}

	virtual vec3 emitted(const ray& r, const hit_record& hrec, float u, float v, const vec3& p) const override {
		return m_emit->value(hrec.u, hrec.v, p);
	}

private:
	texture* m_emit;
};

class metal : public material {
public:
	metal(const vec3& a, float f) : m_albedo(a), m_fuzz(f<1?f:1) {}

	virtual bool scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const override;

private:
	vec3 m_albedo;
	float m_fuzz;
};

bool metal::scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const
{
	vec3 reflected = reflect(normalize(r.direction()), hrec.n);
	srec.specular_ray = ray(hrec.p, reflected + m_fuzz*random_in_unit_sphere());
	srec.albedo = m_albedo;
	srec.is_specular = true;
	srec.pdfp = nullptr;
	return true;
	//return dot(srec.specular_ray.direction(), hrec.n) > 0;
}

class dielectric : public material {
public:
	dielectric(float roi) : m_roi(roi) {}

	virtual bool scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const override;

private:
	float m_roi;
};

bool dielectric::scatter(const ray& r, const hit_record& hrec, scatter_record& srec) const
{
	srec.albedo = vec3(1, 1, 1);
	srec.is_specular = true;
	srec.pdfp = nullptr;
	vec3 outward_normal;
	vec3 reflected = reflect(r.direction(), hrec.n);
	vec3 refracted;
	float ni_over_nt;
	float reflect_prob;
	float cosine;
	if (dot(r.direction(), hrec.n) > 0) {
		outward_normal = -hrec.n;
		ni_over_nt = m_roi;
		cosine = m_roi * dot(r.direction(), hrec.n) / length(r.direction());
	}
	else {
		outward_normal = hrec.n;
		ni_over_nt = 1.0 / m_roi;
		cosine = -dot(r.direction(), hrec.n) / length(r.direction());
	}

	if (refract(r.direction(), outward_normal, ni_over_nt, refracted)) {
		reflect_prob = schlick(cosine, m_roi);
	}
	else {
		reflect_prob = 1;
	}

	if (drand48() < reflect_prob) {
		srec.specular_ray = ray(hrec.p, reflected);
	}
	else {
		srec.specular_ray = ray(hrec.p, refracted);
	}

	return true;
}