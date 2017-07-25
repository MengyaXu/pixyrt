#pragma once
#include "vec.h"
#include "onb.h"

class pdf {
public:
	virtual float value(const vec3& direction) const = 0;
	virtual vec3 generate() const = 0;
};

class cosine_pdf : public pdf {
public:
	cosine_pdf(const vec3& w) { m_uvw.build_from_w(w); }

	virtual float value(const vec3& direction) const {
		float cosine = dot(normalize(direction), m_uvw.w());
		if (cosine > 0) {
			return cosine / M_PI;
		}
		else {
			return 0;
		}
	}

	virtual vec3 generate() const {
		return m_uvw.local(random_cosine_direction());
	}

private:
	onb m_uvw;
};

class hitable_pdf : public pdf {
public:
	hitable_pdf(hitable* p, const vec3& origin) : m_hitable(p), m_origin(origin) {}

	virtual float value(const vec3& direction) const override {
		return m_hitable->pdf_value(m_origin, direction);
	}

	virtual vec3 generate() const override {
		return m_hitable->random(m_origin);
	}

private:
	hitable* m_hitable;
	vec3 m_origin;
};

class mixture_pdf : public pdf {
public:
	mixture_pdf(pdf* p0, pdf* p1) { m_pdfs[0] = p0; m_pdfs[1] = p1; }

	virtual float value(const vec3& direction) const override {
		return 0.5 * m_pdfs[0]->value(direction) + 0.5 * m_pdfs[1]->value(direction);
	}

	virtual vec3 generate() const override {
		if (drand48() < 0.5) {
			return m_pdfs[0]->generate();
		}
		else {
			return m_pdfs[1]->generate();
		}
	}

private:
	pdf* m_pdfs[2];
};