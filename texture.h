#pragma once

#include "vec.h"
#include "perlin.h"

class texture {
public:
	virtual vec3 value(float u, float v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public:
	constant_texture() {}
	constant_texture(const vec3& c) : m_color(c) {}

	virtual vec3 value(float u, float v, const vec3& p) const override {
		return m_color;
	}

private:
	vec3 m_color;
};

class checker_texture : public texture {
public:
	checker_texture() {}
	checker_texture(texture* t0, texture* t1) : m_even(t0), m_odd(t1) {}

	virtual vec3 value(float u, float v, const vec3& p) const override {
		float sines = sin(1*p.getX()) * sin(1*p.getY()) * sin(1*p.getZ());
		if (sines < 0) {
			return m_odd->value(u, v, p);
		}
		else {
			return m_even->value(u,v,p);
		}
	}

private:
	texture* m_odd;
	texture* m_even;
};

class noise_texture : public texture {
public:
	noise_texture() {}
	noise_texture(float s) : m_scale(s) {}
	
	virtual vec3 value(float u, float v, const vec3& p) const override {
		//return vec3(1) * m_noise.noise(m_scale * p);
		//return vec3(1) * m_noise.turb(m_scale * p);
		return vec3(0.5) * (1 + sin(m_scale*p.getX() + 10*m_noise.turb(p)));
	}

private:
	perlin m_noise;
	float m_scale;
};

class image_texture : public texture {
public:
	image_texture() {}
	image_texture(unsigned char* pixels, int w, int h) : m_data(pixels), m_nx(w), m_ny(h) {}

	virtual vec3 value(float u, float v, const vec3& p) const override {
		int i = (  u)*m_nx;
		int j = (1-v)*m_ny-0.001;
		//if (i < 0) i = 0;
		//if (j < 0) j = 0;
		//if (i > m_nx-1) i = m_nx-1;
		//if (j > m_ny-1) j = m_ny-1;

		// bilinear
		float ucoef = m_nx*u - int(m_nx*u);
		float vcoef = m_ny*v - int(m_ny*v);
		vec3 c0 = sample(i,j);
		vec3 c1 = sample(i+1,j);
		vec3 c2 = sample(i,j+1);
		vec3 c3 = sample(i+1,j+1);
		vec3 c01 = lerp(ucoef, c0, c1);
		vec3 c23 = lerp(ucoef, c2, c3);
		return lerp(vcoef, c01, c23);

		//float r = int(m_data[3*i + 3*m_nx*j]  ) / 255.0;
		//float g = int(m_data[3*i + 3*m_nx*j+1]) / 255.0;
		//float b = int(m_data[3*i + 3*m_nx*j+2]) / 255.0;
		//return vec3(r,g,b);
	}

	vec3 sample(int u, int v) const
	{
		u = u<0 ? 0 : u >= m_nx ? m_nx - 1 : u;
		v = v<0 ? 0 : v >= m_ny ? m_ny - 1 : v;
		return vec3(
			int(m_data[3 * u + 3 * m_nx * v]) / 255.0,
			int(m_data[3 * u + 3 * m_nx * v+1]) / 255.0,
			int(m_data[3 * u + 3 * m_nx * v+2]) / 255.0);
	}

private:
	unsigned char* m_data;
	int m_nx;
	int m_ny;
};

class cube_texture {
public:
	enum { 
		kFacePositiveX = 0,
		kFaceNegativeX,
		kFacePositiveY,
		kFaceNegativeY,
		kFacePositiveZ,
		kFaceNegativeZ,
		kFaceNum
	};
	cube_texture() {}
	cube_texture(texture* px, texture* nx, texture* py, texture* ny, texture* pz, texture* nz) {
		m_faces[kFacePositiveX] = px;
		m_faces[kFaceNegativeX] = nx;
		m_faces[kFacePositiveY] = py;
		m_faces[kFaceNegativeY] = ny;
		m_faces[kFacePositiveZ] = pz;
		m_faces[kFaceNegativeZ] = nz;
	};

	vec3 value(const vec3& n) const
	{
		vec3 absNrm = absPerElem(n);
		float nx = n.getX();
		float ny = -n.getY();
		float nz = -n.getZ();
		float x = absNrm.getX();
		float y = absNrm.getY();
		float z = absNrm.getZ();
		if (x >= y && x >= z) {
			if (nx > 0.0f) {
				return m_faces[kFacePositiveX]->value(
					1.0f - (nz / nx + 1.0f) * 0.5f,
					(ny / nx + 1.0f) * 0.5f,
					n);
			}
			else if (nx < 0.0f) {
				return m_faces[kFaceNegativeX]->value(
					1.0f - (nz / nx + 1.0f) * 0.5f,
					1.0f - (ny / nx + 1.0f) * 0.5f,
					n);
			}
		}
		else if (y >= x && y >= z) {
			if (ny > 0.0f) {
				return m_faces[kFacePositiveY]->value(
					(nx/ny+1.0f)*0.5f,
					1.0f - (nz/ny+1.0f)*0.5f,
					n);
			}
			else if (ny < 0.0f) {
				return m_faces[kFaceNegativeY]->value(
					1.0f - (nx/ny+1.0f)*0.5f,
					(nz/ny+1.0f)*0.5f,
					n);
			}
		}
		else if (z >= x && z >= y) {
			if (nz > 0.0f) {
				return m_faces[kFacePositiveZ]->value(
					(nx/nz+1.0f)*0.5f,
					(ny/nz+1.0f)*0.5f,
					n);
			}
			else if (nz < 0.0f) {
				return m_faces[kFaceNegativeZ]->value(
					(nx/nz+1.0f)*0.5f,
					1.0f - (ny/nz+1.0f)*0.5f,
					n);
			}
		}

		return vec3(0);
	}

private:
	texture* m_faces[kFaceNum];
};
