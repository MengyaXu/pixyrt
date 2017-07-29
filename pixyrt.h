/*
* Copyright (c) 2017 mebiusbox software. All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
* OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
* HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
* OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
* SUCH DAMAGE.
*/
#pragma once

#include "vec.h"
#include "camera.h"
#include "texture.h"

#define Max(a,b) ((a)>(b))?(a):(b)
#define Min(a,b) ((a)<(b))?(a):(b)

class ReflectedLight {
public:
	ReflectedLight()
		: directDiffuse(0.0f)
		, directSpecular(0.0f)
		, indirectDiffuse(0.0f)
		, indirectSpecular(0.0f)
	{
	}

	Vector3 getColor() {
		return directDiffuse + directSpecular + indirectDiffuse + indirectSpecular;
	}

	Vector3 directDiffuse;
	Vector3 directSpecular;
	Vector3 indirectDiffuse;
	Vector3 indirectSpecular;
};

struct GeometricContext {
	Vector3 position;
	Vector3 normal;
	Vector3 viewDir;
	Vector2 uv0;
};

struct Material {
	Vector3 diffuseColor;
	Vector3 specularColor;
	float specularRoughness;
	int diffuseMethod;
	Vector3 albedo;
};

struct IncidentLight {
	Vector3 color;
	Vector3 direction;
};

struct Samplers {
	texture* albedo;
	texture* roughness;
	cube_texture* envmap;
};

typedef struct _sphere {
	Vector3 color;
	Vector3 center;
	float radius;
	texture* albedo;
	texture* roughness;
	cube_texture* envmap;
} sphere_t;

class Scene {
public:
	Scene() {}

	void setCamera(const camera& cam) { m_camera = cam; }
	void setBgColor(const vec3& col) { m_bgcolor = col; }
	void setBall(const sphere_t& ball) { m_ball = ball; }
	void setLight(const vec3& pos) { m_light = pos; }
	void setSize(int w, int h) { m_width = w; m_height = h; }

	const camera& getCamera() const { return m_camera; }
	const vec3& getBgColor() const { return m_bgcolor; }
	const sphere_t getBall() const { return m_ball; }
	const vec3& getLight() const { return m_light; }
	int getWidth() const { return m_width; }
	int getHeight() const { return m_height; }

private:
	vec3 m_bgcolor;
	camera m_camera;
	sphere_t m_ball;
	vec3 m_light;
	int m_width;
	int m_height;
};

class Image {
public:
	struct rgb {
		unsigned char r;
		unsigned char g;
		unsigned char b;
		//unsigned char a;
	};

	Image() : m_pixels(nullptr) { resetFilter(); }
	Image(int w, int h) {
		m_width = w;
		m_height = h;
		m_pixels = new rgb[m_width*m_height];
		resetFilter();
	}
	~Image() {
		if (m_pixels) delete[] m_pixels;
	}

	void resetFilter() {
		m_filter_gamma = true;
		m_filter_denan = true;
		m_filter_saturate = true;
	}

	void setFilterGamma(bool flag) { m_filter_gamma = flag; }
	void setFilterDenan(bool flag) { m_filter_denan = flag; }
	void setFilterSaturate(bool flag) { m_filter_saturate = flag; }

	int width() const { return m_width; }
	int height() const { return m_height; }
	void* pixels() const { return m_pixels; }

	void clear(const vec3& color)
	{
		for (int i=0; i<m_height; ++i) {
			for (int j=0; j<m_width; ++j) {
				write(j, i, color);
			}
		}
	}

	void write(int x, int y, const vec3& color)
	{
		vec3 c = filter(color) * 255.99f;
		int index = m_width*y + x;
		m_pixels[index].r = static_cast<unsigned char>(c.getX());
		m_pixels[index].g = static_cast<unsigned char>(c.getY());
		m_pixels[index].b = static_cast<unsigned char>(c.getZ());
	}

	vec3 filter(const vec3& color) {
		vec3 c = color;
		if (m_filter_gamma) {
			c = linear_to_gamma(c, GAMMA_FACTOR);
		}
		if (m_filter_denan) {
			c = de_nan(c);
		}
		if (m_filter_saturate) {
			c = saturate(c);
		}
		return c;
	}

private:
	int m_width;
	int m_height;
	rgb* m_pixels;

	bool m_filter_gamma;
	bool m_filter_denan;
	bool m_filter_saturate;
};

class Workarea {
public:
	Workarea() {}
	Workarea(int x, int y, int w, int h, int outx, int outy) : m_x(x), m_y(y), m_width(w), m_height(h), m_outx(outx), m_outy(outy) {}

	void setPos(int x, int y) { m_x = x; m_y = y; }
	void setOutPos(int x, int y) { m_outx = x; m_outy = y; }
	void setSize(int w, int h) { m_width = w; m_height = h; }

	int x() const { return m_x; }
	int y() const { return m_y; }
	int w() const { return m_width; }
	int h() const { return m_height; }
	int outx() const { return m_outx; }
	int outy() const { return m_outy; }

private:
	int m_x;
	int m_y;
	int m_width;
	int m_height;
	int m_outx;
	int m_outy;
};