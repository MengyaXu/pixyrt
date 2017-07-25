#pragma once
#include "vec.h"

inline float trilinear_interp(float c[2][2][2], float u, float v, float w) {
	float accum = 0;
	for (int i=0; i<2; ++i) {
		for (int j=0; j<2; j++) {
			for (int k=0; k<2; k++) {
				accum += (i*u + (1 - i)*(1 - u)) * 
					     (j*v + (1 - j)*(1 - v)) *
					     (k*w + (1 - k)*(1 - w)) * c[i][j][k];
			}
		}
	}
	return accum;
}

inline float perlin_interp(vec3 c[2][2][2], float u, float v, float w) {
	float uu = smoothstep(0, 1, u);
	float vv = smoothstep(0, 1, v);
	float ww = smoothstep(0, 1, w);
	float accum = 0;
	for (int i = 0; i<2; ++i) {
		for (int j = 0; j<2; j++) {
			for (int k = 0; k<2; k++) {
				vec3 weight_v(u-i, v-j, w-k);
				accum += (i*uu + (1 - i)*(1 - uu)) *
					     (j*vv + (1 - j)*(1 - vv)) *
					     (k*ww + (1 - k)*(1 - ww)) * dot(c[i][j][k], weight_v);
			}
		}
	}
	return accum;
}

class perlin {
public:
	float noise(const vec3& p) const {
		float u = p.getX() - floor(p.getX());
		float v = p.getY() - floor(p.getY());
		float w = p.getZ() - floor(p.getZ());
		int i = floor(p.getX());
		int j = floor(p.getY());
		int k = floor(p.getZ());
		vec3 c[2][2][2];
		for (int di=0; di<2; ++di) {
			for (int dj=0; dj<2; ++dj) {
				for (int dk=0; dk<2; ++dk) {
					c[di][dj][dk] = m_ranvec[
						m_perm_x[(i+di)&255] ^ m_perm_y[(j+dj)&255] ^ m_perm_z[(k+dk)&255]
					];
				}
			}
		}
		return saturate(perlin_interp(c, u, v, w));
	}

	float turb(const vec3& p, int depth=7) const {
		float accum = 0;
		vec3 tmp = p;
		float weight = 1.0;
		for (int i=0; i<depth; ++i) {
			accum += weight * noise(tmp);
			weight *= 0.5;
			tmp *= 2;
		}
		return fabs(accum);
	}

	static vec3* m_ranvec;
	static int* m_perm_x;
	static int* m_perm_y;
	static int* m_perm_z;
};

static vec3* perlin_generate() {
	vec3* p = new vec3[256];
	for (int i=0; i<256; ++i) {
		p[i] = normalize(rndvec()*2+vec3(-1));
	}
	return p;
}

static void permute(int* p, int n) {
	for (int i=n-1; i>0; i--) {
		int target = int(drand48()*(i+1));
		std::swap(p[i], p[target]);
	}
}

static int* perlin_generate_perm() {
	int* p = new int[256];
	for (int i=0; i<256; ++i) {
		p[i] = i;
	}
	permute(p, 256);
	return p;
}

vec3* perlin::m_ranvec = perlin_generate();
int* perlin::m_perm_x = perlin_generate_perm();
int* perlin::m_perm_y = perlin_generate_perm();
int* perlin::m_perm_z = perlin_generate_perm();

