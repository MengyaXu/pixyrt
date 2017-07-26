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
#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <mutex>
//#include <Windows.h>
#include "pixyrt.h"
#include "timer.h"
#include "sphere.h"
#include "rect.h"
#include "box.h"
#include "constant_medium.h"
#include "material.h"
#include "bvh.h"
#include "pdf.h"
#include "texture.h"
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#define OUTPUT "render.bmp"
#define THREAD_NUM (8)

std::mutex guard;

bool intersect(sphere_t* ball, const ray& r, float& t) {
	vec3 oc = r.origin() - ball->center;
	float a = dot(r.direction(), r.direction());
	float b = dot(oc, r.direction());
	float c = dot(oc, oc) - pow2(ball->radius);
	float D = b*b - a*c;
	if (D > 0) {
		float sqrtD = sqrt(D);
		t = (-b - sqrtD) / a;
		if (t < 1.0 && t > 0.0) {
			return true;
		}

		t = (-b + sqrtD) / a;
		if (t < 1.0 && t > 0.0) {
			return true;
		}
	}
	return false;
}

#include "brdfs.h"

void RE_Direct(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, ReflectedLight& reflectedLight)
{
	float dotNL = saturate(dot(geometry.normal, directLight.direction));
	Vector3 irradiance = dotNL * directLight.color;

	irradiance *= PI;

	reflectedLight.directDiffuse += mulPerElem(irradiance, DiffuseBRDF(directLight, geometry, material.diffuseColor, material.specularRoughness, material.diffuseMethod));
	reflectedLight.directSpecular += mulPerElem(irradiance, SpecularBRDF(directLight, geometry, material.specularColor, material.specularRoughness));
}

void render(Image& image, Scene& scene, Workarea& warea, float metallic, float roughness, int type)
{
	{
		std::lock_guard<std::mutex> lock(guard);
		std::cout << "Rendering xy[" << warea.x() << "," << warea.y() << "] wh[" << warea.w() << "," << warea.h() << "] outxy[" << warea.outx() << "," << warea.outy() << "]" << std::endl;
	}

	int x, y;
	Vector3 color;

	int w = warea.w();
	int h = warea.h();
	Vector3 view;
	//Vector3 viewpoint = scene->viewpoint;
	Vector3 light = scene.getLight();
	sphere_t ball = scene.getBall();
	float dv, t;
	int index;

	//float metallic = 0.0f;
	//float roughness = 0.5f;
	Vector3 albedo = Vector3(1.0f);

	IncidentLight directLight;
	directLight.color = Vector3(1.0f);
	directLight.direction = normalize(light);

	GeometricContext geometry;
	Material material;
	material.diffuseColor = lerp(metallic, albedo, Vector3(0.0f));
	material.specularColor = lerp(metallic, Vector3(0.04f), albedo);
	material.specularRoughness = roughness;
	material.diffuseMethod = kDiffuseLambert;

	int sx = warea.x();
	int sy = warea.y();
	int ex = sx + warea.w();
	int ey = sy + warea.h();

	for (y = sy; y < ey; ++y) {
		for (x = sx; x < ex; ++x) {

			float u = float(x) / float(scene.getWidth());
			float v = float(y) / float(scene.getHeight());
			ray r = scene.getCamera().get_ray(u, v);
			if (intersect(&ball, r, t)) {
				Vector3 p, tv, n, L;
				tv = view * t;
				geometry.position = r.point_at_parameter(t);
				geometry.normal = normalize(geometry.position - ball.center);
				geometry.viewDir = normalize(r.origin() - geometry.position);

				float u0, v0;
				get_sphere_uv(geometry.position, u0, v0);
				if (ball.albedo) {
					albedo = gamma_to_linear(ball.albedo->value(u0, v0, geometry.position), GAMMA_FACTOR);
					material.diffuseColor = lerp(metallic, albedo, Vector3(0.0f));
					material.specularColor = lerp(metallic, Vector3(0.04f), albedo);
				}
				if (ball.roughness) {
					material.specularRoughness = gamma_to_linear(ball.roughness->value(u0, v0, geometry.position), GAMMA_FACTOR).getX() * roughness;
				}

				ReflectedLight reflectedLight;
				Vector3 color;

				switch (type)
				{
				case kOutputDefault: {
					RE_Direct(directLight, geometry, material, reflectedLight);

					if (ball.envmap) {
						float dotNL = saturate(dot(geometry.normal, directLight.direction));
						reflectedLight.indirectSpecular += ApproximateSpecularIBL(material.specularColor, material.specularRoughness, geometry.normal, geometry.viewDir, ball.envmap);
					}
				
					color = reflectedLight.getColor();
					break; }
				case kOutputDiffuse:
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directDiffuse;
					break;
				case kOutputDiffuseBurley:
					material.diffuseMethod = kDiffuseBurley;
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directDiffuse;
					break;
				case kOutputDiffuseRenormalizedBurley:
					material.diffuseMethod = kDiffuseRenormalizedBurley;
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directDiffuse;
					break;
				case kOutputDiffuseOrenNayar:
					material.diffuseMethod = kDiffuseOrenNayar;
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directDiffuse;
					break;
				case kOutputDiffuseQualitativeOrenNayar:
					material.diffuseMethod = kDiffuseQualitativeOrenNayar;
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directDiffuse;
					break;
				case kOutputDiffuseImprovedOrenNayar:
					material.diffuseMethod = kDiffuseImprovedOrenNayar;
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directDiffuse;
					break;
				case kOutputDiffuseImprovedFastOrenNayar:
					material.diffuseMethod = kDiffuseImprovedFastOrenNayar;
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directDiffuse;
					break;
				case kOutputSpecular:
					RE_Direct(directLight, geometry, material, reflectedLight);
					color = reflectedLight.directSpecular;
					break;
				case kOutputIndirectSpecular: {
					float dotNL = saturate(dot(geometry.normal, directLight.direction));
					vec3 R = reflect(-geometry.viewDir, geometry.normal);
					vec3 radiance = ball.envmap->value(R);
					color = mulPerElem(radiance, EnvBRDFApprox(material.specularColor, material.specularRoughness, dotNL));
					break; }
				case kOutputIndirectSpecularIBL: {
					float dotNL = saturate(dot(geometry.normal, directLight.direction));
					color = SpecularIBL(material.specularColor, material.specularRoughness, geometry.normal, geometry.viewDir, ball.envmap);
					break; }
				case kOutputIndirectApproximateSpecularIBL: {
					float dotNL = saturate(dot(geometry.normal, directLight.direction));
					color = ApproximateSpecularIBL(material.specularColor, material.specularRoughness, geometry.normal, geometry.viewDir, ball.envmap);
					break; }
				case kOutputEnvmap: {
					vec3 R = reflect(-geometry.viewDir, geometry.normal);
					color = gamma_to_linear(ball.envmap->value(R), GAMMA_FACTOR);
					break; }
				case kOutputAlbedo:
					color = albedo;
					break;
				case kOutputRoughness:
					color = vec3(material.specularRoughness);
					break;
				case kOutputD: {
					Vector3 N = geometry.normal;
					Vector3 V = geometry.viewDir;
					Vector3 L = directLight.direction;
					Vector3 H = normalize(L + V);
					float dotNH = saturate(dot(N, H));
					float a = roughness * roughness;
					float D = D_GGX(a, dotNH);
					color = Vector3(D);
					break; }
				case kOutputG: {
					Vector3 N = geometry.normal;
					Vector3 V = geometry.viewDir;
					Vector3 L = directLight.direction;
					float dotNL = saturate(dot(N, L));
					float dotNV = saturate(dot(N, V));
					float a = roughness * roughness;
					float G = G_Smith_Schlick_GGX(a, dotNV, dotNL);
					color = Vector3(G);
					break; }
				case kOutputF: {
					Vector3 V = geometry.viewDir;
					Vector3 L = directLight.direction;
					Vector3 H = normalize(L + V);
					Vector3 F = F_Schlick(material.specularColor, V, H);
					color = F;
					break; }
				}
				//draw_pixel(&img->buf[index], (geometry.normal + Vector3(2.0f)) * 0.5f);
				//draw_pixel(&img->buf[index], absPerElem(geometry.normal));
				//draw_pixel(&img->buf[index], Vector3(saturate(dot(geometry.normal, directLight.direction))));
				//draw_pixel(&img->buf[index], saturate(reflectedLight.getColor()));
				//draw_pixel(&img->buf[index], saturate(color));
				image.write(warea.outx() + x, warea.outy() + y, color);
			}
			else {
				//draw_pixel(&img->buf[index], scene->bgcolor);
				image.write(warea.outx() + x, warea.outy() + y, scene.getBgColor());
			}
		}
	}
}

void render(Image& image, Scene& scene, Workarea& warea, float metallic, float roughness, int type, int row)
{
#pragma omp parallel for schedule(dynamic, 1) num_threads(THREAD_NUM)
	for (int i = 0; i <= 11; i++) {
	//for (int i = 0; i < 5; i++) {
		float theta = radians(45);
		float phi = radians(i*30);
		//float phi = (i * 30 / 180.0f*PI);
		//float phi = (2 * 30 / 180.0f*PI);
		Quat rotY = Quat::rotationY(phi);
		Quat rotZ = Quat::rotationZ(-theta);

		Scene localScene = scene;
		Workarea localArea = warea;
		localScene.setLight(rotate(rotZ*rotY, Vector3::zAxis()));
		localArea.setOutPos(i*scene.getWidth(), row*scene.getHeight());
		render(image, localScene, localArea, metallic, static_cast<float>(i)/11.0f, type);
		//render(image, localScene, localArea, metallic, roughness, type);
	}
}

texture* load_texture(const char* name)
{
	int nx, ny, nn;
	unsigned char* texels = stbi_load(name, &nx, &ny, &nn, 0);
	return new image_texture(texels, nx, ny);
}

void plot_diffuse_reflectance() 
{
	std::ofstream ostrm("plot.csv");
	ostrm << "Roughness,ViewAngle,Reflectance" << std::endl;
	vec3 diffuse(1);
	vec3 N = vec3(0, 1, 0);
	vec3 L = vec3(0, 1, 0);

	//{
	//	float theta = radians(45);
	//	float phi = radians(45);
	//	Quat rotY = Quat::rotationY(phi);
	//	Quat rotZ = Quat::rotationZ(-theta);
	//	L = normalize(rotate(rotZ*rotY, Vector3::zAxis()));
	//}

	for (int j = 0; j <= 20; ++j) {
		for (int i = 0; i <= 90; i += 5) {
			float theta = radians(i);
			Quat rotX = Quat::rotationX(-theta);
			vec3 V = rotate(rotX, vec3::zAxis());
			float roughness = float(j) / 20.0f;
			//vec3 albedo = LambertDiffuse(diffuse) * PI;
			//vec3 albedo = BurleyDiffuse(diffuse, roughness, N, V, L) * PI;
			//vec3 albedo = RenormalizedBurleyDiffuse(diffuse, roughness, N, V, L) * PI;
			//vec3 albedo = OrenNayarDiffuse(diffuse, roughness, N, V, L) * PI;
			//vec3 albedo = QualitativeOrenNayarDiffuse(diffuse, roughness, N, V, L) * PI;
			//vec3 albedo = ImprovedOrenNayarDiffuse(diffuse, roughness, N, V, L) * PI;
			vec3 albedo = ImprovedFastOrenNayarDiffuse(diffuse, roughness, N, V, L) * PI;
			//ostrm << roughness << "," << i << "," << albedo.getX() << std::endl;
			ostrm << roughness << "," << 90 - i << "," << albedo.getX() << std::endl;
		}
	}
}

void plot_random_cosine_direction() 
{
	onb uvw;
	uvw.build_from_w(vec3(0, 1, 0));

	std::ofstream ostrm("plot.csv");
	ostrm << "X,Y,X" << std::endl;
	int N = 100;
	for (int i = 0; i < 100; i++) 
	{
		//vec3 v = random_cosine_direction();
		vec3 v = uvw.local(random_cosine_direction());
		ostrm << v.getX() << "," << v.getY() << "," << v.getZ() << std::endl;
	}
}

void plot()
{
	//plot_diffuse_reflectance();
	plot_random_cosine_direction();
}

#define SCENE_WIDTH  128
#define SCENE_HEIGHT 128

void render_scene()
{
	int w = SCENE_WIDTH;
	int h = SCENE_HEIGHT;

	vec3 lookfrom = vec3(0.0f, 0.0f, 45.0f);
	vec3 lookat(0);
	vec3 vup(0, 1, 0);
	float focus_dist = length(lookfrom - lookat);
	camera cam(lookfrom, lookat, vup, 45, float(w) / float(h), 0, focus_dist, 0, 1);

	sphere_t ball;
	ball.color = vec3(0, 0.5, 1.0);
	ball.center = vec3(0);
	ball.radius = 16.0f;
	ball.albedo = nullptr;
	ball.roughness = nullptr;
	ball.envmap = nullptr;
	//ball.albedo = new checker_texture(
	//	new constant_texture(vec3(0.2, 0.3, 0.1)),
	//	new constant_texture(vec3(0.9)));
	//ball.albedo = load_texture("assets/brick_diffuse.jpg");
	//ball.roughness = load_texture("assets/brick_roughness.jpg");
	ball.envmap = new cube_texture(
		load_texture("assets/px.bmp"),
		load_texture("assets/nx.bmp"),
		load_texture("assets/py.bmp"),
		load_texture("assets/ny.bmp"),
		load_texture("assets/pz.bmp"),
		load_texture("assets/nz.bmp"));

	//ball.albedo = load_texture("assets/uvchecker.bmp");
	//ball.albedo = new noise_texture(4);

	Scene scene;
	scene.setBgColor(vec3(0));
	scene.setBall(ball);
	scene.setCamera(cam);
	scene.setSize(w, h);

	Workarea warea(0, 0, w, h, 0, 0);

	Image image(w * 12, h * 4);
	//Image image(w*12, h*1);
	//Image image(w * 5, h * 7);
	//image.setFilterGamma(false);

	unsigned int t, time;
	start_timer(&t);

	float metallic = 1.0f;
	float roughness = 0.0f;
	//render(&img, &scene, &warea, metallic, roughness, kOutputD, 0);
	//render(&img, &scene, &warea, metallic, roughness, kOutputG, 1);
	//render(&img, &scene, &warea, metallic, roughness, kOutputF, 2);
	//render(&img, &scene, &warea, metallic, roughness, kOutputDiffuse, 3);
	//render(&img, &scene, &warea, metallic, roughness, kOutputSpecular, 4);

	//render(image, scene, warea, metallic, roughness, kOutputDiffuse, 0);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuseBurley, 1);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuseRenormalizedBurley, 2);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuseOrenNayar, 3);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuseQualitativeOrenNayar, 4);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuseImprovedOrenNayar, 5);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuseImprovedFastOrenNayar, 6);

	render(image, scene, warea, metallic, roughness, kOutputIndirectSpecularIBL, 0);
	render(image, scene, warea, metallic, roughness, kOutputIndirectApproximateSpecularIBL, 1);
	render(image, scene, warea, metallic, roughness, kOutputSpecular, 2);
	render(image, scene, warea, metallic, roughness, kOutputDefault, 3);

	//render(image, scene, warea, metallic, roughness, kOutputAlbedo, 0);
	//render(image, scene, warea, metallic, roughness, kOutputRoughness, 1);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuse, 2);
	//render(image, scene, warea, metallic, roughness, kOutputSpecular, 3);
	//render(image, scene, warea, metallic, roughness, kOutputIndirectSpecular, 4);
	//render(image, scene, warea, metallic, roughness, kOutputDefault, 5);

	time = stop_timer(&t);
	print_timer(time);

	stbi_write_bmp(OUTPUT, image.width(), image.height(), sizeof(Image::rgb), image.pixels());
}

#include "scene.h"
#include "shade.h"
#define MAX_DEPTH 50
void raytrace()
{

	//int nx = 200;
	//int ny = 100;
	//int ns = 100;
	//int nx = 1200;
	//int ny = 800;
	int nx = 500;
	int ny = 500;
	int ns = 300;

	std::unique_ptr<Image> image(new Image(nx,ny));
	//image->setFilterGamma(false);
	//std::unique_ptr<Image::rgb[]> image(new Image::rgb[nx*ny]);

	//hitable* list[5];
	//list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0.8,0.3,0.3)));
	//list[1] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.8,0.8,0.0)));
	//list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 1.0));
	//list[3] = new sphere(vec3(-1, 0, -1), 0.5, new metal(vec3(0.8, 0.8, 0.8), 0.3));
	//hitable* world = new hitable_list(list, 4);

	//list[0] = new sphere(vec3(0, 0, -1), 0.5, new lambertian(vec3(0.1, 0.2, 0.6)));
	//list[1] = new sphere(vec3(0, -100.5, -1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
	//list[2] = new sphere(vec3(1, 0, -1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.5));
	//list[3] = new sphere(vec3(-1, 0, -1), 0.5, new dielectric(1.5));
	//list[4] = new sphere(vec3(-1, 0, -1), -0.45, new dielectric(1.5));
	//hitable* world = new hitable_list(list, 5);

	//float R = cos(M_PI/4);
	//list[0] = new sphere(vec3(-R, 0, -1), R, new lambertian(vec3(0,0,1)));
	//list[1] = new sphere(vec3( R, 0, -1), R, new lambertian(vec3(1,0,0)));
	//hitable* world = new hitable_list(list, 2);

	//hitable* world = random_scene();
	//hitable* world = two_spheres();
	//hitable* world = two_perlin_spheres();
	//hitable* world = earth();
	//hitable* world = simple_light();
	//hitable* world = cornel_box();
	//hitable* world = cornel_smoke();
	//hitable* world = final_scene();
	hitable* world = refract_sphere();

	//hitable* list[2];
	//list[0] = new xz_rect(213, 343, 227, 332, 554, 0);
	//list[1] = new sphere(vec3(190, 90, 190), 90, 0);
	//hitable* lights = new hitable_list(list, 2);
	hitable* lights = nullptr;

	//vec3 lookfrom = vec3(3,3,2);
	//vec3 lookat = vec3(0,0,0);
	//vec3 lookfrom = vec3(13,2,3) * 2;
	//vec3 lookat = vec3(0,0,0);
	vec3 lookfrom = vec3(278, 278, -800);
	vec3 lookat = vec3(278, 278, 0);
	vec3 vup = vec3(0, 1, 0);
	float vfov = 40;
	//float dist_to_focus = length(lookfrom - lookat);
	float dist_to_focus = 10;
	float aperture = 0.0;
	camera cam(lookfrom, lookat, vup, vfov, float(nx) / float(ny), aperture, dist_to_focus, 0.0, 1.0);

#pragma omp parallel for schedule(dynamic, 1) num_threads(8)
	for (int j = 0; j<ny; ++j) {
		std::cerr << "Rendering (y = " << j << ") " << (100.0 * j / (ny - 1)) << "%" << std::endl;
		for (int i = 0; i<nx; ++i) {
			vec3 col(0);
			for (int s = 0; s<ns; ++s) {
				float u = float(i + drand48()) / float(nx);
				float v = float(j + drand48()) / float(ny);
				col += shade(cam.get_ray(u, v), world, lights, 0, MAX_DEPTH);
				//col += shade(cam.get_ray(u, v), world, lights, 0, 1);
			}

			col /= float(ns);
			//image->write(i, ny-j-1, sqrtPerElem(col));
			image->write(i, ny - j - 1, col);
		}
	}

	stbi_write_bmp("render.bmp", nx, ny, sizeof(Image::rgb), image->pixels());
}

int main(int argc, char* argv[]) 
{
	//plot();
	//render_scene();
	raytrace();
	return 0;
}