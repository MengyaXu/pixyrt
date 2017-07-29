#pragma once

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

	std::unique_ptr<Image> image(new Image(nx, ny));
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
			image->write(i, ny - j - 1, col);
		}
	}

	stbi_write_bmp("render.bmp", nx, ny, sizeof(Image::rgb), image->pixels());
}
