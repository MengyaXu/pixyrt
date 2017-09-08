#pragma once

#include "brdfs.h"

void RE_Direct(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, ReflectedLight& reflectedLight)
{
	float dotNL = saturate(dot(geometry.normal, directLight.direction));
	Vector3 irradiance = dotNL * directLight.color;

	irradiance *= PI;

	reflectedLight.directDiffuse += mulPerElem(irradiance, DiffuseBRDF(directLight, geometry, material.diffuseColor, material.specularRoughness, material.diffuseMethod));
	reflectedLight.directSpecular += mulPerElem(irradiance, SpecularBRDF(directLight, geometry, material.specularColor, material.specularRoughness));
}

enum {
	kOutputDefault = 0,
	kOutputD,
	kOutputG,
	kOutputF,
	kOutputDiffuse,
	kOutputDiffuseBurley,
	kOutputDiffuseRenormalizedBurley,
	kOutputDiffuseOrenNayar,
	kOutputDiffuseQualitativeOrenNayar,
	kOutputDiffuseImprovedOrenNayar,
	kOutputDiffuseImprovedFastOrenNayar,
	kOutputDiffuseGotandaOrenNayar,
	kOutputDiffuseGotandaNormalizedOrenNayar,
	kOutputDiffuseGGXApprox,
	kOutputSpecular,
	kOutputIndirectSpecular,
	kOutputIndirectSpecularIBL,
	kOutputIndirectApproximateSpecularIBL,
	kOutputEnvMap,
	kOutputAlbedo,
	kOutputRoughness
};

typedef void (*RenderFunction)(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight);

void RenderDefault(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	RE_Direct(directLight, geometry, material, reflectedLight);

	// ambient
	//reflectedLight.directDiffuse += material.diffuseColor * 0.2f;

	if (samplers.envmap) {
		float dotNL = saturate(dot(geometry.normal, directLight.direction));
		reflectedLight.indirectSpecular += ApproximateSpecularIBL(material.specularColor, material.specularRoughness, geometry.normal, geometry.viewDir, samplers.envmap);
	}
}

void RenderCookTorranceD(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Vector3 N = geometry.normal;
	Vector3 V = geometry.viewDir;
	Vector3 L = directLight.direction;
	Vector3 H = normalize(L + V);
	float dotNH = saturate(dot(N, H));
	float a = pow2(material.specularRoughness);
	float D = D_GGX(a, dotNH);
	reflectedLight.directDiffuse = vec3(D);
}

void RenderCookTorranceG(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Vector3 N = geometry.normal;
	Vector3 V = geometry.viewDir;
	Vector3 L = directLight.direction;
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	float a = pow2(material.specularRoughness);
	float G = G_Smith_Schlick_GGX(a, dotNV, dotNL);
	reflectedLight.directDiffuse = vec3(G);
}

void RenderCookTorranceF(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Vector3 V = geometry.viewDir;
	Vector3 L = directLight.direction;
	Vector3 H = normalize(L + V);
	Vector3 F = F_Schlick(material.specularColor, V, H);
	reflectedLight.directDiffuse = F;
}

void RenderDiffuseLambert(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	RE_Direct(directLight, geometry, material, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseBurley(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseBurley;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseRenormalizedBurley(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseRenormalizedBurley;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseOrenNayar(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseOrenNayar;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseQualitativeOrenNayar(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseQualitativeOrenNayar;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseImprovedOrenNayar(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseImprovedOrenNayar;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseImprovedFastOrenNayar(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseImprovedFastOrenNayar;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseGotandaOrenNayar(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseGotandaOrenNayar;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseGotandaNormalizedOrenNayar(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseGotandaNormalizedOrenNayar;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderDiffuseGGXApprox(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	Material mat = material;
	mat.diffuseMethod = kDiffuseGGXApprox;
	RE_Direct(directLight, geometry, mat, reflectedLight);
	reflectedLight.directSpecular = vec3(0);
}

void RenderSpecular(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	RE_Direct(directLight, geometry, material, reflectedLight);
	reflectedLight.directDiffuse = vec3(0);
}

void RenderIndirectSpecular(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	if (samplers.envmap) {
		float dotNL = saturate(dot(geometry.normal, directLight.direction));
		vec3 R = reflect(-geometry.viewDir, geometry.normal);
		vec3 radiance = samplers.envmap->value(R);
		reflectedLight.indirectSpecular = mulPerElem(radiance, EnvBRDFApprox(material.specularColor, material.specularRoughness, dotNL));
	}
}

void RenderIndirectSpecularIBL(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	if (samplers.envmap) {
		float dotNL = saturate(dot(geometry.normal, directLight.direction));
		reflectedLight.indirectSpecular = SpecularIBL(material.specularColor, material.specularRoughness, geometry.normal, geometry.viewDir, samplers.envmap);
	}
}

void RenderIndirectApproximateSpecularIBL(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	if (samplers.envmap) {
		float dotNL = saturate(dot(geometry.normal, directLight.direction));
		reflectedLight.indirectSpecular = ApproximateSpecularIBL(material.specularColor, material.specularRoughness, geometry.normal, geometry.viewDir, samplers.envmap);
	}
}

void RenderEnvMap(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	if (samplers.envmap) {
		vec3 R = reflect(-geometry.viewDir, geometry.normal);
		reflectedLight.directDiffuse = gamma_to_linear(samplers.envmap->value(R), GAMMA_FACTOR);
	}
}

void RenderAlbedo(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	reflectedLight.directDiffuse = material.albedo;
}

void RenderRoughness(const IncidentLight& directLight, const GeometricContext& geometry, const Material& material, const Samplers& samplers, ReflectedLight& reflectedLight)
{
	reflectedLight.directDiffuse = vec3(material.specularRoughness);
}

RenderFunction renderFunctions[] = {
	RenderDefault,
	RenderCookTorranceD,
	RenderCookTorranceG,
	RenderCookTorranceF,
	RenderDiffuseLambert,
	RenderDiffuseBurley,
	RenderDiffuseRenormalizedBurley,
	RenderDiffuseOrenNayar,
	RenderDiffuseQualitativeOrenNayar,
	RenderDiffuseImprovedOrenNayar,
	RenderDiffuseImprovedFastOrenNayar,
	RenderDiffuseGotandaOrenNayar,
	RenderDiffuseGotandaNormalizedOrenNayar,
	RenderDiffuseGGXApprox,
	RenderSpecular,
	RenderIndirectSpecular,
	RenderIndirectSpecularIBL,
	RenderIndirectApproximateSpecularIBL,
	RenderEnvMap,
	RenderAlbedo,
	RenderRoughness
};

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
	Vector3 light = scene.getLight();
	sphere_t ball = scene.getBall();
	float t;

	IncidentLight directLight;
	directLight.color = Vector3(1.0f);
	directLight.direction = normalize(light);

	GeometricContext geometry;
	Material material;
	material.albedo = vec3(1);
	material.diffuseColor = lerp(metallic, material.albedo, Vector3(0.0f));
	material.specularColor = lerp(metallic, Vector3(0.04f), material.albedo);
	material.specularRoughness = roughness;
	material.diffuseMethod = kDiffuseLambert;

	Samplers samplers;
	samplers.albedo = nullptr;
	samplers.roughness = nullptr;
	samplers.envmap = nullptr;

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
				geometry.position = r.point_at_parameter(t);
				geometry.normal = normalize(geometry.position - ball.center);
				geometry.viewDir = normalize(r.origin() - geometry.position);

				samplers.albedo = ball.albedo;
				samplers.roughness = ball.roughness;
				samplers.envmap = ball.envmap;

				get_sphere_uv(geometry.normal, geometry.uv0.x, geometry.uv0.y);
				if (samplers.albedo) {
					material.albedo = gamma_to_linear(samplers.albedo->value(geometry.uv0.x, geometry.uv0.y, geometry.position), GAMMA_FACTOR);
					material.diffuseColor = lerp(metallic, material.albedo, Vector3(0.0f));
					material.specularColor = lerp(metallic, Vector3(0.04f), material.albedo);
				}
				if (samplers.roughness) {
					material.specularRoughness = gamma_to_linear(samplers.roughness->value(geometry.uv0.x, geometry.uv0.y, geometry.position), GAMMA_FACTOR).getX() * roughness;
				}

				ReflectedLight reflectedLight;
				renderFunctions[type](directLight, geometry, material, samplers, reflectedLight);
				Vector3 color = reflectedLight.getColor();
				//draw_pixel(&img->buf[index], (geometry.normal + Vector3(2.0f)) * 0.5f);
				//draw_pixel(&img->buf[index], absPerElem(geometry.normal));
				//draw_pixel(&img->buf[index], Vector3(saturate(dot(geometry.normal, directLight.direction))));
				//draw_pixel(&img->buf[index], saturate(reflectedLight.getColor()));
				//draw_pixel(&img->buf[index], saturate(color));
				image.write(warea.outx() + x, warea.outy() + y, color);
			}
			else {
				image.write(warea.outx() + x, warea.outy() + y, scene.getBgColor());

				// sky light
				//vec3 n = normalize(r.direction());
				//float t = 0.5f * (n.getY() + 1.0f);
				//vec3 color = (1.0f - t) * vec3(1.0f) + t*vec3(0.5f, 0.7f, 1.0f);
				//image.write(warea.outx() + x, warea.outy() + y, color);
			}
		}
	}
}

void render(Image& image, Scene& scene, Workarea& warea, float metallic, float roughness, int type, int row)
{
#pragma omp parallel for schedule(dynamic, 1) num_threads(THREAD_NUM)
	//for (int i = 0; i <= 11; i++) {
		for (int i = 0; i < 5; i++) {
		float theta = radians(45);
		float phi = radians(i*30);
		//float phi = radians(i * 30);
		//float phi = radians(2 * 30);
		Quat rotY = Quat::rotationY(phi);
		Quat rotZ = Quat::rotationZ(-theta);

		Scene localScene = scene;
		Workarea localArea = warea;
		localScene.setLight(rotate(rotZ*rotY, Vector3::zAxis()));
		localArea.setOutPos(i*scene.getWidth(), row*scene.getHeight());
		//render(image, localScene, localArea, metallic, static_cast<float>(i) / 11.0f, type);
		render(image, localScene, localArea, metallic, roughness, type);
	}
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
	//ball.envmap = new cube_texture(
	//	load_texture("assets/px.bmp"),
	//	load_texture("assets/nx.bmp"),
	//	load_texture("assets/py.bmp"),
	//	load_texture("assets/ny.bmp"),
	//	load_texture("assets/pz.bmp"),
	//	load_texture("assets/nz.bmp"));

	//ball.albedo = load_texture("assets/uvchecker.bmp");
	//ball.albedo = new noise_texture(4);

	Scene scene;
	scene.setBgColor(vec3(0));
	scene.setBall(ball);
	scene.setCamera(cam);
	scene.setSize(w, h);

	Workarea warea(0, 0, w, h, 0, 0);

	//Image image(w * 12, h * 1);
	//Image image(w*12, h*1);
	Image image(w * 5, h * 10);
	//image.setFilterGamma(false);
	//image.clear(vec3(1.0f));

	unsigned int t, time;
	start_timer(&t);

	float metallic = 0.0f;
	float roughness = 0.0f;
	//render(image, scene, warea, metallic, roughness, kOutputD, 0);
	//render(image, scene, warea, metallic, roughness, kOutputG, 0);
	//render(image, scene, warea, metallic, roughness, kOutputF, 0);

	//render(image, scene, warea, metallic, roughness, kOutputD, 0);
	//render(image, scene, warea, metallic, roughness, kOutputG, 1);
	//render(image, scene, warea, metallic, roughness, kOutputF, 2);
	//render(image, scene, warea, metallic, roughness, kOutputDiffuse, 3);
	//render(image, scene, warea, metallic, roughness, kOutputSpecular, 4);
	//render(image, scene, warea, metallic, roughness, kOutputDefault, 5);

	render(image, scene, warea, metallic, roughness, kOutputDiffuse, 0);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseBurley, 1);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseRenormalizedBurley, 2);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseOrenNayar, 3);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseQualitativeOrenNayar, 4);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseImprovedOrenNayar, 5);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseImprovedFastOrenNayar, 6);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseGotandaOrenNayar, 7);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseGotandaNormalizedOrenNayar, 8);
	render(image, scene, warea, metallic, roughness, kOutputDiffuseGGXApprox, 9);

	//render(image, scene, warea, metallic, roughness, kOutputIndirectSpecularIBL, 0);
	//render(image, scene, warea, metallic, roughness, kOutputIndirectApproximateSpecularIBL, 1);
	//render(image, scene, warea, metallic, roughness, kOutputSpecular, 2);
	//render(image, scene, warea, metallic, roughness, kOutputDefault, 3);

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