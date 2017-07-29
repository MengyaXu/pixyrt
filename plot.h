#pragma once

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