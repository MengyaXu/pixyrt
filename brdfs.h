#pragma once

// BRDFs

enum {
	kDiffuseLambert,
	kDiffuseBurley,
	kDiffuseRenormalizedBurley,
	kDiffuseOrenNayar,
	kDiffuseQualitativeOrenNayar,
	kDiffuseImprovedOrenNayar,
	kDiffuseImprovedFastOrenNayar,
	kDiffuseGotandaOrenNayar,
	kDiffuseGotandaNormalizedOrenNayar,
	kDiffuseGGXApprox,
};

Vector3 LambertDiffuse(const Vector3& diffuseColor)
{
	return diffuseColor / PI;
}

Vector3 BurleyDiffuse(const Vector3& diffuseColor, float a, const Vector3& N, const Vector3& V, const Vector3& L)
{
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	Vector3 H = normalize(L + V);
	float dotLH = saturate(dot(L, H));
	float fd90 = 0.5f + 2.0f * dotLH * dotLH * a;
	float nl = schlick(dotNL, 1.0f, fd90);
	float nv = schlick(dotNV, 1.0f, fd90);
	return diffuseColor * ((nl*nv) / PI);
}

// Frostbite
Vector3 RenormalizedBurleyDiffuse(const Vector3& diffuseColor, float roughness, const Vector3& N, const Vector3& V, const Vector3& L)
{
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	Vector3 H = normalize(L + V);
	float dotLH = saturate(dot(L, H));
	float energyBias = mix(0.0f, 0.5f, roughness);
	float energyFactor = mix(1.0f, 1.0f / 1.51f, roughness);
	float fd90 = energyBias + 2.0f * dotLH * dotLH * roughness;
	float nl = schlick(dotNL, 1.0f, fd90);
	float nv = schlick(dotNV, 1.0f, fd90);
	return diffuseColor * ((nl*nv*energyFactor) / PI);
}

// http://ruh.li/GraphicsOrenNayar.html
Vector3 OrenNayarDiffuse(const Vector3& diffuseColor, float a, const Vector3& N, const Vector3& V, const Vector3& L)
{
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	Vector3 H = normalize(L + V);
	float dotLH = saturate(dot(L, H));

#if 0
	float angleVN = acosf(dotNV);
	float angleLN = acosf(dotNL);

	float alpha = max(angleVN, angleLN);
	float beta = min(angleVN, angleLN);
	float gamma = dot(normalize(V - N * dotNV), normalize(L - N * dotNL));

	float roughnessSquared = a;
	float roughnessSquared9 = (roughnessSquared / (roughnessSquared + 0.09f));

	// calculate C1, C2 and C3
	float C1 = 1.0f - 0.5f * (roughnessSquared / (roughnessSquared + 0.33f));
	float C2 = 0.45f * roughnessSquared9;

	if (gamma >= 0.0f) {
		C2 *= sinf(alpha);
	}
	else {
		//C2 *= (sinf(alpha) - pow3((2.0f*beta) / PI));
		C2 *= (sinf(alpha) - powf((2.0f * beta) / PI, 3.0f));
	}

	float powValue = (4.0f * alpha * beta) / (PI*PI);
	float C3 = 0.125f * roughnessSquared9 * powValue * powValue;

	// now calculate both main parts of the formula
	float A = gamma * C2 * tanf(beta);
	float B = (1.0f - abs(gamma)) * C3 * tanf((alpha + beta) / 2.0f);

	// put it all together
	//float L1 = max(0.0f, dotNL) * (C1 + A + B);
	float L1 = (C1 + A + B);

	// also calculate interreflection
	float twoBetaPi = 2.0f * beta / PI;
	//float L2 = 0.17f * max(0.0f, dotNL) * (roughnessSquared / (roughnessSquared + 0.13f)) * (1.0f - gamma * twoBetaPi * twoBetaPi);
	float L2 = 0.17f * (roughnessSquared / (roughnessSquared + 0.13f)) * (1.0f - gamma * twoBetaPi * twoBetaPi);

	Vector3 rho_div_pi = diffuseColor / PI;
	return mulPerElem(diffuseColor, rho_div_pi * (L1 + L2));
	//return diffuseColor * (L1 + L2);
#endif
	// BRDF Viewer
#if 1
	Vector3 v1 = V - N * dotNV;
	Vector3 v2 = L - N * dotNL;
	float theta_r = acosf(dotNV);
	float sigma2 = a;
	float cos_phi_diff = dot(normalize(V - N * dotNV), normalize(L - N * dotNL));
	if (isnan(cos_phi_diff)) cos_phi_diff = 1.0f;
	//if (length(v1) < EPSILON || length(v2) < EPSILON) std::cout << cos_phi_diff << std::endl;
	//if (length(v1) < EPSILON || length(v2) < EPSILON) cos_phi_diff = 1.0f;
	float theta_i = acosf(dotNL);
	float alpha = max(theta_i, theta_r);
	float beta = min(theta_i, theta_r);
	if (alpha > PI / 2) return Vector3(0.0f);

	float C1 = 1.0f - 0.5f * sigma2 / (sigma2 + 0.33f);
	float C2 = 0.45f * sigma2 / (sigma2 + 0.09f);
	if (cos_phi_diff >= 0.0f) C2 *= sinf(alpha);
	else C2 *= (sinf(alpha) - pow(2 * beta / PI, 3.0f));
	float C3 = 0.125f * sigma2 / (sigma2 + 0.09f) * pow((4 * alpha * beta) / (PI*PI), 2.0f);
	float L1 = (C1 + cos_phi_diff * C2 * tanf(beta) + (1.0f - abs(cos_phi_diff)) * C3 * tanf((alpha + beta) / 2));
	float L2 = 0.17f * (sigma2 / (sigma2 + 0.13f)) * (1.0f - cos_phi_diff * (4.0f * beta * beta) / (PI*PI));
	Vector3 rho_div_pi = diffuseColor / PI;
	Vector3 rho_rho_div_pi = mulPerElem(rho_div_pi, diffuseColor);
	return rho_div_pi * L1 + rho_rho_div_pi * L2;
#endif
}

// http://ruh.li/GraphicsOrenNayar.html
Vector3 QualitativeOrenNayarDiffuse(const Vector3& diffuseColor, float a, const Vector3& N, const Vector3& V, const Vector3& L)
{
	// calculate intermediary values
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	float dotLV = saturate(dot(L, V));
	Vector3 H = normalize(L + V);
	float dotLH = saturate(dot(L, H));

	float angleVN = acosf(dotNV);
	float angleLN = acosf(dotNL);

	float alpha = max(angleVN, angleLN);
	float beta = min(angleVN, angleLN);
	float gamma = dot(normalize(V - N * dotNV), normalize(L - N * dotNL));
	if (isnan(gamma)) gamma = 1.0f;

	float roughnessSquared = a;

	// calculate A and B
	float A = 1.0f - 0.5f * (roughnessSquared / (roughnessSquared + 0.57f));
	float B = 0.45f * (roughnessSquared / (roughnessSquared + 0.09f));
	float C = sinf(alpha) * tanf(beta);

	// put it all together
	//float L1 = max(0.0f, dotNL) * (A + B * max(0.0f, gamma) * C);
	float L1 = (A + B * max(0.0f, gamma) * C);

	Vector3 rho_div_pi = diffuseColor / PI;
	return rho_div_pi * L1;
}

// http://mimosa-pudica.net/improved-oren-nayar.html
Vector3 ImprovedOrenNayarDiffuse(const Vector3& diffuseColor, float a, const Vector3& N, const Vector3& V, const Vector3& L)
{
	// calculate intermediary values
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	float dotLV = saturate(dot(L, V));
	Vector3 H = normalize(L + V);
	float dotLH = saturate(dot(L, H));

	float s = dotLV - dotNL * dotNV;
	float t = mix(1.0f, max(dotNL, dotNV), step(0.0f, s));
	float st = s * (1.0f / (t + EPSILON));

	float sigma2 = a;
	Vector3 A = diffuseColor * (0.17f * sigma2 / (sigma2 + 0.13f)) + Vector3(1.0f - 0.5f * sigma2 / (sigma2 + 0.33f));
	float B = 0.45f * sigma2 / (sigma2 + 0.09f);
	//return mulPerElem(diffuseColor * max(0.0f, dotNL), (A + Vector3(B * s / t) / PI));
	return mulPerElem(diffuseColor, (A * RECIPROCAL_PI + Vector3(B * st * RECIPROCAL_PI)));
}
// http://mimosa-pudica.net/improved-oren-nayar.html
Vector3 ImprovedFastOrenNayarDiffuse(const Vector3& diffuseColor, float a, const Vector3& N, const Vector3& V, const Vector3& L)
{
	// calculate intermediary values
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	float dotLV = saturate(dot(L, V));
	Vector3 H = normalize(L + V);
	float dotLH = saturate(dot(L, H));

	float s = dotLV - dotNL * dotNV;
	float t = mix(1.0f, max(dotNL, dotNV), step(0.0f, s));
	float st = s * (1.0f / (t + EPSILON));

	float A = 1.0f / ((PI * 0.5f - 2.0f / 3.0f) * a + PI);
	float B = a * A;
	return diffuseColor * (A + B * st);
}
// GOTANDA
// http://brabl.com/diffuse-approximation-for-consoles/
Vector3 GotandaOrenNayarDiffuse(const Vector3& diffuseColor, float roughness, const Vector3& N, const Vector3& V, const Vector3& L)
{
	// calculate intermediary values
	float f0 = 0.04f;
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	float dotLV = saturate(dot(L, V));
	Vector3 H = normalize(L + V);
	float dotVH = saturate(dot(V, H));
	float a = roughness*roughness;
	float a2 = a*a;

	// float fTheta = f0 + (1.0 - f0) * pow(2.0, -9.60232*pow(VdotH, 8) - 8.58092 * VdotH);
	float Fr1 = (1.0 - (0.542026*a2 + 0.303573*a) / (a2 + 1.36053));
	float Fr2 = (1.0 - (pow(1.0 - dotNV, 5.0 - 4.0*a2)) / (a2 + 1.36053));
	float Fr3 = (-0.733996*a2*a + 1.50912*a2 - 1.16402*a);
	float Fr4 = (pow(1.0 - dotNV, 1.0 + (1.0 / (39.0*a2*a2 + 1.0))));
	float Fr = Fr1*Fr2*(Fr3*Fr4 + 1);
	float Lm1 = (max(1.0 - (2.0*a), 0.0)*(1.0 - pow(1.0 - dotNL, 5.0)) + min(2.0*a, 1.0));
	float Lm2 = ((1.0 - 0.5*a)*dotNL + (0.5*a)*(pow(dotNL, 2.0)));
	float Lm = Lm1 * Lm2;
	float Vd1 = (a2 / ((a2 + 0.09)*(1.31072 + 0.995584*dotNV)));
	float Vd2 = (1.0 - (pow(1.0 - dotNL, (1.0 - 0.3726732*(dotNV*dotNV)) / (0.188566 + 0.38841*dotNV))));
	float Vd = Vd1*Vd2;
	float Bp = dotLV - (dotNV*dotNL);
	if (Bp < 0.0) Bp *= 1.4*dotNL*dotNL;
	float L1 = 1.05f*(1.0 - f0)*(Fr*Lm + Vd*Bp);
	return diffuseColor * (L1*RECIPROCAL_PI);

	//float fTheta = f0 + (1.0f - f0) * powf(2.0f, -9.60232*powf(dotVH, 8.0f) - 8.58092 * dotVH);
	//float Bp = 0.0f;
	//float Bpl = dotLV - (dotNV*dotNL);
	//float Fr = (1.0f - (0.542026f*a + 0.303573f*roughness)/(a+1.36053f)) *
	//	(1.0f - (powf(1.0f-dotNV, 5.0f-4.0f*a))/(a+1.36053f)) *
	//	((-0.733996f*a*roughness+1.50912f*a-1.16402f*roughness)*
	//	(powf(1.0f-dotNV, 1.0f+(1.0f/(39.0f*a*a+1.0f))))+1.0f);
	//float Lm = (max(1.0f-(2.0f*roughness),0.0f)*(1.0f-(1.0f-pow(dotNL,5.0f)))+min(2.0f*roughness,1.0f))*((1.0f-0.5f*a)*dotNL+(0.5*a)*(powf(dotNL,2.0f)));
	//float Vd = (a/((a+0.09f)*(1.31072f+0.995584f*dotNV)))*(1.0f-(powf(1.0f-dotNL,(1.0f-0.3726732*(dotNV*dotNV))/(0.188566f+0.38841*dotNV))));
	//if (Bpl < 0.0f)
	//	Bp = 1.4f*dotNV*dotNL*(dotLV-(dotNV*dotNL));
	//else 
	//	Bp = dotLV-((dotNV*dotNL));
	//float aDiffuse = (21.0f/(20.0f*PI))*(1.0f-f0)*((Fr*Lm)+(Vd*Bp));
	//return diffuseColor * aDiffuse;
}
// GOTANDA
// http://brabl.com/diffuse-approximation-for-consoles/
Vector3 GotandaNormalizedOrenNayarDiffuse(const Vector3& diffuseColor, float roughness, const Vector3& N, const Vector3& V, const Vector3& L)
{
	float f0 = 0.04f;
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	float dotLV = saturate(dot(L, V));
	Vector3 H = normalize(L + V);
	float dotVH = saturate(dot(V, H));
	float a = roughness*roughness;

	//float fTheta = f0 + (1.0f - f0) * powf(2.0f, -9.60232*powf(dotVH, 8.0f) - 8.58092 * dotVH);

	float fdiff = 1.05f * (1.0f - f0) * (1.0f - powf(1.0f - dotNL, 5.0f)) * (1.0f - powf(1.0f - dotNV, 5.0f));
	float A = 1.0f - 0.5f * (a / (a + 0.65f));
	float B = 0.45f * (a / (a + 0.09f));
	float Bp = dotLV - dotNV*dotNL;
	float Bm = min(1.0f, dotNL / (dotNV + EPSILON));
	float L1 = (1.0f - f0) * (fdiff*dotNL*A + B*Bp*Bm);
	return diffuseColor * (L1 * RECIPROCAL_PI);
}

// GGX Approximation
Vector3 GGXApproxDiffuse(const Vector3& diffuseColor, float roughness, const Vector3& N, const Vector3& V, const Vector3& L)
{
	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	float dotLV = saturate(dot(L, V));
	Vector3 H = normalize(L + V);
	float dotVH = saturate(dot(V, H));
	float dotNH = saturate(dot(N, H));
	float a = roughness*roughness;

	float facing = 0.5f + 0.5f * dotLV;
	float rough = facing * (0.9f - 0.4f*facing) * ((0.5f + dotNH)/dotNH);
	float Smooth = 1.05f * (1.0f - powf(1.0f - dotNL, 5.0f)) * (1.0f - powf(1.0f-dotNL, 5.0f));
	float single = mix(Smooth, rough, a) * RECIPROCAL_PI;
	float multi = 0.1159f * a;
	return mulPerElem(diffuseColor, vec3(single) + diffuseColor * multi);
}

Vector3 DiffuseBRDF(const IncidentLight& directLight, const GeometricContext& geometry, const Vector3& diffuseColor, float roughness, int method)
{
	if (method == kDiffuseLambert) {
		return LambertDiffuse(diffuseColor);
	}

	float a = roughness * roughness;
	switch (method) {
	case kDiffuseBurley:
		return BurleyDiffuse(diffuseColor, a, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseRenormalizedBurley:
		return RenormalizedBurleyDiffuse(diffuseColor, roughness, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseOrenNayar:
		return OrenNayarDiffuse(diffuseColor, a, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseQualitativeOrenNayar:
		return QualitativeOrenNayarDiffuse(diffuseColor, a, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseImprovedOrenNayar:
		return ImprovedOrenNayarDiffuse(diffuseColor, a, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseImprovedFastOrenNayar:
		return ImprovedFastOrenNayarDiffuse(diffuseColor, a, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseGotandaOrenNayar:
		return GotandaOrenNayarDiffuse(diffuseColor, roughness, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseGotandaNormalizedOrenNayar:
		return GotandaNormalizedOrenNayarDiffuse(diffuseColor, roughness, geometry.normal, geometry.viewDir, directLight.direction);
	case kDiffuseGGXApprox:
		return GGXApproxDiffuse(diffuseColor, roughness, geometry.normal, geometry.viewDir, directLight.direction);
	}

	return LambertDiffuse(diffuseColor);
}

Vector3 F_Schlick(const Vector3& specularColor, const Vector3& H, const Vector3& V)
{
	return (specularColor + (Vector3(1.0) - specularColor) * pow(1.0f - saturate(dot(V, H)), 5.0f));
}

float D_GGX(float a, float dotNH)
{
	float a2 = a*a;
	float dotNH2 = dotNH*dotNH;
	float d = dotNH2 * (a2 - 1.0f) + 1.0f;
	return a2 / (PI * d * d);
}

float G_Smith_Schlick_GGX(float a, float dotNV, float dotNL)
{
	float k = a*a*0.5 + EPSILON;
	//float k = pow2(sqrtf(a)+1) / 8.0f;
	float gl = dotNL / (dotNL * (1.0f - k) + k);
	float gv = dotNV / (dotNV * (1.0f - k) + k);
	return gl*gv;
}

Vector3 SpecularBRDF(const IncidentLight& directLight, const GeometricContext& geometry, const Vector3& specularColor, float roughness)
{
	Vector3 N = geometry.normal;
	Vector3 V = geometry.viewDir;
	Vector3 L = directLight.direction;

	float dotNL = saturate(dot(N, L));
	float dotNV = saturate(dot(N, V));
	Vector3 H = normalize(L + V);
	float dotNH = saturate(dot(N, H));
	float dotVH = saturate(dot(V, H));
	float dotLV = saturate(dot(L, V));
	float a = roughness * roughness;

	float D = D_GGX(a, dotNH);
	float G = G_Smith_Schlick_GGX(a, dotNV, dotNL);
	Vector3 F = F_Schlick(specularColor, V, H);
	return (F*(G*D)) / (4.0f * dotNL * dotNV + EPSILON);
}

vec3 EnvBRDFApprox(vec3 specularColor, float roughness, float NoV) {
	const vec4 c0 = vec4(-1, -0.0275, -0.572, 0.022);
	const vec4 c1 = vec4(1, 0.0425, 1.04, -0.04);
	vec4 r = roughness * c0 + c1;
	float a004 = min(r.getX() * r.getX(), exp2(-9.28 * NoV)) * r.getX() + r.getY();
	float A = -1.04 * a004 + r.getZ();
	float B =  1.04 * a004 + r.getW();
	//vec2 AB = vec2(-1.04, 1.04) * a004 + r.zw;
	//return specularColor * AB.x + AB.y;
	return specularColor * A + vec3(B);
}

// http://brabl.com/ibl/
float radicalInverse_VdC(unsigned int bits)
{
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);

	return float(bits) * 2.3283064365386963e-10; // / 0x100000000
}

vec2 Hammersley(unsigned int i, unsigned int n)
{
	return vec2(float(i) / float(n), radicalInverse_VdC(i));
}

vec3 ImportanceSampleGGX(vec2 Xi, float Roughness, vec3 N)
{
	float a = Roughness * Roughness;
	//float a = pow(Roughness + 1.0f, 2.0f);
	float Phi = 2.0f * PI * Xi.x;
	float CosTheta = sqrt((1.0f - Xi.y) / (1.0f + (a*a - 1.0f) * Xi.y));
	float SinTheta = sqrt(1.0f - CosTheta * CosTheta);
	vec3 H = vec3(0);
	H.setX(SinTheta * cos(Phi));
	H.setY(SinTheta * sin(Phi));
	H.setZ(CosTheta);

	vec3 UpVector = abs(N.getZ()) < 0.999 ? vec3(0, 0, 1) : vec3(1, 0, 0);
	vec3 TangentX = normalize(cross(UpVector, N));
	vec3 TangentY = cross(N, TangentX);
	// Tangent to world space
	return TangentX * H.getX() + TangentY * H.getY() + N * H.getZ();
}

vec3 SpecularIBL(vec3 SpecularColor, float Roughness, vec3 N, vec3 V, cube_texture* envmap)
{
	vec3 SpecularLighting(0);
	const unsigned int NumSamples = 128;
	for (unsigned int i = 0; i < NumSamples; i++)
	{
		vec2 Xi = Hammersley(i, NumSamples);
		vec3 H = ImportanceSampleGGX(Xi, Roughness, N);
		vec3 L = 2 * dot(V, H) * H - V;
		float NoV = saturate(dot(N, V));
		float NoL = saturate(dot(N, L));
		float NoH = saturate(dot(N, H));
		float VoH = saturate(dot(V, H));
		if (NoL > 0)
		{
			//vec3 SampleColor = EnvMap.SampleLevel(EnvMapSampler, L, 0).rgb;
			vec3 SampleColor = gamma_to_linear(envmap->value(L), GAMMA_FACTOR);
			//float G = G_Smith(Roughness, NoV, NoL);
			float G = G_Smith_Schlick_GGX(Roughness, NoV, NoL);
			float Fc = pow(1 - VoH, 5);
			vec3 F = (1 - Fc) * SpecularColor + vec3(Fc);
			// Incident light = SampleColor * NoL
			// Microfacet specular = D*G*F / (4*NoL*NoV)
			// pdf = D * NoH / (4 * VoH)
			SpecularLighting += mulPerElem(SampleColor, F * G * VoH / (NoH * NoV));
		}
	}
	return SpecularLighting / NumSamples;
}

vec2 IntegrateBRDF(float Roughness, float NoV)
{
	vec3 V(sqrt(1.0f - NoV * NoV), 0, NoV);

	float A = 0;
	float B = 0;
	const unsigned int NumSamples = 128;
	for (unsigned int i = 0; i < NumSamples; i++) {
		vec2 Xi = Hammersley(i, NumSamples);
		vec3 H = ImportanceSampleGGX(Xi, Roughness, vec3(0, 0, 1));
		vec3 L = 2 * dot(V, H) * H - V;
		float NoL = saturate(L.getZ());
		float NoH = saturate(H.getZ());
		float VoH = saturate(dot(V, H));
		if (NoL > 0)
		{
			//float G = G_Smith(Roughness, NoV, NoL);
			float G = G_Smith_Schlick_GGX(Roughness, NoV, NoL);
			float G_Vis = G * VoH / (NoH * NoV);
			float Fc = pow(1 - VoH, 5);
			A += (1 - Fc) * G_Vis;
			B += Fc * G_Vis;
		}
	}
	return vec2(A, B) / NumSamples;
}

vec3 PrefilterEnvMap(float Roughness, vec3 R, cube_texture* envmap)
{
	vec3 N = R;
	vec3 V = R;
	vec3 PrefilteredColor(0);
	float TotalWeight = 0;
	const unsigned int NumSamples = 128;

	for (unsigned int i = 0; i < NumSamples; i++) {
		vec2 Xi = Hammersley(i, NumSamples);
		vec3 H = ImportanceSampleGGX(Xi, Roughness, N); 
		vec3 L = 2.0f * dot(V, H) * H - V;
		float NoL = saturate(dot(N, L));  
		if (NoL > 0)
		{
			PrefilteredColor += gamma_to_linear(envmap->value(L), GAMMA_FACTOR) * NoL;
			TotalWeight += NoL;
		}
	}
	return PrefilteredColor / TotalWeight;
}

vec3 ApproximateSpecularIBL(vec3 SpecularColor, float Roughness, vec3 N, vec3 V, cube_texture* envmap)
{
	float NoV = saturate(dot(N, V));
	vec3 R = 2 * dot(V, N) * N - V;
	vec3 PrefilteredColor = PrefilterEnvMap(Roughness, R, envmap);
	vec2 EnvBRDF = IntegrateBRDF(Roughness, NoV);
	return mulPerElem(PrefilteredColor, (SpecularColor * EnvBRDF.getX() + vec3(EnvBRDF.getY())));
}