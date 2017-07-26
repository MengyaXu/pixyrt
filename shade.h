#pragma once

vec3 shade(const ray& r, hitable* world, hitable* lights, int depth, int max_depth)
{
	hit_record hrec;
	if (world->hit(r, 0.001, FLT_MAX, hrec)) {
		scatter_record srec;
		vec3 emitted = hrec.matp->emitted(r, hrec, hrec.u, hrec.v, hrec.p);
		if (depth < max_depth && hrec.matp->scatter(r, hrec, srec)) {

			if (srec.is_specular) {
				return mulPerElem(srec.albedo, shade(srec.specular_ray, world, lights, depth + 1, max_depth));
			}
			else if (lights) {
				hitable_pdf pdf0(lights, hrec.p);
				mixture_pdf pdf(&pdf0, srec.pdfp);

				srec.specular_ray = ray(hrec.p, pdf.generate(), r.time());
				float pdf_val = pdf.value(srec.specular_ray.direction());
				delete srec.pdfp;

				vec3 albedo = srec.albedo * hrec.matp->scattering_pdf(r, hrec, srec);
				return emitted + mulPerElem(albedo, shade(srec.specular_ray, world, lights, depth + 1, max_depth)) / pdf_val;
			}
			else {
				return emitted + mulPerElem(srec.albedo, shade(srec.specular_ray, world, lights, depth+1, max_depth));
			}
		}
		else {
			return emitted;
		}
	}
	else {
		return vec3(0);
	}
}