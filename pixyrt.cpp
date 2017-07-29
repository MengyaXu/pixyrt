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

texture* load_texture(const char* name)
{
	int nx, ny, nn;
	unsigned char* texels = stbi_load(name, &nx, &ny, &nn, 0);
	return new image_texture(texels, nx, ny);
}

#include "render.h"
#include "trace.h"
#include "plot.h"

int main(int argc, char* argv[]) 
{
	//plot();
	render_scene();
	//raytrace();
	return 0;
}