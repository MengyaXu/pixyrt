#pragma once

typedef struct _thread_arg {
	int cpuid;
	image_t* img;
	scene_data_t* scene;
	workarea_t warea;
} thread_arg_t;

void trace_ray_thread(thread_arg_t* arg)
{
	image_t* img = arg->img;
	scene_data_t* scene = arg->scene;
	workarea_t* warea = &arg->warea;
	trace_ray(img, scene, warea, 0.0f, 0.5f, kOutputDefault);
}

void init_trace_thread(int thread_num, thread_arg_t* targ, workarea_t* warea, image_t* img, scene_data_t* scene, int w, int h) {
	int i;
	for (i = 0; i < thread_num; ++i) {
		set_workarea(&warea[i], 0, w, i*h / thread_num, (i + 1)*h / thread_num);
		targ[i].img = img;
		targ[i].scene = scene;
		targ[i].warea = warea[i];
	}
}

void start_trace_thraed(int thread_num, HANDLE* handles, thread_arg_t* targ) {
	int i;
	for (i = 0; i < thread_num; ++i) {
		handles[i] = ::CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)trace_ray_thread, (void*)&targ[i], 0, NULL);
		SetThreadIdealProcessor(handles[i], targ[i].cpuid);
	}
}

void wait_trace_thread(int thread_num, HANDLE* handles) {
	WaitForMultipleObjects(thread_num, handles, TRUE, INFINITE);
}