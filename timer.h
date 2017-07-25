#pragma once

#include <cstdio>
#include <windows.h>

#define TIME_FROM_1601_TO_1970 (116444736000000000L)

ULONGLONG gettimeval() {
	union {
		ULONGLONG ns100;// time since 1 Jan 1601 in 100ns units
		FILETIME ft;
	} now;
	GetSystemTimeAsFileTime(&now.ft);
	return (now.ns100 - TIME_FROM_1601_TO_1970) / 10;
}

void start_timer(unsigned int* startt) {
	*startt = (unsigned int)gettimeval();
}

unsigned int stop_timer(unsigned int* startt) {
	unsigned int stopt = (unsigned int)gettimeval();
	return (stopt >= *startt) ? (stopt - *startt) : (stopt);
}

//#define print_timer(te) { printf("time of %s:%f[msec]\n", #te, te*1.0e-3); }
#define print_timer(te) { char buf[1024]; sprintf(buf, "time of %s:%f[msec]\n", #te, te*1.0e-3); OutputDebugStringA(buf); }
