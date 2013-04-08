#include <sys/time.h>
#include "common.h"
Real f_funk(int x, int y, Real h, Real pi){
	
	Real ix, iy;
	ix = (x+1)*h;
	iy = (y+1)*h;
	//double r = exp(x)*sin(2*PI*x)*sin(PI*y);
	Real r = 5*pi*pi*sin(pi*ix)*sin(2*pi*iy);
	return r;
}
Real g_funk(int x, int y, Real h, Real pi){
	Real ix, iy;
	ix = (x+1)*h;
	iy = (y+1)*h;
	Real r = sin(pi*ix)*sin(2*pi*iy);
	return r;
}
Real absolute(Real i){
	if (i < 0.0) return -i;
	else return i;
}
Real wallTime(){
	struct timeval t;
	double msec = 0;
	gettimeofday(&t,0x0);
	msec = t.tv_usec / 1000.0 + t.tv_sec*1000.0;
	
	return msec;

}