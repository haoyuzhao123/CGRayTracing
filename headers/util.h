#ifndef _UTIL_H_
#define _UTIL_H_

#include <cmath>

// MAX

int MAX(int a, int b) {
    return (a > b) ? a : b;
}

double MAX(double a, double b) {
    return (a > b) ? a : b;
}

// tone mapping and gamma correction
int gammaCorr(double x){
	return int(pow(1-exp(-x),1/2.2)*255+.5);
} 

#endif