#ifndef _UTIL_H_
#define _UTIL_H_

#include <cmath>

// MAX

int max(int a, int b) {
    return (a > b) ? a : b;
}

double max(double a, double b) {
    return (a > b) ? a : b;
}

double max(double a, double b, double c) {
    if (a > b && a > c) {
        return a;
    }
    else if (b > c) {
        return b;
    }
    else {
        return c;
    }
}

// tone mapping and gamma correction
int gammaCorr(double x){
	return int(pow(1-exp(-x),1/2.2)*255+.5);
} 

#endif