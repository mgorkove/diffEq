#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std; 

double f(double tn, double yn) {
	return (3 + tn - yn); 
}

//improved euler
double ie(double tn, double h, double tIncrement) {
	if (tn - tIncrement < 0) { 
		double yn = 1.0; 
		return yn; 
	} else {
		double yn = recursion(tn - tIncrement, h, tIncrement) 
		+ (f((tn - tIncrement), recursion(tn - tIncrement, h, tIncrement)) 
		+ f((tn - tIncrement) + h, 
			recursion(tn - tIncrement, h, tIncrement) 
			+ h*f((tn - tIncrement), recursion(tn - tIncrement, h, tIncrement))) )*(h/2.0);
		 
		return yn;  
	}
}
  
//runge-kutta
double rk(double tn, double h) {
	if (tn - h < 0) { 
		double yn = 1.0; 
		return yn; 
	} else { 
		double kn1 = f(tn - h, rk(tn - h, h)); // y_n-1 = f(tn - h, rk(tn - h, h))
		double kn2 = f(tn - h + 0.5*h, rk(tn - h, h) + 0.5*h*kn1); 
		double kn3 = f(tn - h + 0.5*h, rk(tn - h, h) + 0.5*h*kn2);
		double kn4 = f(tn - h + h,   rk(tn - h, h) + h*kn3);  
		double yn = rk(tn - h, h) +  (h/6)*(kn1 + 2*kn2 + 2*kn3 + kn4); 
		return yn; 
	} 
}

int main() {
	double yn = 0; 
	double h = 0.05; double h1 = 0.025; double h2 = 0.0125; double h3 = 0.1; 
	double tn = 0.1; double tIncrement = 0.1; double tStart = 0.1; double tMax = 0.4; 
	int nTimes = (tMax - tStart)/tIncrement + 1; 
	cout << setw(2) << "t" << setw(22) << "Euler with h = 0.05" << setw(22) << "Euler with h = 0.025" << setw(22) << "Euler with h = 0.0125" 
	<< setw(22) << "RK with h = 0.1" << setw(22) << "RK with h = 0.05" <<endl; 
	for (int i = 0; i<nTimes; i++) {
		cout << setw(2) << tn << setw(22) << ie(tn, h, h)<< setw(22) << ie(tn, h1, h1) << setw(22) <<ie(tn, h2, h2) 
		<< setw(22) << rk(tn, h3) << setw(22) << rk(tn, h)<< endl; 
		tn+= tIncrement; 
	} 
}



