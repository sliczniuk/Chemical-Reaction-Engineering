#define BZZ_COMPILER 3
#include <BzzMath.hpp>
#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

void SSDIST(BzzVector &x, BzzVector &f);
void DYNDIST(BzzVector &x, double t, BzzVector &f);

//DATA
double alpha = 1.5, mt = 0.5, zf = 0.5, R = 2.7, Vs = 3.2, D = 0.5;
int n = 30, nf = 15, F = 1, q = 1, md = 5, mb = 5, t0 = 1, tf = 801;

// Dynamic variables
BzzVector RVec(tf), FVec(tf), VSVec(tf), ZFVec(tf);

void main(void)
{
	// Vector required to get a solution vector of SS
	BzzVector f;

	// Initial Vector of concentrations
	BzzVector x0(n);
	for (int i = 1; i <= n; ++i) { x0[i] = 0.5; }
	
	// Solution of Steady State
	BzzNonLinearSystem nls(x0, SSDIST);
	nls();
	nls.GetSolution(&x0, &f);
	nls.BzzPrint("NLS");
	
	// Vector of time
	BzzVector tspan(tf);
	tspan.Linearize(t0, tf);
	double delta_t = tspan[2] - tspan[1];

	// Dynamic Variables
	for (int i = 1; i <= tf; i++) { RVec(i) = R, FVec(i) = F, VSVec(i) = Vs, ZFVec(i) = zf; }
	
	// for (int i = 350; i <= tf; i++) { ZFVec(i) = 0.8; }

	// Dynamic response
	BzzOdeNonStiff o(x0, t0, DYNDIST);
	o(tf);
	o.BzzPrint("Dyn");
	
	// Storage for recording values
	BzzMatrix xd(n, tf); xd.InsertColumn(1, x0);
	BzzVector op(tf), pv(tf), e(tf), ie(tf), dpv(tf), P(tf), I(tf), D(tf), sp(tf);
	for (int i = 1; i <= n; ++i) { op(i) = RVec(i); }
	for (int i = 1; i <= n; ++i) { sp(i) = 0.92; }
	sp(1) = 0;

	// PID Controller - tunning parameters
	double Kc = 50, tauI = 15, tauD = 0.0;

	// Upper and Lower limits on OP / RR
	double op_hi = 10;
	double op_lo = 1.0;

	// Dormammu, I've come to bargain | #TimeLoop #DrStrange
		
	int j = 1;
	while (j<=tf-1)
	{
		e(j) = sp(j) - pv(j);
		if (j>=2)
		{
			dpv(j) = (pv(j) - pv(j - 1)) / delta_t;
			ie(j) = ie(j - 1) + e(j) * delta_t;
		}

		P(j)  =  Kc * e(j);
		I(j)  =  Kc / tauI * ie(j);
		D(j)  = -Kc * tauD * dpv(j);
		op(j) =  op(1) + P(j) + I(j) + D(j);

		if (op(j)>op_hi)		// Check upper limit
		{
			op(j) = op_hi;
			ie(j) = ie(j) - e(j) * delta_t;
		}

		if (op(j) < op_lo)		// Check lower limit
		{
			op(j) = op_lo;
			ie(j) = ie(j) - e(j) * delta_t;
		}

		RVec(j) = op(j);

		BzzOdeNonStiff o(x0, tspan(j), DYNDIST);
		f = o(tspan(j+1));
		f.BzzPrint("f");
		
		// xd.InsertColumn(j+1, f);
		x0 = f;
		
		if (j<tf)
		{
			pv(j + 1) = f(1);
		}
		
		cout << "t = "<< j << "\txd = "<<f(1) << "\tOP = " << op(j) << endl;
		j++;
	}
}


void SSDIST(BzzVector &x, BzzVector &f)
{
	double Lr = R, B = F - D, Ls = R + F*q, Vs = Ls - B, Vr = Vs + F*(1 - q);
	BzzVector y(n);

	for (int i = 1; i <= n; i++) {
		y[i] = alpha*x[i] / (1 + (alpha - 1)*x[i]);
	}

	f(1) = (Vr*y(2) - (D + R)*x(1)); // condenser

	for (int i = 2; i <= nf - 1; i++) {
		f(i) = Lr*x(i - 1) + Vr*y(i + 1) - Lr*x(i) - Vr*y(i);  //rec sec
	}

	f(nf) = Lr*x(nf - 1) + Vs*y(nf + 1) - Ls*x(nf) - Vr*y(nf) + F*zf;  //feed

	for (int i = nf + 1; i <= n - 1; i++) {
		f(i) = Ls*x(i - 1) + Vs*y(i + 1) - Ls*x(i) - Vs*y(i);  //str sec
	}

	f(n) = (Ls*x(n - 1) - B*x(n) - Vs*y(n)); // reboiler
}

void DYNDIST(BzzVector &x, double t, BzzVector &f)
{
	// fstream myfile;
	// myfile.open("Results1.txt", ios_base::app);
	

	R  = RVec[t];
	F  = FVec[t];
	Vs = VSVec[t];
	zf = ZFVec[t];

	double Lr = R, B = F - D, Ls = R + F*q, Vs = Ls - B, Vr = Vs + F*(1 - q);

	BzzVector y(n);

	for (int i = 1; i <= n; i++) {
		y[i] = alpha*x[i] / (1 + (alpha - 1)*x[i]);
	}

	f(1) = Vr*(y(2) - x(1)) / md;	// Condenser

	for (int i = 2; i <= nf - 1; i++) {
		f(i) = (Lr*x(i - 1) + Vr*y(i + 1) - Lr*x(i) - Vr*y(i)) / mt;  // Rec Sec
	}

	f(nf) = (Lr*x(nf - 1) + Vs*y(nf + 1) + F*zf - Ls*x(nf) - Vr*y(nf)) / mt; // Feed

	for (int i = nf + 1; i <= n - 1; i++) {
		f(i) = (Ls*x(i - 1) + Vs*y(i + 1) - Ls*x(i) - Vs*y(i)) / mt;		// Strip Sec
	}

	f(n) = (Ls*x(n - 1) - B*x(n) - Vs*y(n)) / mb; // Boiler

	// cout << x(1) << endl;
	// myfile << t << "a" << x(1) << endl;
	// myfile.close();
}
