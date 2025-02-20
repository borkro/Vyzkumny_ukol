#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <random>
// #include <omp.h>
#include "boundaries.h"

// D2Q9

class LBM
{
private:
	// D2Q9 params
	const int Q = 9;
	const double w[9] = {4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.};
	const int dir_x[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	const int dir_y[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	const int oppositeOf[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
	const double cs2inv = 3.;

	// lattice params
	int NX = 0;
	int NY = 0;
	double *f_old;
	double *f_new;
	double *f_eq;
	double *rho;
	double rho0 = 1.0;
	double *ux;
	double *uy;
	int dt = 1;
	double tau;
	int t_end;
	double omega = (double)dt / tau;
	double omega_l = 1. - omega;
	double c = 1.0 / (double)dt;
	double *S; // force term

	// Boundary
	BounceBack *boundary;

	// conversion factors
	double C_u = 1.0;
	double C_rho = 1.0;

	// driving scheme
	bool WN = false;	   // wet-node inlet
	bool periodic = false; // periodic inlet/outlet
	double rhoInlet = 1.0;
	double rhoOutlet = 1.0;
	double *velocityX;
	double *velocityY;
	double forceX = 0.0;
	double forceY = 0.0;

	// PRIVATE CALCs
	void updMacro(int x, int y);
	void calcEq(int x, int y, int k);
	void calcS(int x, int y, int k);
	void collide(int x, int y, int k);
	void streamOn(int x, int y, int k);

	// random noise generation
	std::random_device rd;
	std::mt19937 gen;
	std::normal_distribution<double> normal_d{1.0, 0.0};

	// index
	int calcIndex(int x, int y, int q = 0)
	{
		return NX * (NY * q + y) + x;
	}

public:
	// CONSTRUCTOR & DESTRUCTOR
	LBM(int x, int y, int T, BounceBack *b);
	~LBM();
	double random_noise();

	// GETTERS
	double *getRho()
	{
		return rho;
	};
	double *getUx() { return ux; };
	double *getUy() { return uy; };
	int getDt() { return dt; };
	int getNX() { return NX; };
	int getNY() { return NY; };
	int getT_end() { return t_end; };
	double getTau() { return tau; };
	double getRho0() { return rho0; };
	double getCs2inv() { return cs2inv; };
	BounceBack *getBoundary() { return boundary; };
	bool getWN() { return WN; };
	bool getPeriodic() { return periodic; };

	// SETTERS
	void setVelocityProfile(double *vx, double *vy)
	{
		velocityX = vx;
		velocityY = vy;
	};
	void setRhoInlet(double newRhoInlet) { rhoInlet = newRhoInlet; };
	void setRhoOutlet(double newRhoOutlet) { rhoOutlet = newRhoOutlet; };
	void setForce(double Fx, double Fy);
	void setWetNode(bool WetNode) { WN = WetNode; }; // WetNode = true for wet node inlet and outlet
	void setPeriodic(bool p) { periodic = p; };
	void setParams(double newTau);
	void printParams();
	void setC_u_rho(double realDx, double realDt, double realRho);
	void set_noise_gen(double noiseStdDev);

	// files
	void initTimesFile();

	// SIM calculations
	void init();
	void update_macroscopic();
	void calc_eq();
	void collision();
	void stream();
	void streamBulk();	// not done
	void streamInOut(); // not done
	void swap();
	void upd_calc_coll();

	// SAVE MACROSCOPIC PROPERTIES
	void save_macroscopic(int n);
	void save_momentum(int dx, int dy, int n); // IN LATTICE UNITS !!!
	void save_momentum_to(int x, int y, int n, double *rhoTrg, double *uxTrg, double *uyTrg, int index);
	void save_momentum_noise(int dx, int dy, int n, int save_every, double forget);
	void save_velocityProfile();
	void save_velocity_at(int x);

	void save_screenshot_for_error(int n);

	// RUN SIM
	void runModel();
	void iterateModel();
};
