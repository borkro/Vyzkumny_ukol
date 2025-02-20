#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include "lbm.h"
#include "boundaries.h"

class DiscreteAdjointLBM
{
private:
	// adjoint variables
	double eps = 1.0;
	LBM *LBMmodel;
	double cost = 0.0;
	double tempCost = 0.0;
	// adjoint with velocity profile
	double *guessX;
	double *guessY;

	// lagrangian multiplicators and terms for calculation
	double *pLagrange;
	double *pLagrange_new;
	double *pLagrange_pc; // post-collision
	/* double *V; // "force" term */

	// gradient
	double *Gx = 0;
	double *Gy = 0;

	// measured data
	int *times = 0;
	int timesSize;
	int *xPos = 0;
	int *yPos = 0;
	double *rhoDes = 0;
	double *uxDes = 0;
	double *uyDes = 0;
	int sizeOfKnownMomentumArray;

	// lbm sim params
	double *rho = 0;
	double *ux = 0;
	double *uy = 0;
	double rhoOutlet = 1.0;

	// D2Q9 params
	const int Q = 9;
	const double w[9] = {4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.};
	const int dir_x[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	const int dir_y[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	const int oppositeOf[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
	const double cs2inv = 3.0;

	// lattice params
	int NX = 0;
	int NY = 0;
	int dt = 1;
	double tau;
	int t_end;
	double omega = (double)dt / tau;
	double omega_l = 1.0 - omega;
	double c = 1.0 / (double)dt;

	// Boundary
	BounceBack *boundary;

	// conversion factors
	double C_u = 1.0;
	double C_rho = 10e3;

	// driving scheme
	bool WN = false;	   // wet-node inlet
	bool periodic = false; // periodic inlet/outlet

	// private functions file managment
	void loadTimes();
	void loadMomentumData(int n);
	void loadMomentum(int n);

	// private functions other
	/* void fillZeros(double *pointer,   int size); */
	void initTempModel(LBM *tempModel);
	void calcRhoUxUy(int n);
	void collideP(int x, int y, int q);
	void collidePwithData(int x, int y, int q, int momentumIndex);
	void streamP(int x, int y, int q);
	void addG(int x, int y);
	void newGuess(int n, bool line_search = true);

	// index
	int calcIndex(int x, int y, int q = 0)
	{
		return NX * (NY * q + y) + x;
	}

public:
	DiscreteAdjointLBM(LBM *model, double epsilon, double *initGuessX, double *initGuessY);
	~DiscreteAdjointLBM();

	double getCost() { return cost; }
	double getTempCost() { return tempCost; }
	double getNormOfGradient();
	double *getGDx() { return Gx; }
	double *getGDy() { return Gy; }

	double getEps() { return eps; }
	void setEps(double e) { eps = e; }

	void initIterationFile();
	void saveIterationResults(double iterationCost);
	void saveNormOfGD(double norm);
	// void initGuessFile();
	void saveGuessResults(int n);

	void initP();
	void collideP();
	void streamP();
	void streamInOutP();
	void swapP();

	void saveState(int n);

	void iterate(int n, bool line_search = true, bool save_state = false);

	void adjustGuess(double newEps, double ratio = 1.0);

	void plotOverLine(double *initGuessX, double *initGuessY, double *dirX, double *dirY, int n_of_steps = 100);

	void printGuess();
	void printNormOfGradient();
	void printCostFunction();

	void costFunction();
	double costFunction(int momentumIndex);
	void tempCostFunction(int t);
};