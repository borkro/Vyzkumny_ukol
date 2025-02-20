#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include "lbm.h"
#include "boundaries.h"

class AdjointLBM
{
private:
	// adjoint variables
	double eps = 1.0;
	LBM *LBMmodel;
	double cost;
	// adjoint with forces
	double initGuessX;
	double initGuessY;
	double guessX;
	double guessY;

	// lagrangian multiplicators and terms for calculation
	double *pLagrange;
	double *pLagrange_new;
	double *pLagrange_eq;
	double *V; // "force" term

	// gradient
	double Gx;
	double Gy;

	// collision operator
	unsigned int *times = 0;
	unsigned int timesSize;
	unsigned int *xPos = 0;
	unsigned int *yPos = 0;
	double *rhoDes = 0;
	double *uxDes = 0;
	double *uyDes = 0;
	unsigned int sizeOfKnownMomentumArray;

	// lbm sim params
	double *rho = 0;
	double *ux = 0;
	double *uy = 0;

	// D2Q9 params
	const int Q = 9;
	const double w[9] = {4. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 9., 1. / 36., 1. / 36., 1. / 36., 1. / 36.};
	const int dir_x[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	const int dir_y[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	const int oppositeOf[9] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
	const double cs2inv = 3.0;

	// lattice params
	unsigned int NX = 0;
	unsigned int NY = 0;
	unsigned int dt = 1;
	double tau;
	unsigned int t_end;
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

	// private functions
	void initTempModel(LBM *tempModel);
	void calcRhoUxUy();
	void calcPEq(unsigned int x, unsigned int y, unsigned int q, unsigned int momentumIndex);
	void calcV(unsigned int x, unsigned int y, unsigned int q, unsigned int momentumIndex);
	void collideP(unsigned int x, unsigned int y, unsigned int q, unsigned int momentumIndex);
	void streamP(unsigned int x, unsigned int y, unsigned int q);
	void addG(unsigned int x, unsigned int y, unsigned int momentumIndex);
	void addG_overAll(unsigned int time);

	void loadTimes();
	void loadMomentumData(unsigned int n);
	void loadMomentum(unsigned int n);

	// index
	unsigned int calcIndex(unsigned int x, unsigned int y, unsigned int q = 0)
	{
		return NX * (NY * q + y) + x;
	}

public:
	AdjointLBM(LBM *model, double epsilon, double initGuessX, double initGuessY);
	~AdjointLBM();

	double getCost() { return cost; }

	void initIterationFile();
	void saveIterationResults();
	void initGuessFile();
	void saveGuessResults();

	void initP();
	void streamP();
	void swapP();

	void iterate(unsigned int n);

	void printGuess();
	void printNormOfGradient();
	void printCostFunction();

	void costFunction();
	double costFunction(unsigned int momentumIndex);
};