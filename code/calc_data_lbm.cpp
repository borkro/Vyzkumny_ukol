#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
// #include <omp.h>
#include "boundaries.h"
#include "boundaries.cpp"
#include "lbm.h"
#include "lbm.cpp"

const double TAU = 0.93; // dimensionless
const unsigned int FRAMES = 50;

// REAL WORLD PARAMS
const double X_DIM = 0.02;								// in metres
const double Y_DIM = 0.01;								// in metres
const double TIME = 5.0;								// in seconds
const unsigned int N_NODES = 202;						// number of nodes on left boundary
const double N_DENSITY = (double)(N_NODES - 2) / Y_DIM; // in nodes per meter
const double VISCOSITY = 1e-6;							// in m2s-1
const double ACCELERATION = 0.01;						// in ms-2, acceleration (from force, eg. gravity)
const double SAVE_SEC = TIME / (double)FRAMES;			// in seconds (save rho,u every SAVE_SEC seconds)
const double T_DENSITY = 3.0 * VISCOSITY * (N_DENSITY * N_DENSITY) / (TAU - 0.5);
const double RHO = 1e3;

const bool PERIODIC = false; // true for periodic inlet/outlet
const bool WN = true;

// const double T_DENSITY = alfa * N_DENSITY * N_DENSITY;	// in time steps per second
// const double U_NORM = 0.3;								// in ms-1

// LATTICE PARAMS
const unsigned int NX = (unsigned int)(X_DIM * N_DENSITY);
const unsigned int NY = N_NODES;
const unsigned int T = (unsigned int)(TIME * T_DENSITY);
const unsigned int SAVE = (unsigned int)(T / FRAMES);
const unsigned int START = 0;
const unsigned int START_OF_M_SAVE = 1;
const double RHO_DIFF = 0.001;
const double RHO_IN = 1.0 + 0.5 * RHO_DIFF;
const double RHO_OUT = 1.0 - 0.5 * RHO_DIFF;
const double FX = ACCELERATION / (T_DENSITY * T_DENSITY / N_DENSITY); // lattice acceleration
const double FY = 0.0;

unsigned int calcIndex(unsigned int x, unsigned int y, unsigned int q = 0)
{
	return NX * (NY * q + y) + x;
}
bool isInCircle(unsigned int x, unsigned int y, unsigned int sx, unsigned int sy, unsigned int r)
{
	return ((sx - x) * (sx - x) + (sy - y) * (sy - y) <= r * r);
}

bool *boundaryMap = (bool *)malloc(sizeof(bool) * NX * NY);
void initBoundaryBounceBack()
{
	for (unsigned int x = 0; x < NX; x++)
	{
		for (unsigned int y = 0; y < NY; y++)
		{
			if (y == 0 || y == NY - 1)
			// if (y == 0 || y == NY - 1 || isInCircle(x, y, NY / 2, NY / 3, NY / 8))
			//   if (y == 0 || y == NY - 1 || (x <= NX / 3 && y <= 2 * NY / 3) || (x >= 2 * NX / 3 && y >= NY / 3))
			//   if (y == 0 || y == NY - 1 || (x <= NX / 6 && y <= 3 * NY / 4) || (x >= 2 * NX / 6 && x <= 3 * NX / 6 && y >= NY / 4) || (x >= 4 * NX / 6 && x <= 5 * NX / 6 && (y >= NY / 2 || y <= NY / 4)) /* || (x >= NX - 4 && y <= 3 * NY / 4) */)
			//   if (y == 0 || y == NY - 1 || (x == NX / 3 && (y <= 2 * NY / 5 || y >= 3 * NY / 5)))
			//   if (y == 0 || y == NY - 1 || (x == NX / 6 && (y - NY / 4) * (y - NY / 4) >= NY * NY / 10) || isInCircle(x, y, NY / 2, NY / 2, NY / 8))
			//   if (y == 0 || y == NY - 1 || (x == NX / 6 && y <= NY / 8) || isInCircle(x, y, NY / 2, NY / 2, NY / 8))
			//   if (y == 0 || y == NY - 1 || (x >= NX / 6 && x <= 2 * NX / 6 && y <= NY / 6))
			// if (y == 0 || y == NY - 1 || x == NX - 1 || (y == NY / 2 && x <= 3 * NX / 4))
			{
				boundaryMap[calcIndex(x, y)] = true;
			}
			else
			{
				boundaryMap[calcIndex(x, y)] = false;
			}
		}
	}
}

double *velocityProfileInX = (double *)malloc(sizeof(double) * (NY));
double *velocityProfileInY = (double *)malloc(sizeof(double) * (NY));
void velocityProfile()
{
	for (unsigned int y = 0; y < NY; y++)
	{
		// velocityProfileInX[y - 1] = 0.1 * sin((-1. + y) * M_PI / (double)(NY / 10));
		// velocityProfileInX[y - 1] = 0.1 * sin((-1. + (double)(y)) * M_PI / (double)(NY - 3));
		// velocityProfileInX[y - 1] = 0.1;
		// velocityProfileInY[y - 1] = 0.1 * sin((-1. + y) * M_PI / (double)(NY - 2));
		velocityProfileInX[y] = 0.0;
		velocityProfileInY[y] = 0.0;
		/* if (y >= 1 * NY / 5 && y <= 2 * NY / 5)
			velocityProfileInX[y] = 0.1;
		if (y >= 3 * NY / 5 && y <= 4 * NY / 5)
			velocityProfileInX[y] = -0.1; */
		if (y >= 1 * NY / 3 && y <= 2 * NY / 3)
		{
			velocityProfileInX[y] = 0.1;
			velocityProfileInY[y] = 0.1;
		}
	}
}
void saveVelocityProfile()
{
	std::ofstream velocityProfileFile("realVelocityProfile.txt");

	for (unsigned int y = 0; y < NY; y++)
		velocityProfileFile << '\n'
							<< velocityProfileInX[y] << ' ' << velocityProfileInY[y];

	velocityProfileFile.close();
}
void loadVelocityProfile()
{
	std::ifstream velocityProfileFile("realVelocityProfile.txt");

	unsigned int y = 0;
	while (!velocityProfileFile.eof())
	{
		velocityProfileFile >> velocityProfileInX[y];
		velocityProfileFile >> velocityProfileInY[y];
		y++;
	}

	velocityProfileFile.close();
}

void saveParams()
{
	std::ofstream params("data/params.txt");

	params << "X = " << X_DIM << " m\n"
		   << "Y = " << Y_DIM << " m\n"
		   << "T = " << TIME << " s\n"
		   << "dx = " << 1.0 / N_DENSITY << " m\n"
		   << "dt = " << 1.0 / T_DENSITY << " s\n"
		   << "acceleration = " << ACCELERATION << " ms-2\n"
		   << "Fx = " << FX << "\n"
		   << "Fy = " << FY << "\n"
		   << "NX = " << NX << "\n"
		   << "NY = " << NY << "\n"
		   << "T = " << T << "\n"
		   << "N_DENSITY = " << N_DENSITY << "\n"
		   << "T_DENSITY = " << T_DENSITY << "\n"
		   << "viscosity = " << VISCOSITY << " m2s-1\n"
		   << "TAU = " << TAU << "\n"
		   << "saving every " << SAVE_SEC << " seconds";

	params.close();
}

void saveAnalyticalPoiseuille()
{
	const int resolution = 1000;
	std::ofstream poiseuille("poiseuille.txt");

	// double h = (Y_DIM - 1.0 / N_DENSITY);
	const double latticeViscosity = (TAU - 0.5) / 3.0;
	const double latticeH = ((double)NY - 2.0);
	for (unsigned int i = 0; i < resolution; i++)
	{
		double y = (latticeH / (double)(resolution - 1)) * i + 0.5;
		double latticeY = -latticeH / 2.0 + y - 0.5;
		poiseuille << y << " " << T_DENSITY / N_DENSITY * FX / (2.0 * latticeViscosity) * (latticeH * latticeH / 4.0 - latticeY * latticeY) << "\n";
	}

	poiseuille.close();
}

int main()
{
	initBoundaryBounceBack();
	velocityProfile();
	saveParams();

	// saveAnalyticalPoiseuille();

	BounceBack *boundary = new BounceBack(boundaryMap, NX, NY);
	LBM *model = new LBM(NX, NY, T, boundary);
	model->setParams(TAU);
	model->printParams();
	model->setC_u_rho(1.0 / N_DENSITY, 1.0 / T_DENSITY, RHO);
	// model->setForce(FX, FY);
	model->setPeriodic(PERIODIC);
	model->setWetNode(WN);

	// model->setRhoInlet(RHO_IN);
	// model->setRhoOutlet(RHO_OUT);

	model->setVelocityProfile(velocityProfileInX, velocityProfileInY);

	model->init();
	// model->initMomentumFile();
	model->save_velocityProfile();

	double dt = model->getDt();
	for (unsigned int t = 0; t < T; t += dt)
	{
		if (!(t % SAVE) && t >= START)
		{
			model->save_macroscopic(t);
			std::cout << "saved #" << t << std::endl;
		}
		model->iterateModel();
		/* if (t >= START_OF_M_SAVE)
			model->save_momentum(1, 1, t); */
	}
	// model->save_velocity_at(NX / 2);
	// std::cout << "saved #" << T << std::endl;

	delete model;
	delete boundary;

	free(boundaryMap);
	free(velocityProfileInX);
	free(velocityProfileInY);

	return 0;
}