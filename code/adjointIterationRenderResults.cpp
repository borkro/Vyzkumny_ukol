#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <filesystem>
// #include <omp.h>
#include "boundaries.h"
#include "boundaries.cpp"
#include "lbm.h"
#include "lbm.cpp"
#include "d_adjoint_lbm.h"
#include "d_adjoint_lbm.cpp"

const unsigned int ITERATION = 111;

const double PRECISION = 1e1;
const unsigned int MAX_ITER = 1000;

const bool calculateMomentumData = true; // set to false if data already calculated
const bool calculateAdjoint = false;
const bool continueCalculation = false;
const unsigned int N = 38; // # iteration if continueCalculation = true

const unsigned int SAVE_EVERY = 1; // saves the desired rho,u only once in every SAVE_EVERY*SAVE_EVERY square

// LBM PARAMS
const double NOISE_POWER = 0.05;
const bool PERIODIC = false;
const bool WN = false; // WET NODE - INLET EQ

// LATTICE PARAMS - dimensionless
const double TAU = 1.0;
const unsigned int NX = 100;
const unsigned int NY = 45;
const unsigned int T = 500;
const unsigned int SAVE = 10;
const unsigned int START = 0;
const unsigned int START_OF_M_SAVE = 0;
const double FX = 0.0;
const double FY = 0.0;

// ADJOINT VARIABLES
const double EPSILON = 5e-6;
double *initGuessX = (double *)malloc(NY);
double *initGuessY = (double *)malloc(NY);

unsigned int calcIndex(unsigned int, unsigned int);

bool *boundaryMap = (bool *)malloc(sizeof(bool) * NX * NY);
void initBoundaryBounceBack();
bool isInCircle(unsigned int x, unsigned int y, unsigned int sx, unsigned int sy, unsigned int r);

double *velocityProfileInX = (double *)malloc(sizeof(double) * (NY));
double *velocityProfileInY = (double *)malloc(sizeof(double) * (NY));
void velocityProfile();
void saveVelocityProfile();
void loadVelocityProfile(unsigned int);
void loadInitGuess();

void saveParams();

int main()
{

	std::filesystem::create_directories("iterationRes");

	initBoundaryBounceBack();
	BounceBack *boundary = new BounceBack(boundaryMap, NX, NY);

	for (unsigned int iteration = 0; iteration <= ITERATION; iteration++)
	{
		loadVelocityProfile(iteration);

		LBM *model = new LBM(NX, NY, T, boundary);
		model->setParams(TAU);
		// model->printParams();
		// model->setForce(FX, FY);
		model->setPeriodic(PERIODIC);
		model->setWetNode(WN);

		model->setVelocityProfile(velocityProfileInX, velocityProfileInY);

		model->init();

		double dt = model->getDt();
		for (unsigned int t = 0; t < T; t += dt)
			model->iterateModel();

		char filename_rho[128];
		char filename_u[128];
		char format[29];
		int ndigits = floor(log10((double)ITERATION) + 1.0);
		sprintf(format, "iterationRes/%%s%%0%dd.vtk", ndigits);
		sprintf(filename_rho, format, "ro", iteration);
		sprintf(filename_u, format, "u", iteration);

		std::ofstream file_rho(filename_rho);
		std::ofstream file_u(filename_u);

		file_rho << "# vtk DataFile Version 2.0\ndensity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
				 << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
				 << "\nSCALARS density_scalars double 1\nLOOKUP_TABLE default";
		file_u << "# vtk DataFile Version 2.0\nvelocity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
			   << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
			   << "\nVECTORS velocity_vectors double";

		double *rho = model->getRho();
		double *ux = model->getUx();
		double *uy = model->getUy();

		for (unsigned int y = 0; y < NY; ++y)
		{
			for (unsigned int x = 0; x < NX; ++x)
			{
				if (!boundary->isOnBoundary(calcIndex(x, y)))
				{
					file_rho << "\n"
							 << (rho[calcIndex(x, y)]);
					file_u << "\n"
						   << (ux[calcIndex(x, y)]) << " " << (uy[calcIndex(x, y)]) << " 0";
				}
				else
				{
					file_rho << "\n"
							 << "0";
					file_u << "\n"
						   << "0 0 0";
				}
			}
		}

		file_rho.close();
		file_u.close();

		std::cout << "saved #" << iteration << std::endl;

		delete model;
	}

	delete boundary;

	free(boundaryMap);
	free(velocityProfileInX);
	free(velocityProfileInY);
	free(initGuessX);
	free(initGuessY);

	return 0;
}

unsigned int calcIndex(unsigned int x, unsigned int y)
{
	return NX * (NY * 0 + y) + x;
}
bool isInCircle(unsigned int x, unsigned int y, unsigned int sx, unsigned int sy, unsigned int r)
{
	return ((sx - x) * (sx - x) + (sy - y) * (sy - y) <= r * r);
}

void initBoundaryBounceBack()
{
	for (unsigned int x = 0; x < NX; x++)
	{
		for (unsigned int y = 0; y < NY; y++)
		{
			// if (y == 0 || y == NY - 1)
			// if (y == 0 || y == NY - 1 || isInCircle(x, y, NY / 2, NY / 3, NY / 8))
			//  if (y == 0 || y == NY - 1 || (x <= NX / 3 && y <= 2 * NY / 3) || (x >= 2 * NX / 3 && y >= NY / 3))
			//  if (y == 0 || y == NY - 1 || (x <= NX / 6 && y <= 3 * NY / 4) || (x >= 2 * NX / 6 && x <= 3 * NX / 6 && y >= NY / 4) || (x >= 4 * NX / 6 && x <= 5 * NX / 6 && (y >= NY / 2 || y <= NY / 4)) || (x >= NX - 4 && y <= 3 * NY / 4))
			//  if (y == 0 || y == NY - 1 || (x == NX / 3 && (y <= 2 * NY / 5 || y >= 3 * NY / 5)))
			//  if (y == 0 || y == NY - 1 || (x == NX / 6 && (y - NY / 4) * (y - NY / 4) >= NY * NY / 10) || isInCircle(x, y, NY / 2, NY / 2, NY / 8))
			//  if (y == 0 || y == NY - 1 || (x == NX / 6 && y <= NY / 8) || isInCircle(x, y, NY / 2, NY / 2, NY / 8))
			//  if (y == 0 || y == NY - 1 || (x >= NX / 6 && x <= 2 * NX / 6 && y <= NY / 6))
			// if (y == 0 || y == NY - 1 || x == NX - 1 || (y == NY / 2 && x <= 3 * NX / 4))
			if (y == 0 || y == NY - 1 || (x >= NX / 2 && y < 17 * NY / 30 && x <= 3 * NX / 4))
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
			velocityProfileInX[y] = 0.1; */
		if (y >= 1 * NY / 3 && y <= 2 * NY / 3)
		{
			velocityProfileInX[y] = 0.1;
			velocityProfileInY[y] = 0.1;
		}
		initGuessX[y] = 0.0;
		initGuessY[y] = 0.0;
	}
}
void saveVelocityProfile()
{
	std::ofstream velocityProfileFile("guessResults/realVelocityProfile.txt");

	for (unsigned int y = 0; y < NY; y++)
		velocityProfileFile << '\n'
							<< velocityProfileInX[y] << ' ' << velocityProfileInY[y];

	velocityProfileFile.close();
}
void loadVelocityProfile(unsigned int n)
{
	std::ifstream velocityProfileFile("guessResults/guessResults" + std::to_string(n) + ".txt");

	unsigned int y = 0;
	do
	{
		velocityProfileFile >> velocityProfileInX[y];
		velocityProfileFile >> velocityProfileInY[y];
		y++;
	} while (!velocityProfileFile.eof() && y < NY);

	velocityProfileFile.close();
}
void loadInitGuess()
{
	std::ifstream velocityProfileFile("guessResults/guessResults" + std::to_string(N) + ".txt");
	if (!velocityProfileFile)
		std::cout << "velocity profile file not found" << std::endl;

	unsigned int y = 0;
	do
	{
		velocityProfileFile >> initGuessX[y];
		velocityProfileFile >> initGuessY[y];
		y++;
	} while (!velocityProfileFile.eof() && y < NY);

	velocityProfileFile.close();
}

void saveParams()
{
	std::ofstream adjointParams("adjointData/adjointParams.txt");

	adjointParams << "NX = " << NX << "\n"
				  << "NY = " << NY << "\n"
				  << "T = " << T << "\n"
				  << "TAU = " << TAU << "\n"
				  << "PERIODIC = " << PERIODIC << "\n"
				  << "EPSILON = " << EPSILON << "\n"
				  << "Fx = " << FX << "\n"
				  << "Fy = " << FY << "\n"
				  /* << "initGuessX = " << initGuessX << "\n"
				  << "initGuessY = " << initGuessY << "\n" */
				  << "precision = " << PRECISION << "\n"
				  << "save every " << SAVE_EVERY << "th node"
				  << "\n"
				  << "noisePower = " << NOISE_POWER;

	adjointParams.close();
}