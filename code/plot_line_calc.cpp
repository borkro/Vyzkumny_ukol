#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <random>
// #include <omp.h>
#include "boundaries.h"
#include "boundaries.cpp"
#include "lbm.h"
#include "lbm.cpp"
#include "d_adjoint_lbm.h"
#include "d_adjoint_lbm.cpp"

//* ADJOINT VARIABLES
const int MAX_ITER = 1;
const double PRECISION = 0.0;	  //! if cost function is smaller than PRECISION, end GD
const double NORMCHANGE = 1e-2;	  //! if norm of gradient is smaller than NORMCHANGE, end iterations
const int HISTORY = 10;			  //! how many iterations back should norm change check look and average
const double EPSILON = 3e-3;	  //! eps * gradient = step in GD
const double FORGET = 0.2;		  //! on which percentage of the domain (on the x-axis) is data not saved in "measurement"
const double NOISE_STD_DEV = 0.0; //! noise in pre-calculated data

const bool CALCULATEMOMENTUMDATA = true; // set to false if data already calculated
const bool CALCULATEADJOINT = true;

//* LBM PARAMS
const bool PERIODIC = false;
const bool WN = false; // WET NODE

//* LATTICE PARAMS - dimensionless
const double TAU = 3.5; // affects viscosity (with resolution and time step))
const int NX = 100;
const int NY = 100;
const int T = 2000;
const int SAVE = 100;
const int START = 0;
const int START_OF_M_SAVE = 1000;
const int SAVE_EVERY = 1; // saves the desired rho,u only once in every SAVE_EVERY*SAVE_EVERY square

double *initGuessX = (double *)malloc(sizeof(double) * NY);
double *initGuessY = (double *)malloc(sizeof(double) * NY);

double normOfChangeArray[HISTORY] = {0.0};
double updateNormOfChangeArray(double);

int calcIndex(int, int, int);

void folders();

bool *boundaryMap = (bool *)malloc(sizeof(bool) * NX * NY);
void initBoundaryBounceBack(); // define bounce back boundaries
bool isInCircle(int, int, int, int, int);

double *velocityProfileInX = (double *)malloc(sizeof(double) * (NY));
double *velocityProfileInY = (double *)malloc(sizeof(double) * (NY));
void velocityProfile();
void saveVelocityProfile();
void loadVelocityProfile();
void loadInitGuess();

double *dirX = (double *)malloc(sizeof(double) * (NY));
double *dirY = (double *)malloc(sizeof(double) * (NY));
void calcGuessDir();

void saveParams();

int main()
{
	folders();

	// loadVelocityProfile();
	velocityProfile();
	loadInitGuess();
	saveVelocityProfile();

	initBoundaryBounceBack();
	saveParams();

	/***
	 ** LBM CALCULATION PART
	 *
	 */

	BounceBack *boundary = new BounceBack(boundaryMap, NX, NY);
	LBM *model = new LBM(NX, NY, T, boundary);
	model->setParams(TAU);
	model->setPeriodic(PERIODIC);
	model->setWetNode(WN);
	model->set_noise_gen(NOISE_STD_DEV);

	model->setVelocityProfile(velocityProfileInX, velocityProfileInY);

	if (CALCULATEMOMENTUMDATA)
	{
		model->init();
		model->initTimesFile();

		double dt = model->getDt();
		for (int t = 0; t < T; t += dt)
		{
			if (!(t % SAVE) && t >= START)
			{
				model->save_macroscopic(t);
				std::cout << "saved #" << t << std::endl;
			}
			model->iterateModel();
			if (t >= START_OF_M_SAVE)
				model->save_momentum_noise(SAVE_EVERY, SAVE_EVERY, t, SAVE, FORGET);
		}
		model->save_macroscopic(T);
		std::cout << "saved #" << T << std::endl;

		std::cout << std::endl
				  << std::endl;
	}

	/***
	 ** ADJOINT CALCULATION PART
	 *
	 */
	if (CALCULATEADJOINT)
	{
		std::cout << "eps = " << EPSILON << std::endl;

		DiscreteAdjointLBM *adjoint = new DiscreteAdjointLBM(model, EPSILON, initGuessX, initGuessY);

		// adjoint->printGuess();
		// adjoint->printNormOfGradient();

		adjoint->initIterationFile();
		adjoint->saveGuessResults(0);

		for (int i = 1; adjoint->getTempCost() >= PRECISION && i <= MAX_ITER && adjoint->getEps() > 1e-10; i++)
		{
			/* adjoint->iterate(i, false, true);

			// adjoint->printGuess();
			adjoint->printNormOfGradient();
			if (i == 1)
			{
				adjoint->printCostFunction();
			}
			std::ostringstream s;
			s << std::setprecision(10) << adjoint->getTempCost();
			std::cout << "Cost function = " << s.str() << std::endl;

			if (i == 1)
				adjoint->saveIterationResults(adjoint->getCost());
			adjoint->saveIterationResults(adjoint->getTempCost());
			adjoint->saveNormOfGD(adjoint->getNormOfGradient());
			adjoint->saveGuessResults(i);

			double noc = updateNormOfChangeArray(adjoint->getNormOfGradient());
			std::cout << "Norm of gradient (avg over last " << HISTORY << " iterations) = " << noc << std::endl;
			if (noc < NORMCHANGE && i >= HISTORY)
				break; */

			//! plot line search cost function
			{
				velocityProfile();
				loadInitGuess();
				calcGuessDir();
				adjoint->plotOverLine(initGuessX, initGuessY, dirX, dirY, 10);
			}
		}

		delete adjoint;
	}

	delete model;
	delete boundary;

	free(velocityProfileInX);
	free(velocityProfileInY);
	free(initGuessX);
	free(initGuessY);

	return 0;
}

int calcIndex(int x, int y, int q = 0)
{
	return NX * (NY * q + y) + x;
}
bool isInCircle(int x, int y, int sx, int sy, int r)
{
	return ((sx - x) * (sx - x) + (sy - y) * (sy - y) <= r * r);
}

void initBoundaryBounceBack()
{
	for (int x = 0; x < NX; x++)
	{
		for (int y = 0; y < NY; y++)
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
			// if (y == 0 || y == NY - 1 || (x >= NX / 2 && y < 17 * NY / 30 && x <= 3 * NX / 4))
			if (y == 0 || y == NY - 1 || x == NX - 1)
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
	double A = 0.1;
	double ampl = -4.0 * A / ((double)((NY - 1) * (NY - 1)));

	for (int y = 0; y < NY; y++)
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
		// if (y >= 1 * NY / 4 && y <= 3 * NY / 4)
		// if (y >= 3 * NY / 7 && y <= 4 * NY / 7)
		/* {
			velocityProfileInX[y] = 0.15;
			velocityProfileInY[y] = 0.15;
		} */

		velocityProfileInX[y] = 0.0;
		// if (y < NY / 2)
		// velocityProfileInY[y] = ampl * ((double)y * ((double)y - (double)(NY / 2)));

		velocityProfileInY[y] = A;
		if (y < (int)(0.15 * (NY - 1)) + 1)
			velocityProfileInY[y] = A * std::sin(M_PI * (double)y / 0.3 / (double)(NY - 1)) * std::sin(M_PI * (double)y / 0.3 / (double)(NY - 1));
		if (y > (int)(0.85 * (NY - 1)))
			velocityProfileInY[y] = A * std::sin(M_PI * ((double)y - (double)(NY - 1)) / 0.3 / (double)(NY - 1)) * std::sin(M_PI * ((double)y - (double)(NY - 1)) / 0.3 / (double)(NY - 1));
		/* if (y <= 1 || y >= NY - 2)
			velocityProfileInY[y] = 0.0; */

		initGuessX[y] = 0.0;
		initGuessY[y] = 0.0;
	}
}
void saveVelocityProfile()
{
	std::ofstream velocityProfileFile("guessResults/realVelocityProfile.txt");

	for (int y = 0; y < NY; y++)
		velocityProfileFile << '\n'
							<< velocityProfileInX[y] << ' ' << velocityProfileInY[y];

	velocityProfileFile.close();
}
void loadVelocityProfile()
{
	std::ifstream velocityProfileFile("guessResults/realVelocityProfile.txt");

	int y = 0;
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
	std::ifstream velocityProfileFile("initGuess.txt");
	if (!velocityProfileFile)
		std::cout << "velocity profile file not found" << std::endl;

	int y = 0;
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
	std::ofstream adjointParams("data/adjointParams.txt");

	adjointParams << "NX = " << NX << "\n"
				  << "NY = " << NY << "\n"
				  << "T = " << T << "\n"
				  << "TAU = " << TAU << "\n"
				  << "PERIODIC = " << PERIODIC << "\n"
				  << "EPSILON = " << EPSILON << "\n"
				  /* << "Fx = " << FX << "\n"
				  << "Fy = " << FY << "\n" */
				  /* << "initGuessX = " << initGuessX << "\n"
				  << "initGuessY = " << initGuessY << "\n" */
				  << "precision = " << PRECISION << "\n"
				  << "save every " << SAVE_EVERY << "th node"
				  << "\n"
				  << "start of momentum saving = " << START_OF_M_SAVE << "\n"
				  << "noise std dev = " << NOISE_STD_DEV << "\n"
				  << "forget ratio = " << FORGET << "\n";

	adjointParams.close();
}

void folders()
{
	system("rm -rf errorData");
	system("rm guessResults/guessResults*");
	system("rm -rf adjDataVis");
	if (CALCULATEMOMENTUMDATA)
	{
		system("rm -rf data");
		system("rm -rf noiseData");
		system("rm -rf adjointData");
		system("rm -rf adjointTempDataAll");
	}

	std::filesystem::create_directories("noiseData");
	std::filesystem::create_directories("data");
	std::filesystem::create_directories("adjDataVis");
	std::filesystem::create_directories("errorData");
	std::filesystem::create_directories("adjointData");
	std::filesystem::create_directories("adjointTempDataAll");
}

double updateNormOfChangeArray(double value)
{
	double result = value;
	for (int i = 1; i < HISTORY; i++)
	{
		result += normOfChangeArray[i];
		normOfChangeArray[i - 1] = normOfChangeArray[i];
	}
	normOfChangeArray[HISTORY - 1] = value;
	return result / (double)HISTORY;
}

void calcGuessDir()
{
	for (int y = 0; y < NY; y++)
	{
		dirX[y] = initGuessX[y] - velocityProfileInX[y];
		dirY[y] = initGuessY[y] - velocityProfileInY[y];
	}
}