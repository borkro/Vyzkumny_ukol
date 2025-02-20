#include "adjoint_lbm.h"

AdjointLBM::AdjointLBM(LBM *model, double epsilon, double initGuessX, double initGuessY)
	: eps(epsilon), LBMmodel(model), guessX(initGuessX), guessY(initGuessY)
{
	NX = model->getNX();
	NY = model->getNY();
	dt = model->getDt();
	tau = model->getTau();
	t_end = model->getT_end();
	omega = (double)dt / tau;
	omega_l = 1. - omega;
	c = 1.0 / (double)dt;
	boundary = model->getBoundary();
	WN = model->getWN();
	periodic = model->getPeriodic();

	Gx = Gy = 0.0;
	cost = 1e10;

	unsigned const int sizeVector = sizeof(double) * NX * NY * Q;
	pLagrange = (double *)malloc(sizeVector);
	pLagrange_eq = (double *)malloc(sizeVector);
	pLagrange_new = (double *)malloc(sizeVector);
	V = (double *)malloc(sizeVector);

	loadTimes();
}

AdjointLBM::~AdjointLBM()
{
	free(pLagrange);
	free(pLagrange_new);
	free(pLagrange_eq);
	free(V);

	free(times);
	free(xPos);
	free(yPos);
	free(rhoDes);
	free(uxDes);
	free(uyDes);

	free(rho);
	free(ux);
	free(uy);
}

void AdjointLBM::loadTimes()
{
	std::ifstream timesData("adjointData/times.txt");
	unsigned int counter = 0;
	times = (unsigned int *)malloc(sizeof(unsigned int) * counter);
	while (!timesData.eof())
	{
		counter++;

		// realloc
		times = (unsigned int *)realloc(times, sizeof(unsigned int) * counter);

		// save
		timesData >> times[counter - 1];
	}
	timesSize = counter;
	timesData.close();
}

void AdjointLBM::loadMomentumData(unsigned int n)
{
	char filename_adj[128];
	char format[28];
	int ndigits = floor(log10((double)t_end) + 1.0);
	sprintf(format, "adjointData/%%s%%0%dd.txt", ndigits);
	sprintf(filename_adj, format, "adj", n);

	std::ifstream momentumFile(filename_adj);
	if (!momentumFile)
		std::cout << "momentum file not found" << std::endl;

	free(xPos);
	free(yPos);
	free(rhoDes);
	free(uxDes);
	free(uyDes);

	unsigned int counter = 0;
	// times = (unsigned int *)malloc(sizeof(unsigned int) * counter);
	xPos = (unsigned int *)malloc(sizeof(unsigned int) * counter);
	yPos = (unsigned int *)malloc(sizeof(unsigned int) * counter);
	rhoDes = (double *)malloc(sizeof(double) * counter);
	uxDes = (double *)malloc(sizeof(double) * counter);
	uyDes = (double *)malloc(sizeof(double) * counter);
	while (!momentumFile.eof())
	{
		counter++;

		// realloc
		// times = (unsigned int *)realloc(times, sizeof(unsigned int) * counter);
		xPos = (unsigned int *)realloc(xPos, sizeof(unsigned int) * counter);
		yPos = (unsigned int *)realloc(yPos, sizeof(unsigned int) * counter);
		rhoDes = (double *)realloc(rhoDes, sizeof(double) * counter);
		uxDes = (double *)realloc(uxDes, sizeof(double) * counter);
		uyDes = (double *)realloc(uyDes, sizeof(double) * counter);

		// save
		// momentumFile >> times[counter - 1];
		momentumFile >> xPos[counter - 1];
		momentumFile >> yPos[counter - 1];
		momentumFile >> rhoDes[counter - 1];
		momentumFile >> uxDes[counter - 1];
		momentumFile >> uyDes[counter - 1];

		// std::cout << times[counter - 1] << " " << xPos[counter - 1] << " " << yPos[counter - 1] << std::endl;
	}
	momentumFile.close();

	sizeOfKnownMomentumArray = counter;
	/* std::cout << "momentum data loaded SUCCESSFULLY: total " << sizeOfKnownMomentumArray << " points loaded" << std::endl
	<< std::endl; */

	rho = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
	ux = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
	uy = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
}

void AdjointLBM::loadMomentum(unsigned int n)
{
	char filename_adj[128];
	char format[32];
	int ndigits = floor(log10((double)t_end) + 1.0);
	sprintf(format, "adjointTempData/%%s%%0%dd.txt", ndigits);
	sprintf(filename_adj, format, "adj", n);

	std::ifstream momentumFile(filename_adj);
	if (!momentumFile)
		std::cout << "temp model momentum file not found" << std::endl;

	free(rho);
	free(ux);
	free(uy);

	unsigned int counter = 0;
	rho = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
	ux = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
	uy = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
	while (!momentumFile.eof())
	{
		counter++;

		// save
		momentumFile >> rho[counter - 1];
		momentumFile >> ux[counter - 1];
		momentumFile >> uy[counter - 1];

		// std::cout << times[counter - 1] << " " << xPos[counter - 1] << " " << yPos[counter - 1] << std::endl;
	}
	momentumFile.close();

	/* std::cout << "momentum of temp model loaded SUCCESSFULLY: total " << counter << " points loaded" << std::endl
	<< std::endl; */
}

void AdjointLBM::initIterationFile()
{
	std::ofstream costFunctionFile;
	costFunctionFile.open("costFunction.txt");
	costFunctionFile.close();
}

void AdjointLBM::saveIterationResults()
{
	std::ofstream costFunctionFile;
	costFunctionFile.open("costFunction.txt", std::ofstream::app);

	costFunctionFile << getCost() << "\n";

	costFunctionFile.close();
}

void AdjointLBM::initGuessFile()
{
	std::ofstream costFunctionFile;
	costFunctionFile.open("guessResults.txt");
	costFunctionFile.close();
}

void AdjointLBM::saveGuessResults()
{
	std::ofstream costFunctionFile;
	costFunctionFile.open("guessResults.txt", std::ofstream::app);

	costFunctionFile << guessX << " " << guessY << "\n";

	costFunctionFile.close();
}

void AdjointLBM::initTempModel(LBM *tempModel)
{
	tempModel->setParams(tau);
	tempModel->setPeriodic(periodic);
	tempModel->setWetNode(WN);
	tempModel->setForce(guessX, guessY);
	tempModel->init();

	/* for (unsigned int counter = 0; counter < sizeOfKnownMomentumArray; counter++)
	{
		tempModel->save_momentum_to(xPos[counter], yPos[counter], times[counter], rho, ux, uy, counter);
	} */
}

void AdjointLBM::calcRhoUxUy()
{
	LBM *tempModel = new LBM(NX, NY, t_end, boundary);
	initTempModel(tempModel);
	loadMomentumData(times[0]);

	free(rho);
	free(ux);
	free(uy);
	rho = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
	ux = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);
	uy = (double *)malloc(sizeof(double) * sizeOfKnownMomentumArray);

	unsigned int counter = 0;
	for (unsigned int t = 0; t < t_end; t += dt)
	{
		tempModel->iterateModel();

		if (times[counter] == t)
		{
			// loadMomentumData(times[t]);

			char filename_adj[128];
			char format[32];
			int ndigits = floor(log10((double)t_end) + 1.0);
			sprintf(format, "adjointTempData/%%s%%0%dd.txt", ndigits);
			sprintf(filename_adj, format, "adj", t);

			std::ofstream momentumData(filename_adj);

			for (unsigned int counterData = 0; counterData < sizeOfKnownMomentumArray; counterData++)
			{
				// unsigned int index = calcIndex(xPos[counterData], yPos[counterData]);

				tempModel->save_momentum_to(xPos[counterData], yPos[counterData], times[counterData], rho, ux, uy, counterData);

				// momentumData << xPos[counterData] << " " << yPos[counterData] << " ";
				momentumData << "\n";
				momentumData << rho[counterData] << " " << ux[counterData] << " " << uy[counterData];
			}

			momentumData.close();
			counter++;
		}
	}
	delete tempModel;
}

void AdjointLBM::calcPEq(unsigned int x, unsigned int y, unsigned int q, unsigned int momentumIndex)
{
	unsigned int index = calcIndex(x, y);
	if (!boundary->isOnBoundary(index))
	{
		pLagrange_eq[calcIndex(x, y, q)] = 0.0;
		for (unsigned int k = 0; k < Q; k++)
		{
			double dot = ux[momentumIndex] * (double)dir_x[k] + uy[momentumIndex] * (double)dir_y[k];
			double power_u = ux[momentumIndex] * ux[momentumIndex] + uy[momentumIndex] * uy[momentumIndex];

			// derivatives
			double dfeq_drho = w[k] * (1.0 + dot * cs2inv + (dot * dot) * (cs2inv * cs2inv) / 2.0 - power_u * cs2inv / 2.0);
			double dfeq_dux = w[k] * ((double)dir_x[k] * cs2inv + dot * (double)dir_x[k] * (cs2inv * cs2inv) - ux[momentumIndex] * cs2inv);
			double dfeq_duy = w[k] * ((double)dir_y[k] * cs2inv + dot * (double)dir_y[k] * (cs2inv * cs2inv) - uy[momentumIndex] * cs2inv);

			double dotGuess = guessX * (double)dir_x[k] + guessY * (double)dir_y[k];

			/* double dS_drho = (1.0 - 0.5 * (double)dt / tau) * w[k] * cs2inv * ((double)dir_x[k] * guessX + (double)dir_y[k] * guessY + cs2inv * ((double)dir_x[k] * guessX + (double)dir_y[k] * guessY) * dot - ux[index] * guessX - uy[index] * guessY);
			double dS_dux = (1.0 - 0.5 * (double)dt / tau) * w[k] * cs2inv * (-1.0 * guessX + cs2inv * (double)dir_x[k] * dotGuess);
			double dS_duy = (1.0 - 0.5 * (double)dt / tau) * w[k] * cs2inv * (-1.0 * guessY + cs2inv * (double)dir_y[k] * dotGuess); */

			double dF_drho = w[k] * cs2inv * (dotGuess + cs2inv * dotGuess * dot - ux[index] * guessX - uy[index] * guessY);
			double dF_dux = w[k] * cs2inv * (-1.0 * guessX + cs2inv * (double)dir_x[k] * dotGuess);
			double dF_duy = w[k] * cs2inv * (-1.0 * guessY + cs2inv * (double)dir_y[k] * dotGuess);

			// add to p_eq
			pLagrange_eq[calcIndex(x, y, q)] += pLagrange[calcIndex(x, y, k)] * (1.0 * (dfeq_drho + tau * dF_drho) + (double)dir_x[q] * (dfeq_dux + tau * dF_dux) + (double)dir_y[q] * (dfeq_duy + tau * dF_duy));
		}
	}
}

void AdjointLBM::calcV(unsigned int x, unsigned int y, unsigned int q, unsigned int momentumIndex)
{
	unsigned int index = calcIndex(x, y);
	if (!boundary->isOnBoundary(index))
	{
		V[calcIndex(x, y, q)] = -((rho[momentumIndex] - rhoDes[momentumIndex]) * 1.0 + (ux[momentumIndex] - uxDes[momentumIndex]) * (double)dir_x[q] / rho[momentumIndex] + (uy[momentumIndex] - uyDes[momentumIndex]) * (double)dir_y[q] / rho[momentumIndex]);
	}
}

void AdjointLBM::collideP(unsigned int x, unsigned int y, unsigned int q, unsigned int momentumIndex)
{
	unsigned int index = calcIndex(x, y);
	if (!boundary->isOnBoundary(index))
	{
		pLagrange[calcIndex(x, y, q)] = omega_l * pLagrange[calcIndex(x, y, q)] + omega * pLagrange_eq[calcIndex(x, y, q)] + /* (1.0 - 0.5 / tau) * */ V[calcIndex(x, y, q)];
	}
}

void AdjointLBM::streamP(unsigned int x, unsigned int y, unsigned int q)
{
	unsigned int x_dir = (NX + x + dir_x[q]) % NX; // + because opposite streaming
	unsigned int y_dir = (NY + y + dir_y[q]) % NY; // + because opposite streaming

	if (boundary->isOnBoundary(calcIndex(x_dir, y_dir)))
	{
		pLagrange_new[calcIndex(x, y, q)] = 0.0;
		// pLagrange_new[calcIndex(x, y, q)] = pLagrange[calcIndex(x, y, oppositeOf[q])];
	}
	else
	{
		pLagrange_new[calcIndex(x, y, q)] = pLagrange[calcIndex(x_dir, y_dir, q)];
	}
}

void AdjointLBM::initP()
{
	for (unsigned int x = 0; x < NX; x++)
	{
		for (unsigned int y = 0; y < NY; y++)
		{
			for (unsigned int q = 0; q < Q; q++)
			{
				pLagrange[calcIndex(x, y, q)] = 0.0;
				pLagrange_new[calcIndex(x, y, q)] = 0.0;
			}
		}
	}
}

void AdjointLBM::streamP()
{
	for (unsigned int x = 0; x < NX; x++)
	{
		for (unsigned int y = 0; y < NY; y++)
		{
			for (unsigned int q = 0; q < Q; q++)
			{
				streamP(x, y, q);
			}
		}
	}
}

void AdjointLBM::addG(unsigned int x, unsigned int y, unsigned int momentumIndex)
{
	for (unsigned int k = 0; k < Q; k++)
	{
		double dot = ux[momentumIndex] * (double)dir_x[k] + uy[momentumIndex] * (double)dir_y[k];

		/* double dotGuess = guessX * (double)dir_x[k] + guessY * (double)dir_y[k];

		double b1x = cs2inv * (-1.0 * guessX + cs2inv * (double)dir_x[k] * dotGuess) * 0.5;
		double b1y = cs2inv * (-1.0 * guessY + cs2inv * (double)dir_y[k] * dotGuess) * 0.5; */

		double b2x = cs2inv * ((double)dir_x[k] - ux[momentumIndex] + cs2inv * dot * (double)dir_x[k]);
		double b2y = cs2inv * ((double)dir_y[k] - uy[momentumIndex] + cs2inv * dot * (double)dir_y[k]);

		/* double dS_dAx = (1.0 - 0.5 / tau) * w[k] * rho[momentumIndex] * (b1x + b2x);
		double dS_dAy = (1.0 - 0.5 / tau) * w[k] * rho[momentumIndex] * (b1y + b2y); */

		double dF_dAx = w[k] * rho[momentumIndex] * (/* b1x + */ b2x);
		double dF_dAy = w[k] * rho[momentumIndex] * (/* b1y + */ b2y);

		/* double b3x = cs2inv * ((double)dir_x[k] + cs2inv * dot * (double)dir_x[k] - ux[momentumIndex]) * 0.5;
		double b3y = cs2inv * ((double)dir_y[k] + cs2inv * dot * (double)dir_y[k] - uy[momentumIndex]) * 0.5;

		double dfeq_dAx = w[k] * rho[momentumIndex] * b3x;
		double dfeq_dAy = w[k] * rho[momentumIndex] * b3y; */

		Gx += -pLagrange[calcIndex(x, y, k)] * (dF_dAx /* + dfeq_dAx / tau */);
		Gy += -pLagrange[calcIndex(x, y, k)] * (dF_dAy /* + dfeq_dAy / tau */);
	}
}

void AdjointLBM::addG_overAll(unsigned int time)
{
	for (unsigned int counter = 0; counter < sizeOfKnownMomentumArray; counter++)
	{
		addG(xPos[counter], yPos[counter], counter);
	}
	Gx /= (double)sizeOfKnownMomentumArray;
	Gy /= (double)sizeOfKnownMomentumArray;
}

void AdjointLBM::swapP()
{
	double *temp = pLagrange_new;
	pLagrange_new = pLagrange;
	pLagrange = temp;
}

void AdjointLBM::costFunction()
{
	for (unsigned int counter = 0; counter < sizeOfKnownMomentumArray; counter++)
	{
		cost += costFunction(counter);
	}
}

double AdjointLBM::costFunction(unsigned int momentumIndex)
{
	double result = 0.0;
	double rhoDiff = rho[momentumIndex] - rhoDes[momentumIndex];
	double uxDiff = ux[momentumIndex] - uxDes[momentumIndex];
	double uyDiff = uy[momentumIndex] - uyDes[momentumIndex];
	result = 0.5 * (rhoDiff * rhoDiff + uxDiff * uxDiff + uyDiff * uyDiff);
	return result;
}

void AdjointLBM::iterate(unsigned int n)
{
	cost = 0.0;

	std::cout << std::endl
			  << "ITERATION #" << n << std::endl;

	calcRhoUxUy();
	initP();

	unsigned int buffer = 0;
	for (unsigned int t = 0; t < t_end; t += dt)
	{
		// load data
		if (t_end - 1 - t == times[timesSize - 1 - buffer])
		{
			loadMomentumData(t_end - 1 - t);
			loadMomentum(t_end - 1 - t);
			buffer++;

			// add to cost
			costFunction();

			// collide pLagrange where we know CO
			for (unsigned int counterData = 0; counterData < sizeOfKnownMomentumArray; counterData++)
			{
				// std::cout << bufferIndex << std::endl;
				for (unsigned int k = 0; k < Q; k++)
				{
					// calc_eq
					calcPEq(xPos[counterData], yPos[counterData], k, counterData);

					// calc force terms
					calcV(xPos[counterData], yPos[counterData], k, counterData);

					// collide
					collideP(xPos[counterData], yPos[counterData], k, counterData);
				}
			}
		}

		// stream pLagrange
		streamP();

		// swap
		swapP();

		if (t == t_end - 1)
		{
			// add to G
			Gx = 0.0;
			Gy = 0.0;
			/* unsigned int momentumIndex = calcIndex(NX / 2, NY / 2) + NX * NY * 0;
			double xPosition = xPos[momentumIndex];
			double yPosition = yPos[momentumIndex];
			double time = times[momentumIndex];
			addG(xPosition, yPosition, momentumIndex); */
			addG_overAll(t_end - 1);
		}
	}

	// new guess
	std::cout << "Gx = " << Gx << "\tGy = " << Gy << std::endl;
	guessX -= eps * Gx;
	guessY -= eps * Gy;

	// std::cout << "KONEC ITERACE" << std::endl;
}

void AdjointLBM::printGuess()
{
	std::cout << "New guess:" << std::endl
			  << "Ax \t\t Ay" << std::endl
			  << guessX << " \t " << guessY << std::endl;
}

void AdjointLBM::printNormOfGradient()
{
	std::cout << "Norm of change: " << eps * (std::sqrt(Gx * Gx + Gy * Gy)) << std::endl;
}

void AdjointLBM::printCostFunction()
{
	std::cout << "Cost function = " << getCost() << std::endl
			  << "Cost function in x = " << costFunction(calcIndex(xPos[sizeOfKnownMomentumArray / 2], yPos[sizeOfKnownMomentumArray / 2])) << std::endl;
}