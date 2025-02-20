#include "d_adjoint_lbm.h"
#include <string>
#include <sstream>
#include <iomanip>

DiscreteAdjointLBM::DiscreteAdjointLBM(LBM *model, double epsilon, double *initGuessX, double *initGuessY)
	: eps(epsilon), LBMmodel(model), NX(model->getNX()), NY(model->getNY()), dt(model->getDt()), tau(model->getTau()), t_end(model->getT_end()), boundary(model->getBoundary()), WN(model->getWN()), periodic(model->getPeriodic())
{
	omega = (double)dt / tau;
	omega_l = 1.0 - omega;
	c = 1.0 / (double)dt;

	Gx = (double *)malloc(sizeof(double) * NY);
	Gy = (double *)malloc(sizeof(double) * NY);
	guessX = (double *)malloc(sizeof(double) * NY);
	guessY = (double *)malloc(sizeof(double) * NY);
	for (int y = 0; y < NY; y++)
	{
		guessX[y] = initGuessX[y];
		guessY[y] = initGuessY[y];
		Gx[y] = 0.0;
		Gy[y] = 0.0;
	}
	cost = tempCost = 1e10;

	const int sizeVector = sizeof(double) * NX * NY * Q;
	pLagrange = (double *)malloc(sizeVector);
	pLagrange_new = (double *)malloc(sizeVector);
	pLagrange_pc = (double *)malloc(sizeVector);

	loadTimes();
}

DiscreteAdjointLBM::~DiscreteAdjointLBM()
{
	free(pLagrange);
	free(pLagrange_new);
	free(pLagrange_pc);

	free(Gx);
	free(Gy);

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

double DiscreteAdjointLBM::getNormOfGradient()
{
	double result = 0.0;
	for (int y = 0; y < NY; y++)
		result += (/* Gx[y] * Gx[y] + */ Gy[y] * Gy[y]);
	return std::sqrt(result);
}

void DiscreteAdjointLBM::initIterationFile()
{
	std::ofstream costFunctionFile("costFunction.txt");
	costFunctionFile.close();
	std::ofstream normGDFile("normGD.txt");
	normGDFile.close();
}

void DiscreteAdjointLBM::saveIterationResults(double iterationCost)
{
	std::ofstream costFunctionFile;
	costFunctionFile.open("costFunction.txt", std::ofstream::app);

	std::ostringstream s;
	s << std::setprecision(10) << iterationCost;

	costFunctionFile << s.str() << "\n";

	costFunctionFile.close();
}

void DiscreteAdjointLBM::saveNormOfGD(double norm)
{
	std::ofstream normGDFile;
	normGDFile.open("normGD.txt", std::ofstream::app);

	std::ostringstream s;
	s << std::setprecision(10) << norm;

	normGDFile << s.str() << "\n";

	normGDFile.close();
}

/* void DiscreteAdjointLBM::initGuessFile()
{
	std::ofstream costFunctionFile;
	costFunctionFile.open("guessResults.txt");
	costFunctionFile.close();
} */

void DiscreteAdjointLBM::saveGuessResults(int n)
{
	std::ofstream costFunctionFile;
	costFunctionFile.open("guessResults/guessResults" + std::to_string(n) + ".txt");

	for (int i = 0; i < NY; i++)
	{
		std::ostringstream s;
		s << std::setprecision(10) << guessX[i];
		costFunctionFile << s.str();
		s.str("");
		s.clear();
		s << std::setprecision(10) << guessY[i];
		costFunctionFile << " " << s.str() << "\n";
	}

	costFunctionFile.close();
}

void DiscreteAdjointLBM::initP()
{
	for (int i = 0; i < NX * NY * Q; i++)
	{
		pLagrange[i] = 0.0;
		pLagrange_new[i] = 0.0;
		pLagrange_pc[i] = 0.0;
	}

	// calc primary problem
	// calcRhoUxUy();

	// initial condition only where data is measured
	loadMomentumData(t_end - 1);
	loadMomentum(t_end - 1);
	// std::cout << rho[calcIndex(0, 7)] << " " << ux[calcIndex(0, 7)] << " " << uy[calcIndex(0, 7)] << " " << rhoDes[700] << " " << uxDes[700] << " " << uyDes[700] << std::endl;
	for (int i = 0; i < sizeOfKnownMomentumArray; i++)
	{
		int index = calcIndex(xPos[i], yPos[i]);
		if (!boundary->isOnBoundary(index))
		{
			if (rho[index] == 0.0)
				std::cout << xPos[i] << " " << yPos[i] << std::endl;
			for (int q = 0; q < Q; q++)
			{
				/* double result = (rho[index] - rhoDes[i]) * 1.0 + (ux[index] - uxDes[i]) * (double)dir_x[q] / rho[index] + (uy[index] - uyDes[i]) * (double)dir_y[q] / rho[index];
				  int x = xPos[i];
				  int y = yPos[i]; */
				pLagrange[calcIndex(xPos[i], yPos[i], q)] = (rho[index] - rhoDes[i]) * 1.0 + (ux[index] - uxDes[i]) * ((double)dir_x[q] * rho[index] - ux[index]) / rho[index] / rho[index] + (uy[index] - uyDes[i]) * ((double)dir_y[q] * rho[index] - uy[index]) / rho[index] / rho[index];
			}

			// add to cost
			cost += costFunction(i);
		}
	}
}

void DiscreteAdjointLBM::collideP()
{
	for (int x = 0; x < NX; x++)
	{
		for (int y = 0; y < NY; y++)
		{
			if (!boundary->isOnBoundary(calcIndex(x, y)))
			{
				for (int q = 0; q < Q; q++)
				{
					collideP(x, y, q);
				}
			}
		}
	}
}

void DiscreteAdjointLBM::streamP()
{
	for (int x = 0; x < NX; x++)
	{
		for (int y = 0; y < NY; y++)
		{
			if (!boundary->isOnBoundary(calcIndex(x, y)))
			{
				for (int q = 0; q < Q; q++)
				{
					streamP(x, y, q);
				}
			}
		}
	}
	if (WN)
	{
		streamInOutP();
	}
}

void DiscreteAdjointLBM::streamInOutP()
{
	for (int y = 0; y < NY; y++)
	{
		int x = 0;
		if (!boundary->isOnBoundary(calcIndex(x, y)))
		{
			/* for (  int k = 0; k < Q; k++)
			{
				pLagrange[calcIndex(x, y, k)] = pLagrange_new[calcIndex(x, y, k)];
			} */
			for (int q = 0; q < Q; q++)
			{
				for (int k = 0; k < Q; k++)
				{
					double dot = guessX[y] * (double)dir_x[k] + guessY[y] * (double)dir_y[k];
					double power_u = guessX[y] * guessX[y] + guessY[y] * guessY[y];
					double dfeq_drho_w = w[k] * (1.0 + dot * cs2inv + (dot * dot) * (cs2inv * cs2inv) / 2.0 - power_u * cs2inv / 2.0);

					if (q == 0 || q == 2 || q == 4)
					{
						double drho_w_df = 1.0 / (1.0 - guessX[y]);
						pLagrange_new[calcIndex(x, y, q)] += pLagrange[calcIndex(x, y, k)] * dfeq_drho_w * drho_w_df;
					}
					else if (q == 3 || q == 6 || q == 7)
					{
						double drho_w_df = 2.0 / (1.0 - guessX[y]);
						pLagrange_new[calcIndex(x, y, q)] += pLagrange[calcIndex(x, y, k)] * dfeq_drho_w * drho_w_df;
					}
				}
			}
		}
	}
}

void DiscreteAdjointLBM::swapP()
{
	double *temp = pLagrange_new;
	pLagrange_new = pLagrange;
	pLagrange = temp;
}

void DiscreteAdjointLBM::saveState(int n)
{

	char filename[128];
	char filename_u[128];
	char format[27];
	int ndigits = floor(log10((double)t_end) + 1.0);
	sprintf(format, "adjDataVis/%%s%%0%dd.vtk", ndigits);
	sprintf(filename, format, "adj", n);
	sprintf(filename_u, format, "adj_u", n);

	std::ofstream file(filename);
	std::ofstream file_u(filename_u);

	file << "# vtk DataFile Version 2.0\nadj\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
		 << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
		 << "\nSCALARS adj_scalars double 1\nLOOKUP_TABLE default";
	file_u << "# vtk DataFile Version 2.0\nadj_u\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
		   << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
		   << "\nVECTORS adj_dir double";

	for (int y = 0; y < NY; y++) // USED TO BE ++y  -  idk why
	{
		for (int x = 0; x < NX; x++)
		{
			if (!boundary->isOnBoundary(calcIndex(x, y)))
			{
				double res = 0.0;
				double res_uX = 0.0;
				double res_uY = 0.0;
				for (int q = 0; q < Q; q++)
				{
					res += pLagrange[calcIndex(x, y, q)];
					res_uX += (double)dir_x[q] * pLagrange[calcIndex(x, y, q)];
					res_uY += (double)dir_y[q] * pLagrange[calcIndex(x, y, q)];
				}
				file << "\n"
					 << res;
				/* if(res != 0.0)
				{
					res_uX /= res;
					res_uY /= res;
				} */
				file_u << "\n"
					   << res_uX << " " << res_uY << " 0";
			}
			else
			{
				file << "\n"
					 << "0";
				file_u << "\n"
					   << "0 0 0";
			}
		}
	}

	file.close();
	file_u.close();
}

void DiscreteAdjointLBM::iterate(int n, bool line_search, bool save_state)
{
	cost = 0.0;
	for (int y = 0; y < NY; y++)
	{
		Gx[y] = 0.0;
		Gy[y] = 0.0;
	}

	std::cout << std::endl
			  << "ITERATION #" << n << std::endl;

	if (!(n - 1))
	{
		calcRhoUxUy(n);
	}
	initP();

	if (save_state)
		saveState(0);

	int buffer = 1;
	for (int t = 1; t < t_end; t += dt)
	{
		// load data
		loadMomentum(t_end - 1 - t);

		// add to gradient
		for (int y = 0; y < NY; y++)
			addG(0, y);

		// stream pLagrange
		streamP();

		// collide
		collideP();

		// std::cout << t << std::endl;

		// +colide where we know
		if (buffer < timesSize && t_end - 1 - t == times[timesSize - 1 - buffer])
		{
			// load data when known
			loadMomentumData(t_end - 1 - t);
			buffer++;

			// costFunction();

			// collide pLagrange where we measured
			for (int counterData = 0; counterData < sizeOfKnownMomentumArray; counterData++)
			{
				// add to cost
				cost += costFunction(counterData);

				// std::cout << bufferIndex << std::endl;
				for (int q = 0; q < Q; q++)
				{
					// collide
					collidePwithData(xPos[counterData], yPos[counterData], q, counterData);
				}
			}
		}

		// swap
		// swapP();

		if (save_state && t % 10 == 0)
			saveState(t);
	}

	if (save_state)
		saveState(t_end);

	newGuess(n, line_search);
}

void DiscreteAdjointLBM::adjustGuess(double newEps, double ratio)
{
	for (int y = 0; y < NY; y++)
	{
		// guessX[y] += ratio * newEps * Gx[y];
		guessY[y] += ratio * newEps * Gy[y];
	}
}

void DiscreteAdjointLBM::plotOverLine(double *initGuessX, double *initGuessY, double *dirX, double *dirY, int n_of_steps)
{
	{
		guessX = initGuessX;
		guessY = initGuessY;
		Gx = dirX;
		Gy = dirY;
	}
	{
		std::ofstream linePlotFile("linePlot.txt");
		linePlotFile.close();
	}
	std::ofstream linePlotFile;
	linePlotFile.open("linePlot.txt", std::ofstream::app);

	for (int i = 0; i <= n_of_steps; i++)
	{
		calcRhoUxUy(10000 + i);

		saveGuessResults(10000 + i);

		std::ostringstream s;
		s << std::setprecision(10) << tempCost;

		linePlotFile << i << " " << s.str() << "\n";

		std::cout << "line step " << i << ": cost = " << s.str() << std::endl;

		adjustGuess(eps, -1.0 / (double)n_of_steps);
	}

	linePlotFile.close();
}

void DiscreteAdjointLBM::printGuess()
{
	std::cout << "guessX\t\tguessY" << std::endl;
	for (int y = 0; y < NY; y++)
		std::cout << guessX[y] << '\t' << guessY[y] << std::endl;
}

void DiscreteAdjointLBM::printNormOfGradient()
{
	std::ostringstream s;
	s << std::setprecision(10) << getNormOfGradient();
	std::cout << "Norm of gradient: " << s.str() << std::endl;
}

void DiscreteAdjointLBM::printCostFunction()
{
	std::ostringstream s;
	s << std::setprecision(10) << getCost();
	std::cout << "Cost function before guess = " << s.str() << std::endl;
}

void DiscreteAdjointLBM::costFunction()
{
	cost = 0.0;
	int buffer = 0;
	while (buffer < timesSize)
	{
		int t = times[buffer];
		buffer++;
		loadMomentum(t);
		loadMomentumData(t);
		for (int i = 0; i < sizeOfKnownMomentumArray; i++)
			cost += costFunction(i);
	}
}

double DiscreteAdjointLBM::costFunction(int momentumIndex)
{
	double result = 0.0;
	int index = calcIndex(xPos[momentumIndex], yPos[momentumIndex]);
	double rhoDiff = rho[index] - rhoDes[momentumIndex];
	double uxDiff = ux[index] - uxDes[momentumIndex];
	double uyDiff = uy[index] - uyDes[momentumIndex];
	result = 0.5 * (rhoDiff * rhoDiff + uxDiff * uxDiff + uyDiff * uyDiff);
	return result;
}

void DiscreteAdjointLBM::tempCostFunction(int t)
{
	loadMomentum(t);
	loadMomentumData(t);
	for (int i = 0; i < sizeOfKnownMomentumArray; i++)
		tempCost += costFunction(i);
}

void DiscreteAdjointLBM::loadTimes()
{
	std::ifstream timesData("adjointData/times.txt");
	if (!timesData)
		std::cout << "times file not found" << std::endl;
	int counter = 0;
	times = (int *)malloc(sizeof(int) * counter);
	while (!timesData.eof())
	{
		counter++;

		// realloc
		times = (int *)realloc(times, sizeof(int) * counter);

		// save
		timesData >> times[counter - 1];
	}
	timesSize = counter;
	timesData.close();
}

void DiscreteAdjointLBM::loadMomentumData(int n)
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

	int counter = 0;
	// times = (  int *)malloc(sizeof(  int) * counter);
	xPos = (int *)malloc(sizeof(int) * counter);
	yPos = (int *)malloc(sizeof(int) * counter);
	rhoDes = (double *)malloc(sizeof(double) * counter);
	uxDes = (double *)malloc(sizeof(double) * counter);
	uyDes = (double *)malloc(sizeof(double) * counter);
	while (!momentumFile.eof() && counter <= NX * NY)
	{
		counter++;

		// realloc
		// times = (  int *)realloc(times, sizeof(  int) * counter);
		xPos = (int *)realloc(xPos, sizeof(int) * counter);
		yPos = (int *)realloc(yPos, sizeof(int) * counter);
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
}

void DiscreteAdjointLBM::loadMomentum(int n)
{
	char filename_adj[128];
	char format[32];
	int ndigits = floor(log10((double)t_end) + 1.0);
	sprintf(format, "adjointTempDataAll/%%s%%0%dd.txt", ndigits);
	sprintf(filename_adj, format, "adj", n);

	std::ifstream momentumFile(filename_adj);
	if (!momentumFile)
		std::cout << "temp model momentum file not found" << std::endl;

	free(rho);
	free(ux);
	free(uy);

	int counter = 0;
	rho = (double *)malloc(sizeof(double) * NX * NY);
	ux = (double *)malloc(sizeof(double) * NX * NY);
	uy = (double *)malloc(sizeof(double) * NX * NY);
	while (!momentumFile.eof() && counter <= NX * NY)
	{
		// save
		momentumFile >> rho[counter];
		momentumFile >> ux[counter];
		momentumFile >> uy[counter];

		counter++;
		// std::cout << times[counter - 1] << " " << xPos[counter - 1] << " " << yPos[counter - 1] << std::endl;
	}
	momentumFile.close();

	/* std::cout << "momentum of temp model loaded SUCCESSFULLY: total " << counter << " points loaded" << std::endl
	<< std::endl; */
}

/* void DiscreteAdjointLBM::fillZeros(double *pointer,   int size)
{
	for (  int i = 0; i < size; i++)
		pointer[i] = 0.0;
} */

void DiscreteAdjointLBM::initTempModel(LBM *tempModel)
{
	tempModel->setParams(tau);
	tempModel->setPeriodic(periodic);
	tempModel->setWetNode(WN);
	tempModel->setVelocityProfile(guessX, guessY);
	tempModel->init();
}

void DiscreteAdjointLBM::calcRhoUxUy(int n)
{
	LBM *tempModel = new LBM(NX, NY, t_end, boundary);
	initTempModel(tempModel);
	loadMomentumData(times[0]);

	free(rho);
	free(ux);
	free(uy);
	rho = (double *)malloc(sizeof(double) * NX * NY);
	ux = (double *)malloc(sizeof(double) * NX * NY);
	uy = (double *)malloc(sizeof(double) * NX * NY);

	//   int counter = 0;
	int buffer = 0;
	tempCost = 0.0;
	for (int t = 0; t < t_end; t += dt)
	{
		tempModel->iterateModel();

		// save momentum everywhere, needed for adjoint collision
		char filename_adjAll[128];
		char formatAll[32];
		int ndigitsAll = floor(log10((double)t_end) + 1.0);
		sprintf(formatAll, "adjointTempDataAll/%%s%%0%dd.txt", ndigitsAll);
		sprintf(filename_adjAll, formatAll, "adj", t);
		std::ofstream momentumDataAll(filename_adjAll);
		for (int y = 0; y < NY; y++)
		{
			for (int x = 0; x < NX; x++)
			{

				int index = calcIndex(x, y);
				/* if (!boundary->isOnBoundary(index))
				{ */
				tempModel->save_momentum_to(x, y, t, rho, ux, uy, index);
				// momentumData << x << " " << y << " ";
				momentumDataAll << "\n";
				momentumDataAll << rho[index] << " " << ux[index] << " " << uy[index];
				/* } */
			}
		}
		momentumDataAll.close();

		if (buffer < timesSize && t == times[buffer])
		{
			buffer++;
			tempCostFunction(t);
		}
	}
	tempModel->save_screenshot_for_error(n);

	delete tempModel;
}

void DiscreteAdjointLBM::collideP(int x, int y, int q)
{
	pLagrange[calcIndex(x, y, q)] = omega_l * pLagrange_new[calcIndex(x, y, q)];
	int index = calcIndex(x, y);
	for (int k = 0; k < Q; k++)
	{
		double dot = ux[index] * (double)dir_x[k] + uy[index] * (double)dir_y[k];
		double power_u = ux[index] * ux[index] + uy[index] * uy[index];

		double dfeq_drho = w[k] * (1.0 + dot * cs2inv + (dot * dot) * (cs2inv * cs2inv) / 2.0 - power_u * cs2inv / 2.0);
		double dfeq_dux = w[k] * ((double)dir_x[k] * cs2inv + dot * (double)dir_x[k] * (cs2inv * cs2inv) - ux[index] * cs2inv);
		double dfeq_duy = w[k] * ((double)dir_y[k] * cs2inv + dot * (double)dir_y[k] * (cs2inv * cs2inv) - uy[index] * cs2inv);

		pLagrange[calcIndex(x, y, q)] += omega * pLagrange_new[calcIndex(x, y, k)] * (1.0 * dfeq_drho + ((double)dir_x[q] - ux[index]) * dfeq_dux + ((double)dir_y[q] - uy[index]) * dfeq_duy);
	}
}

void DiscreteAdjointLBM::collidePwithData(int x, int y, int q, int momentumIndex)
{
	int index = calcIndex(x, y);
	if (rho[index] == 0.0)
		std::cout << x << " " << y << std::endl;
	if (!boundary->isOnBoundary(index))
	{
		pLagrange[calcIndex(x, y, q)] += (rho[index] - rhoDes[momentumIndex]) * 1.0 + (ux[index] - uxDes[momentumIndex]) * ((double)dir_x[q] * rho[index] - ux[index]) / rho[index] / rho[index] + (uy[index] - uyDes[momentumIndex]) * ((double)dir_y[q] * rho[index] - uy[index]) / rho[index] / rho[index];
	}
}

void DiscreteAdjointLBM::streamP(int x, int y, int q)
{
	int x_dir = (NX + x + dir_x[q]) % NX;
	int y_dir = (NY + y + dir_y[q]) % NY; // +dir means opposite streaming from LBM (primary problem)

	if (boundary->isOnBoundary(calcIndex(x_dir, y_dir)) && !(x == 0 && (q == 3 || q == 6 || q == 7)))
	{
		pLagrange_new[calcIndex(x, y, q)] = pLagrange[calcIndex(x, y, oppositeOf[q])];
	}
	else if (periodic)
	{
		pLagrange_new[calcIndex(x, y, q)] = pLagrange[calcIndex(x_dir, y_dir, q)];
	}
	else if (WN && x == 0 && (q == 3 || q == 6 || q == 7)) // WN EQ INLET on LS
	{
		// BOUNCE-BACK VELOCITY
		pLagrange_new[calcIndex(x, y, q)] = 0.0;
	}
	else if (x == 0 && (q == 3 || q == 6 || q == 7)) // INLET on LS
	{
		// BOUNCE-BACK VELOCITY
		pLagrange_new[calcIndex(x, y, q)] = pLagrange[calcIndex(x, y, oppositeOf[q])];
	}
	else if (x == NX - 1 && (q == 1 || q == 5 || q == 8)) // OUTLET on RS (not used if instead boundary)
	{
		// ABB
		pLagrange_new[calcIndex(x, y, q)] = -pLagrange[calcIndex(x, y, oppositeOf[q])];
	}
	else
	{
		pLagrange_new[calcIndex(x, y, q)] = pLagrange[calcIndex(x_dir, y_dir, q)];
	}

	/* if (x == NX - 1) // ABB OUTLET VELOCITY EXTRAPOLATION
	{
		for (  int j = 0; j < Q; j++)
		{
			if (j == 3 || j == 6 || j == 7)
			{
				double vx = ux[calcIndex(x, y)] + 1. / 2. * (ux[calcIndex(x, y)] - ux[calcIndex(x - 1, y)]);
				double vy = uy[calcIndex(x, y)] + 1. / 2. * (uy[calcIndex(x, y)] - uy[calcIndex(x - 1, y)]);
				double dotABB = vx * (double)dir_x[oppositeOf[j]] + vy * (double)dir_y[oppositeOf[j]];
				double bx = w[oppositeOf[j]] * rhoOutlet * (cs2inv * (cs2inv * dotABB * (double)dir_x[oppositeOf[j]] - vx));
				double by = w[oppositeOf[j]] * rhoOutlet * (cs2inv * (cs2inv * dotABB * (double)dir_y[oppositeOf[j]] - vy));

				pLagrange_new[calcIndex(x, y, q)] += pLagrange[calcIndex(x, y, j)] * (3.0 * bx * (double)dir_x[q] / rho[calcIndex(x, y)] + 3.0 * by * (double)dir_y[q] / rho[calcIndex(x, y)]);

				pLagrange_new[calcIndex(x - 1, y, q)] += pLagrange[calcIndex(x, y, j)] * (-bx * (double)dir_x[q] / rho[calcIndex(x - 1, y)] - by * (double)dir_y[q] / rho[calcIndex(x - 1, y)]);
			}
		}
	} */
}

void DiscreteAdjointLBM::addG(int x, int y)
{
	if (WN)
	{
		int index = calcIndex(x, y);
		if (!boundary->isOnBoundary(index))
		{
			for (int k = 0; k < Q; k++)
			{
				double dot = guessX[y] * (double)dir_x[k] + guessY[y] * (double)dir_y[k];
				double power_u = guessX[y] * guessX[y] + guessY[y] * guessY[y];

				double dfeq_drho_w = w[k] * (1.0 + dot * cs2inv + (dot * dot) * (cs2inv * cs2inv) / 2.0 - power_u * cs2inv / 2.0);
				double dfeq_dux = w[k] * ((double)dir_x[k] * cs2inv + dot * (double)dir_x[k] * (cs2inv * cs2inv) - guessX[y] * cs2inv);
				double dfeq_duy = w[k] * ((double)dir_y[k] * cs2inv + dot * (double)dir_y[k] * (cs2inv * cs2inv) - guessY[y] * cs2inv);

				double drho_w_dux = -1.0 / (1.0 - guessX[y]) * rho[index];

				// std::cout << pLagrange[calcIndex(x, y, q)] << std::endl;
				Gx[y] += pLagrange[calcIndex(x, y, k)] * (dfeq_drho_w * drho_w_dux + dfeq_dux);
				Gy[y] += pLagrange[calcIndex(x, y, k)] * (dfeq_duy);
			}
		}
	}
	else
	{
		int index = calcIndex(x, y);
		if (!boundary->isOnBoundary(index))
		{
			for (int q = 0; q < Q; q++)
			{
				if (q == 1 || q == 5 || q == 8)
				{
					double dN_dux = -2.0 * w[oppositeOf[q]] * rho[index] * cs2inv * (double)dir_x[oppositeOf[q]];
					double dN_duy = -2.0 * w[oppositeOf[q]] * rho[index] * cs2inv * (double)dir_y[oppositeOf[q]];
					// std::cout << pLagrange[calcIndex(x, y, q)] << std::endl;
					/* if (!((y == 1 && q == 5) || (y == NY - 2 && q == 8)))
					{ */
					Gx[y] += pLagrange[calcIndex(x, y, q)] * dN_dux;
					Gy[y] += pLagrange[calcIndex(x, y, q)] * dN_duy;
					/* } */
				}
			}
		}
	}
}

void DiscreteAdjointLBM::newGuess(int n, bool line_search)
{
	adjustGuess(eps, -1.0);

	// LINE SEARCH - move new guess 1/2 closer to old guess if new guess is worse
	double newEps = eps;
	calcRhoUxUy(n);
	if (line_search)
	{
		while (tempCost > cost && newEps > 1e-10)
		{
			newEps *= 0.5;
			std::cout << "adjusted eps, new eps = " << newEps << std::endl;

			std::ofstream paramsFile;
			paramsFile.open("data/adjointParams.txt", std::ofstream::app);
			std::ostringstream s;
			s << std::setprecision(10) << newEps;
			paramsFile << "new eps = " << s.str() << " on iteration " << n << "\n";
			paramsFile.close();

			adjustGuess(newEps, 1.0);
			calcRhoUxUy(n);
		}
		setEps(newEps);
	} //? MAYBE COMMENT TO GO OVER SMALL GRADIENTS (or use momentum methods)
}
