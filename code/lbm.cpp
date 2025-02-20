#include "lbm.h"

LBM::LBM(int x, int y, int T, BounceBack *b)
	: NX(x), NY(y), t_end(T), boundary(b)
{
	const int sizeScalar = sizeof(double) * NX * NY;
	const int sizeVector = sizeof(double) * NX * NY * Q;
	f_old = (double *)malloc(sizeVector);
	f_new = (double *)malloc(sizeVector);
	f_eq = (double *)malloc(sizeVector);
	rho = (double *)malloc(sizeScalar);
	ux = (double *)malloc(sizeScalar);
	uy = (double *)malloc(sizeScalar);
	S = (double *)malloc(sizeVector);
	// srand(time(NULL));
}

LBM::~LBM()
{
	free(f_old);
	free(f_new);
	free(f_eq);
	free(rho);
	free(ux);
	free(uy);
	free(S);
}

void LBM::setParams(double newTau)
{
	tau = newTau;
	// tau = newViscosity * cs2inv * realDt / (realDx * realDx) + 1.0 / 2.0;
	omega = (double)dt / tau;
	omega_l = 1.0 - omega;
}
void LBM::printParams()
{
	std::cout << "NX = " << NX << std::endl
			  << "NY = " << NY << std::endl
			  << "t_end = " << t_end << std::endl;
	/* // double Reynolds = (double)NX * ;
	std::cout << "NX = " << NX << std::endl
			  << "NY = " << NY << std::endl
			  << "t_end = " << t_end << std::endl
			  << "tau = " << tau << std::endl;
	// << "Re = " << Reynolds; */
}

void LBM::setC_u_rho(double realDx, double realDt, double realRho)
{
	C_u = realDx / realDt;
	std::cout << "C_u = " << C_u << std::endl;
	C_rho = realRho;
}

void LBM::set_noise_gen(double noiseStdDev)
{
	gen.seed(rd());
	normal_d.param(std::normal_distribution<>::param_type{1.0, noiseStdDev});
}

void LBM::initTimesFile()
{
	std::ofstream timesData;
	timesData.open("adjointData/times.txt");
	timesData.close();
}

void LBM::setForce(double Fx, double Fy)
{
	forceX = Fx;
	forceY = Fy;
}

void LBM::init()
{
	for (int x = 0; x < NX; x++)
	{
		for (int y = 0; y < NY; y++)
		{
			rho[calcIndex(x, y)] = rho0;
			ux[calcIndex(x, y)] = 0.0 - 0.5 * rho0 * forceX * (double)dt / rho0;
			uy[calcIndex(x, y)] = 0.0 - 0.5 * rho0 * forceY * (double)dt / rho0;

			/* if (WN && x == 0)
			{
				ux[calcIndex(x, y)] = velocityX[y]; // velocityX[y]
				uy[calcIndex(x, y)] = velocityY[y]; // velocityY[y]
				rho[calcIndex(x, y)] *= c / (c - velocityX[y]);
			} */

			if (boundary->isOnBoundary(calcIndex(x, y)))
			{
				rho[calcIndex(x, y)] = 0.0;
				ux[calcIndex(x, y)] = 0.0;
				uy[calcIndex(x, y)] = 0.0;
			}

			double power_u = ux[calcIndex(x, y)] * ux[calcIndex(x, y)] + uy[calcIndex(x, y)] * uy[calcIndex(x, y)];
			for (int k = 0; k < Q; k++)
			{
				double dot = ux[calcIndex(x, y)] * (double)dir_x[k] + uy[calcIndex(x, y)] * (double)dir_y[k];
				f_old[calcIndex(x, y, k)] = w[k] * rho[calcIndex(x, y)] * (1.0 + dot * cs2inv + (dot * dot) * (cs2inv * cs2inv) / 2.0 - power_u * cs2inv / 2.0);
				f_new[calcIndex(x, y, k)] = 0.0;
			}
		}
	}
}

void LBM::updMacro(int x, int y)
{
	int index = calcIndex(x, y);
	rho[index] = 0.;
	ux[index] = 0.;
	uy[index] = 0.;
	if (!(boundary->isOnBoundary(index)))
	{
		for (int k = 0; k < Q; k++)
		{
			rho[index] += f_old[calcIndex(x, y, k)];
			ux[index] += (double)dir_x[k] * f_old[calcIndex(x, y, k)];
			uy[index] += (double)dir_y[k] * f_old[calcIndex(x, y, k)];
		}
		ux[index] /= rho[index];
		ux[index] += 0.5 * forceX * (double)dt;
		uy[index] /= rho[index];
		uy[index] += 0.5 * forceY * (double)dt;

		/* if (WN && x == 0)
		{
			rho[index] = c / (c - velocityX[y]) * (f_old[calcIndex(x, y, 0)] + f_old[calcIndex(x, y, 2)] + f_old[calcIndex(x, y, 4)] + 2.0 * (f_old[calcIndex(x, y, 3)] + f_old[calcIndex(x, y, 6)] + f_old[calcIndex(x, y, 7)]) - 0.5 * forceX * rho[index] * (double)dt);
			ux[index] = velocityX[y];
			uy[index] = velocityY[y];
			f_old[calcIndex(x, y, 1)] += 2. / 3. * rho[index] * velocityX[y];
			f_old[calcIndex(x, y, 5)] += 1. / 2. * rho[index] * velocityY[y] + 1. / 6. * rho[index] * velocityX[y];
			f_old[calcIndex(x, y, 8)] += -1. / 2. * rho[index] * velocityY[y] + 1. / 6. * rho[index] * velocityX[y];
		} */

		/* if (!periodic && (x == 0 || x == NX - 1))
		{
			ux[calcIndex(x, y)] = velocityX[y];
			uy[calcIndex(x, y)] = velocityY[y];
		} */

		if (boundary->isOnBoundary(calcIndex(x, y)))
		{
			rho[calcIndex(x, y)] = 0.0;
			ux[calcIndex(x, y)] = 0.0;
			uy[calcIndex(x, y)] = 0.0;
		}
	}
}

void LBM::calcEq(int x, int y, int k)
{
	int index = calcIndex(x, y);
	double power_u = ux[index] * ux[index] + uy[index] * uy[index];
	double dot = ux[index] * (double)dir_x[k] + uy[index] * (double)dir_y[k];
	f_eq[calcIndex(x, y, k)] = w[k] * rho[index] * (1. + dot * cs2inv + (dot * dot) * (cs2inv * cs2inv) / 2. - power_u * cs2inv / 2.);
	if (boundary->isOnBoundary(index))
	{
		f_eq[calcIndex(x, y, k)] = 0.;
	}
}

void LBM::calcS(int x, int y, int k)
{
	int index = calcIndex(x, y);
	double dot = ux[index] * (double)dir_x[k] + uy[index] * (double)dir_y[k];
	double Fx = rho[index] * forceX;
	double Fy = rho[index] * forceY;
	double dotF = (double)dir_x[k] * Fx + (double)dir_y[k] * Fy;
	double F_k = w[k] * cs2inv * (dotF + cs2inv * dotF * dot - ux[index] * Fx - uy[index] * Fy);
	S[calcIndex(x, y, k)] = (1.0 - 0.5 * (double)dt / tau) * F_k;
}

void LBM::collide(int x, int y, int k)
{
	f_old[calcIndex(x, y, k)] = omega_l * f_old[calcIndex(x, y, k)] + omega * f_eq[calcIndex(x, y, k)] + (double)dt * S[calcIndex(x, y, k)];
	if (boundary->isOnBoundary(calcIndex(x, y)))
	{
		f_old[calcIndex(x, y, k)] = 0.;
	}
}

void LBM::streamOn(int x, int y, int k)
{
	int x_dir = (NX + x - dir_x[k]) % NX;
	int y_dir = (NY + y - dir_y[k]) % NY;

	if (boundary->isOnBoundary(calcIndex(x_dir, y_dir)) && !(x == 0 && (k == 1 || k == 5 || k == 8)))
	{
		f_new[calcIndex(x, y, k)] = f_old[calcIndex(x, y, oppositeOf[k])];
	}
	else if (periodic)
	{
		f_new[calcIndex(x, y, k)] = f_old[calcIndex(x_dir, y_dir, k)];
	}
	else if (WN && x == 0)
	{
		f_new[calcIndex(x, y, k)] = f_old[calcIndex(x_dir, y_dir, k)];
		/* double diff = 0.5 * (f_old[calcIndex(x, y, 2)] - f_old[calcIndex(x, y, 4)]);
		f_new[calcIndex(x, y, k)] = f_old[calcIndex(x, y, oppositeOf[k])] + (double)dir_y[oppositeOf[k]] * diff; */
	}
	else if (x == 0 && (k == 1 || k == 5 || k == 8)) // INLET on LS
	{
		// BOUNCE-BACK VELOCITY
		double dot = velocityX[y] * (double)dir_x[oppositeOf[k]] + velocityY[y] * (double)dir_y[oppositeOf[k]];
		f_new[calcIndex(x, y, k)] = f_old[calcIndex(x, y, oppositeOf[k])] - 2 * w[oppositeOf[k]] * rho[calcIndex(x, y)] * cs2inv * dot;
	}
	else if (x == NX - 1 && (k == 3 || k == 6 || k == 7)) // OUTLET on RS
	{
		// ANTI-BOUNCE-BACK EXTRAPOLATED VELOCITY
		double vx = ux[calcIndex(x, y)] + 1. / 2. * (ux[calcIndex(x, y)] - ux[calcIndex(x - 1, y)]);
		double vy = uy[calcIndex(x, y)] + 1. / 2. * (uy[calcIndex(x, y)] - uy[calcIndex(x - 1, y)]);
		double dotABB = vx * (double)dir_x[oppositeOf[k]] + vy * (double)dir_y[oppositeOf[k]];
		double power_u = vx * vx + vy * vy;

		f_new[calcIndex(x, y, k)] = -f_old[calcIndex(x, y, oppositeOf[k])] + 2. * w[oppositeOf[k]] * rhoOutlet * (1. + dotABB * dotABB * cs2inv * cs2inv / 2. - power_u * cs2inv / 2.);
	}
	else
	{
		f_new[calcIndex(x, y, k)] = f_old[calcIndex(x_dir, y_dir, k)];
	}
}

double LBM::random_noise()
{
	return normal_d(gen);
}

void LBM::streamInOut() // post-stream WN operations
{
	for (int y = 0; y < NY; y++)
	{
		int x = 0;
		if (!boundary->isOnBoundary(calcIndex(x, y)))
		{
			// calc rho on wall
			double rho_w = c / (c - velocityX[y]) * (f_new[calcIndex(x, y, 0)] + f_new[calcIndex(x, y, 2)] + f_new[calcIndex(x, y, 4)] + 2.0 * (f_new[calcIndex(x, y, 3)] + f_new[calcIndex(x, y, 6)] + f_new[calcIndex(x, y, 7)])) /* - 0.5 * forceX * rho[index] * (double)dt) */;
			double power_u = velocityX[y] * velocityX[y] + velocityY[y] * velocityY[y];

			for (int k = 0; k < Q; k++)
			{
				double dot = velocityX[y] * (double)dir_x[k] + velocityY[y] * (double)dir_y[k];

				// set post streaming f to eq
				f_new[calcIndex(x, y, k)] = w[k] * rho_w * (1. + dot * cs2inv + (dot * dot) * (cs2inv * cs2inv) / 2. - power_u * cs2inv / 2.);
			}
		}
	}
}

void LBM::update_macroscopic()
{
	for (int x = 0; x < NX; x++)
	{
		for (int y = 0; y < NY; y++)
		{
			updMacro(x, y);
		}
	}
}

void LBM::stream()
{
	for (int x = 0; x < NX; x++)
	{
		for (int y = 0; y < NY; y++)
		{
			for (int k = 0; k < Q; k++)
			{
				if (!(boundary->isOnBoundary(calcIndex(x, y))))
				{
					streamOn(x, y, k);
				}
			}
		}
	}
	if (WN) // WN EQ INLET on LS
	{
		streamInOut();
	}
}

void LBM::swap()
{
	double *temp = f_new;
	f_new = f_old;
	f_old = temp;
}

void LBM::save_macroscopic(int n)
{
	char filename_rho[128];
	char filename_u[128];
	char filename_mm[128];
	char format[21];
	int ndigits = floor(log10((double)t_end) + 1.0);
	sprintf(format, "data/%%s%%0%dd.vtk", ndigits);
	sprintf(filename_rho, format, "ro", n);
	sprintf(filename_u, format, "u", n);
	sprintf(filename_mm, format, "mm", n);

	std::ofstream file_rho(filename_rho);
	std::ofstream file_u(filename_u);
	std::ofstream file_mm(filename_mm);

	file_rho << "# vtk DataFile Version 2.0\ndensity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
			 << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
			 << "\nSCALARS density_scalars double 1\nLOOKUP_TABLE default";
	file_u << "# vtk DataFile Version 2.0\nvelocity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
		   << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
		   << "\nVECTORS velocity_vectors double";
	file_mm << "# vtk DataFile Version 2.0\nmomentum\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
			<< NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
			<< "\nVECTORS momentum_vectors double";

	for (int y = 0; y < NY; y++) // USED TO BE ++y  -  idk why
	{
		for (int x = 0; x < NX; x++)
		{
			if (!boundary->isOnBoundary(calcIndex(x, y)))
			{
				file_rho << "\n"
						 << (C_rho * rho[calcIndex(x, y)]);
				file_u << "\n"
					   << (C_u * ux[calcIndex(x, y)]) << " " << (C_u * uy[calcIndex(x, y)]) << " 0";
				file_mm << "\n"
						<< (C_rho * rho[calcIndex(x, y)]) * (C_u * ux[calcIndex(x, y)]) << " " << (C_rho * rho[calcIndex(x, y)]) * (C_u * uy[calcIndex(x, y)]) << " 0";
			}
			else
			{
				file_rho << "\n"
						 << "0";
				file_u << "\n"
					   << "0 0 0";
				file_mm << "\n"
						<< "0 0 0";
			}
		}
	}

	file_rho.close();
	file_u.close();
	file_mm.close();
}

void LBM::upd_calc_coll()
{
	// #pragma omp parallel for collapse(2)
	for (int x = 0; x < NX; x++)
	{
		for (int y = 0; y < NY; y++)
		{
			if (!(boundary->isOnBoundary(calcIndex(x, y))))
			{
				// update_macroscopic
				updMacro(x, y);
				for (int k = 0; k < Q; k++)
				{
					// calc_eq
					calcEq(x, y, k);

					// calc force terms
					calcS(x, y, k);

					// collide
					collide(x, y, k);
				}
			}
		}
	}
}

void LBM::save_momentum(int dx, int dy, int n)
{
	char filename_adj[128];
	char format[28];
	int ndigits = floor(log10((double)t_end) + 1.0);
	sprintf(format, "adjointData/%%s%%0%dd.txt", ndigits);
	sprintf(filename_adj, format, "adj", n);

	std::ofstream momentumData(filename_adj);

	std::ofstream timesData("adjointData/times.txt", std::ofstream::app);

	timesData << '\n'
			  << n;

	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			if (!boundary->isOnBoundary(calcIndex(x, y)))
			{
				if (x % dx == 0 && y % dy == 0)
				{
					momentumData << "\n";
					momentumData << x << " " << y << " ";
					momentumData << rho[calcIndex(x, y)] << " " << ux[calcIndex(x, y)] << " " << uy[calcIndex(x, y)];
				}
			}
		}
	}

	momentumData.close();
	timesData.close();
}

void LBM::save_momentum_to(int x, int y, int n, double *rhoTrg, double *uxTrg, double *uyTrg, int index)
{
	rhoTrg[index] = rho[calcIndex(x, y)];
	uxTrg[index] = ux[calcIndex(x, y)];
	uyTrg[index] = uy[calcIndex(x, y)];
}

void LBM::runModel()
{
	init();
	for (int t = 0; t < t_end; t += dt)
	{
		upd_calc_coll();
		stream();
		swap();
		/* if (!(t % SAVE) && t >= START)
		{
			model.save_macroscopic(t);
		} */
	}
	// model.save_macroscopic(T);
}

void LBM::iterateModel()
{
	upd_calc_coll();
	stream();
	swap();
}

void LBM::save_velocityProfile()
{
	std::ofstream file("velocityProfile.txt");

	for (int y = 0; y < NY; y++)
	{
		file << y << " " << velocityX[y] /* << " " << guessY[y - 1] */ << "\n";
	}

	file.close();
}

void LBM::save_velocity_at(int x)
{
	std::ofstream file("velocityProfileAtX.txt");

	for (int y = 1; y < NY - 1; ++y)
	{
		file << y << " " << (C_u * ux[calcIndex(x, y)]) /* << " " << (C_u * uy[calcIndex(x, y)]) */ << "\n";
	}

	file.close();
}

void LBM::save_screenshot_for_error(int n)
{
	char filename_u[128];
	char filename_rho[128];
	char format[26];
	int ndigits = floor(log10((double)100) + 1.0);
	sprintf(format, "errorData/%%s%%0%dd.vtk", ndigits);
	sprintf(filename_u, format, "u", n);
	sprintf(filename_rho, format, "ro", n);

	std::ofstream file_rho(filename_rho);
	std::ofstream file_u(filename_u);

	file_rho << "# vtk DataFile Version 2.0\ndensity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
			 << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
			 << "\nSCALARS density_scalars double 1\nLOOKUP_TABLE default";
	file_u << "# vtk DataFile Version 2.0\nvelocity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
		   << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
		   << "\nVECTORS velocity_vectors double";

	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			if (!boundary->isOnBoundary(calcIndex(x, y)))
			{
				file_rho << "\n"
						 << (C_rho * rho[calcIndex(x, y)]);
				file_u << "\n"
					   << (C_u * ux[calcIndex(x, y)]) << " " << (C_u * uy[calcIndex(x, y)]) << " 0";
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
}

void LBM::save_momentum_noise(int dx, int dy, int n, int save_every, double forget)
{
	char filename_adj[128];
	char format[28];
	int ndigits = floor(log10((double)t_end) + 1.0);
	sprintf(format, "adjointData/%%s%%0%dd.txt", ndigits);
	sprintf(filename_adj, format, "adj", n);

	std::ofstream file_rho;
	std::ofstream file_u;

	if (!(n % save_every) || n == t_end - 1)
	{
		char filename_rho[128];
		char filename_u[128];
		char formatVTK[26];
		sprintf(formatVTK, "noiseData/%%s%%0%dd.vtk", ndigits);
		sprintf(filename_rho, formatVTK, "ro_m", n);
		sprintf(filename_u, formatVTK, "u_m", n);

		file_rho.open(filename_rho);
		file_u.open(filename_u);

		file_rho << "# vtk DataFile Version 2.0\ndensity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
				 << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
				 << "\nSCALARS density_scalars double 1\nLOOKUP_TABLE default";
		file_u << "# vtk DataFile Version 2.0\nvelocity\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS "
			   << NX << " " << NY << " 1\nORIGIN 0 0 0\nSPACING 1 1 0\nPOINT_DATA " << NX * NY
			   << "\nVECTORS velocity_vectors double";
	}

	std::ofstream momentumData(filename_adj);

	std::ofstream timesData("adjointData/times.txt", std::ofstream::app);

	timesData << '\n'
			  << n;

	for (int y = 0; y < NY; y++)
	{
		for (int x = 0; x < NX; x++)
		{
			if (!boundary->isOnBoundary(calcIndex(x, y)))
			{
				if (x % dx == 0 && y % dy == 0 && x >= (int)(forget * NX) /* && x <= 1 * NX / 2 */)
				{

					momentumData << "\n";
					momentumData << x << " " << y << " ";

					double noise = random_noise();

					// double noise = 1.0 + (((double)rand() * (2.0 * noisePower) / (double)RAND_MAX) - noisePower);
					if (!(n % save_every) || n == t_end - 1)
						file_rho << "\n"
								 << rho[calcIndex(x, y)] * noise;
					momentumData << rho[calcIndex(x, y)] * noise << " ";

					noise = random_noise();
					// noise = 1.0 + (((double)rand() * (2.0 * noisePower) / (double)RAND_MAX) - noisePower);
					if (!(n % save_every) || n == t_end - 1)
						file_u << "\n"
							   << ux[calcIndex(x, y)] * noise << " ";
					momentumData << ux[calcIndex(x, y)] * noise << " ";

					noise = random_noise();
					// noise = 1.0 + (((double)rand() * (2.0 * noisePower) / (double)RAND_MAX) - noisePower);
					file_u << uy[calcIndex(x, y)] * noise << " 0";
					momentumData << uy[calcIndex(x, y)] * noise;
				}
				else if (!(n % save_every) || n == t_end - 1)
				{
					file_rho << "\n"
							 << "0";
					file_u << "\n"
						   << "0 0 0";
				}
			}
			else if (!(n % save_every) || n == t_end - 1)
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
	momentumData.close();
	timesData.close();
}
