/*
 * parallel.cpp
 *
 *  Created on: Dec 6, 2016
 *      Author: siddhant
 */

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cassert>
#include "grid.hpp"
#include "cgkrylov.hpp"
#include <cstdlib>

#define N 290

#define root 0

bool boundary(int x, int y);

int main(int argc, char* argv[])
{
	int err;
	int nprocs, my_rank, sz, N_square;
	int chunk = 0, start = 0;
	double time1, time2;

	sz = N - 2;
	N_square = sz*sz;


	//Initializing the grid vectors
	//std::vector<double> x;
	//std::vector<double> y;
	grid G(N, 1, 0, 1, 0);

	//x = G.GetXGrid();
	//y = G.GetYGrid();


/*************************************THE COMMENTED PART BELOW IS ONLY FOR WRITING THE ANALYTICAL SOLUTION. IT IS COMMENTED TO REDUCE THE USAGE OF MEMORY**************************************/
/*	//Write the grid data to a file
	std::ofstream write_grid("x.dat");
	assert(write_grid.is_open());
	for (int i = 0; i < N; i++)
	{
		write_grid << x[i] << "\t" << y[i] << std::endl;
	}

	//Analytical solution
	std::vector<std::vector<double> > u_sol(N, std::vector<double>(N, 0.0));

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			u_sol[i][j] = sin(M_PI*x[i])*exp(-M_PI*y[j]);
		}
	}

	//Writing the analytical solution

	std::ofstream write_usol("usol.dat");
	assert(write_usol.is_open());

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			write_usol << u_sol[i][j] << " ";
		}
		write_usol << "\n";
	}
*/
	std::vector<std::vector<double> > u(N, std::vector<double>(N, 0.0));
	u = G.DirichletBC();

	err = MPI_Init(&argc, &argv);
	err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	err = MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	time1 = MPI_Wtime();

	//Partitioning
	chunk = N_square/nprocs;
	//std::cout << "chunk = " << chunk << "\n";

	start = (N_square/nprocs)*my_rank;
	//std::cout << "start = " << start << "\n";


	if (my_rank < (N_square)%nprocs)
	{
		chunk++;
		start = start + my_rank;
	}
	else if(my_rank == (N_square)%nprocs)
	{
		start = start + my_rank;
	}
	else
	{
		start = start + N_square%nprocs;
	}

	//Initializing the A matrix and B vector
//	std::vector<std::vector<double> > A(N_square, std::vector<double>(N_square, 0.0));
	std::vector<double> B(chunk, 0.0);

	//A = G.createAMatrix();
	B = G.createBVector(N, start, (start+chunk));
	cgkrylov CG(N_square, start, nprocs, chunk, my_rank, B);

	std::vector<double> x_sol(N, 0.0);

	x_sol = CG.GetCGSolution();

	time2 = MPI_Wtime();



/******************************************THE COMMENTED PART BELOW IS TO WRITE THE NEW SOLUTION TO A FILE. IT IS COMMENTED TO AVOID THE EXCESSIVE USAGE OF MEMORY**********************************/
/*	if (my_rank == 0)
{

		std::vector<std::vector<double> > unew(N, std::vector<double>(N, 0.0));

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				unew[i][j] = u[i][j];
			}
		}

		//std::vector<double> x_i(N_square, 0.0);
		//std::vector<double> y_j(N_square, 0.0);
		double xx, yy;
		for (int i = 0; i < N_square; i++)
		{
			xx = i%(N-2)+1;
			yy = i/(N-2)+1;

			unew[xx][yy] = x_sol[i];
		}

		std::ofstream write_unew("unew.dat");
		assert(write_unew.is_open());

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				write_unew << unew[i][j] << " ";
			}
			write_unew << "\n";
		}
	}
*/
	if (my_rank == 0)
	{
		std::cout << "Running time = " << time2 - time1 << "\n";
	}

	MPI_Finalize();

//	u.shrink_to_fit();
	//u_sol.shrink_to_fit();
	//unew.shrink_to_fit();
	/*x_sol.shrink_to_fit();
//	x_i.shrink_to_fit();
//	y_j.shrink_to_fit();*/


	return 0;
}
