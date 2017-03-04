/*
 * cgkrylov.hpp
 *
 *  Created on: Dec 2, 2016
 *      Author: siddhant
 */

//Defining a class that solves the linear system of equations using the Conjugate Gradient Krylov solver
//Jacobi preconditioning is implemented

#ifndef CGKRYLOVHEADERDEF
#define CGKRYLOVHEADERDEF

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <mpi.h>


class cgkrylov
{
private:
	int N_square, nprocs, my_rank;
	//std::vector<std::vector<double> > A;
	std::vector<double> B;
	std::vector<double> x;
	std::vector<double> x_sol;
	/*std::vector<double> r;
	std::vector<double> res;

	std::vector<double> p;
	std::vector<double> z;*/
	std::vector<double> A_nnz;
	std::vector<int> col_index;
	std::vector<int> row_pointer;
	double alpha, beta, alpha_den, alpha_num, beta_num, chunk, start;
	double norm_vec(std::vector<double> &vec, int chunk);  	//Function to calculate the norm of a vector
	void sparse_mat(int start, int end, int chunk, std::vector<double> &A_nnz, std::vector<int> &col_index, std::vector<int> &row_pointer);
	void sparse_preconditioner(int start, int end, int chunk, std::vector<double> &M_nnz, std::vector<int> &M_col_ind, std::vector<int> &M_row_ptr);
	std::vector<double> sparse_mat_vec(int Nrow, int nprocs, int my_rank, std::vector<double> &A_nnz, std::vector<int> &col_ind, std::vector<int> &row_ptr, std::vector<double> &vec);
	double DotProduct(std::vector<double> &r, std::vector<double> &p, int chunk);


public:
	cgkrylov()
	{
		N_square = 0.0;
		nprocs = 0;
		my_rank = 0;
		chunk = 0.0;
		start = 0.0;
		//A.resize(N_square, std::vector<double>(N_square, 0.0));
		B.resize(chunk, 0.0);
		x.resize(chunk, 0.0);
		x_sol.resize(N_square, 0.0);
		alpha = 0.0;
		beta = 0.0;
		alpha_den = 0.0;
		alpha_num = 0.0;
		beta_num = 0.0;
		A_nnz.resize(chunk*5, 0.0);
		col_index.resize(chunk*5, 0.0);
		row_pointer.resize(chunk+1, 0.0);
	}

	cgkrylov(int size, int start, int nprocs, int chunk, int my_rank, /*std::vector<std::vector<double> > &A_cg,*/ std::vector<double> &B_cg) : N_square(size), start(start), nprocs(nprocs), chunk(chunk), my_rank(my_rank)
	{

		//A.resize(N_square, std::vector<double>(N_square, 0.0));
		B.resize(chunk, 0.0);
		x.resize(chunk, 0.0);
		x_sol.resize(N_square, 0.0);
		std::vector<double> r(chunk, 0.0);
		std::vector<double> res(chunk, 0.0);
		std::vector<double> p(chunk, 0.0);
		std::vector<double> z(chunk, 0.0);
		std::vector<double> M_nnz(chunk, 0.0);
		std::vector<int> M_col_ind(chunk, 0.0);
		std::vector<int> M_row_ptr(chunk+1, 0.0);
		//A = A_cg;
		B = B_cg;

		//Initialization
		//Computing the residual
		A_nnz.resize(chunk*5, 0.0);
		col_index.resize(chunk*5, 0.0);
		row_pointer.resize(chunk+1, 0.0);
		sparse_mat(start, start+chunk, chunk, A_nnz, col_index, row_pointer);

		res = sparse_mat_vec(N_square, nprocs, my_rank, A_nnz, col_index, row_pointer, x);

		for (int i = 0; i < chunk; i++)
		{
			r[i] = B[i] - res[i];
		}

		sparse_preconditioner(start, start+chunk, chunk, M_nnz, M_col_ind, M_row_ptr);


		//Now z = (M-inverse)*r = p
		z = sparse_mat_vec(N_square, nprocs, my_rank, M_nnz, M_col_ind, M_row_ptr, r);

		for (int i = 0; i < chunk; i++)
		{
			p[i] = z[i];
		}

		//Iterating
		//while ( norm_vec(r, chunk)/norm_vec(B, chunk) > 0.00001)				//Convergence criteria
		for (int k = 0; k < 2000; k++)
		{
			//k++;
			//std::cout << k << "\n";
			//Computing alpha
			for (int i = 0; i < chunk; i++)
			{
				res[i] = 0;
			}
			res = sparse_mat_vec(N_square, nprocs, my_rank, A_nnz, col_index, row_pointer, p);

			alpha_num = 0;
			alpha_den = 0;
			alpha = 0.0;
			beta = 0.0;
			//Computing the numerator and denominator of alpha
			alpha_num += DotProduct(r, z, chunk);

			alpha_den += DotProduct(p, res, chunk);

				alpha = alpha_num/alpha_den;

			//Updating the iterate
			for (int i = 0; i < chunk; i++)
			{
				x[i] = x[i] + alpha*p[i];
				r[i] = r[i] - alpha*res[i];
			}

			//Updating the z
			for (int i = 0; i < chunk; i++)
			{
				z[i] = 0;
			}

			z = sparse_mat_vec(N_square, nprocs, my_rank, M_nnz, M_col_ind, M_row_ptr, r);

			beta_num = 0;
			//Computing the numerator of beta

			beta_num += DotProduct(z, r, chunk);

			beta = beta_num/alpha_num;

			//Updating the values of p
			for (int i = 0; i < chunk; i++)
			{
				p[i] = z[i] + beta*p[i];
			}

			for (int i = 0; i < chunk; i++)
			{
				if (abs(r[i]) < 0.000001)
					break;
			}
			/*if (my_rank == 0)
			{
				if (norm_vec(r,chunk)/norm_vec(B,chunk) > 0.000001)
					break;
			}*/


		}


		MPI_Gather(&x.front(), chunk, MPI_DOUBLE, &x_sol.front(), chunk, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		M_nnz.shrink_to_fit();
		M_col_ind.shrink_to_fit();
		M_row_ptr.shrink_to_fit();
		r.shrink_to_fit();
		res.shrink_to_fit();
		z.shrink_to_fit();
		p.shrink_to_fit();
	}

	std::vector<double> GetCGSolution()
	{
		return x_sol;
		//return x;
	}

	~cgkrylov()
	{
		//std::vector<std::vector<double> >().swap(A);
		std::vector<double>().swap(B);
		x_sol.shrink_to_fit();
		A_nnz.shrink_to_fit();
		col_index.shrink_to_fit();
		row_pointer.shrink_to_fit();
	}
};

//Function to compute the second norm of the vector
double cgkrylov::norm_vec(std::vector<double> &vec, int chunk)
{
	double sum = 0.0, sum_global = 0.0;
	for (int i = 0; i < chunk; i++)
	{
		sum = sum + vec[i]*vec[i];
	}
	MPI_Reduce(&sum, &sum_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return sqrt(sum_global);
}

//Function to store the sparse matrix in Compressed Row Storage format, and computing the matrix-vector multiplication
//Generate sparse matrix
void cgkrylov::sparse_mat(int start, int end, int chunk, std::vector<double> &A_nnz, std::vector<int> &col_index, std::vector<int> &row_pointer)
{
	int N = sqrt(N_square);
	unsigned int nnz = 0;
	A_nnz.resize(chunk*5, 0.0);
	col_index.resize(chunk*5, 0.0);
	row_pointer.resize(chunk + 1, 0.0);
	row_pointer[0] = 0;
	for (int i = start; i < (start+chunk); i++)
	{
		if ((i - N) >= 0)
		{
			A_nnz[nnz] = 1.0;
			col_index[nnz++] = i - N;
		}
		if ((i - 1) >= 0 && i%N != 0)
		{
			A_nnz[nnz] = 1.0;
			col_index[nnz++] = i - 1;
		}
		A_nnz[nnz] = -4.0;
		col_index[nnz++] = i;
		if ((i + 1) >= 0 && (i+1)%N != 0)
		{
			A_nnz[nnz] = 1.0;
			col_index[nnz++] = i + 1;
		}
		if ((i+N) < N_square)
		{
			A_nnz[nnz] = 1.0;
			col_index[nnz++] = i + N;
		}
	row_pointer[i + 1 - start] = nnz;
	}
}

void cgkrylov::sparse_preconditioner(int start, int end, int chunk, std::vector<double> &M_nnz, std::vector<int> &M_col_ind, std::vector<int> &M_row_ptr)
{
	int nnz = 0;
	M_nnz.resize(chunk, 0.0);
	M_col_ind.resize(chunk, 0.0);
	M_row_ptr.resize(chunk+1, 0.0);
	M_row_ptr[0] = 0;
	for (int i = start; i < start + chunk; i++)
	{
		M_nnz[i - start] = 0.25;
		M_col_ind[i - start] = i;
		M_row_ptr[i - start + 1] = ++nnz;
	}
}

std::vector<double> cgkrylov::sparse_mat_vec(int Nrow, int nprocs, int my_rank, std::vector<double> &mat_nnz, std::vector<int> &col_ind, std::vector<int> &row_ptr, std::vector<double> &vec)
{
	std::vector<double> B_vec(Nrow, 0.0);
	std::vector<double> product((Nrow/nprocs), 0.0);
	double s;
	int rc[nprocs], dis[nprocs];

	for (int i = 0, start = 0; i < nprocs; i++)
	{
		dis[i] = start;
		rc[i] = Nrow/nprocs;
		if (i < Nrow%nprocs)

			rc[i]++;
			start += rc[i];

	}
	MPI_Allgatherv(&vec[0], rc[my_rank], MPI_DOUBLE, &B_vec[0], rc, dis, MPI_DOUBLE, MPI_COMM_WORLD);


	for (int row = 0; row < rc[my_rank]; row++)
	{
		s = 0.0;
		for (int icol = row_ptr[row]; icol < row_ptr[row+1]; icol++)
		{
			int col = col_ind[icol];
			s += mat_nnz[icol]*B_vec[col];
		}
		product[row] = s;
	}
	return product;
}

double cgkrylov::DotProduct(std::vector<double> &r, std::vector<double> &p, int chunk)
{
	double prod_local = 0.0, prod_global = 0.0;
	for (int i = 0; i < chunk; i++)
	{
		prod_local += r[i]*p[i];
	}
	MPI_Allreduce(&prod_local, &prod_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return prod_global;
}


#endif /* CGKRYLOV_HPP_ */
