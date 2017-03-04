//Header file to generate a structured grid

#ifndef GRIDHEADERDEF
#define GRIDHEADERDEF

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>

class grid
{
private:
  int N; //Number of nodes in both x and y direction
  double xmax, xmin; //Upper and Lower bound in x direction
  double ymax, ymin; //Upper and lower bound in y direction
  double dx, dy;
  std::vector<double> grid_x;
  std::vector<double> grid_y;
  std::vector<std::vector<double> > u;
  bool boundary(int x, int y);

public:
  
  grid()
  {
    N = 0.0;
    xmax = 0.0; xmin = 0.0;
    ymax = 0.0; ymin = 0.0;
    dx = 0.0; dy = 0.0;
    grid_x.resize(N, 0.0);
    grid_y.resize(N, 0.0);
    u.resize(N, std::vector<double>(N, 0.0));
  }

  grid(int nodes, double x_max, double x_min, double y_max, double y_min) : N(nodes), xmax(x_max), xmin(x_min), ymax(y_max), ymin(y_min)
  {
    grid_x.resize(N, 0.0);
    grid_y.resize(N, 0.0);

    dx = (xmax - xmin) / (N - 1);
    dy = (ymax - ymin) / (N - 1);

    grid_x[0] = xmin;
    grid_y[0] = ymin;

    for (int i = 1; i < N; i++)
      {
	grid_x[i] = grid_x[i-1] + dx;
	grid_y[i] = grid_y[i-1] + dy;
      }
  }

  std::vector<double> GetXGrid()
  {
    return grid_x;
  }

  std::vector<double> GetYGrid()
  {
    return grid_y;
  }

  double GridSpacingX()
  {
    return (grid_x[1] - grid_x[0]);
  }

  double GridSpacingY()
  {
    return (grid_y[1] - grid_y[0]);
  }

  std::vector<std::vector<double> > DirichletBC()
  {
	  u.resize(N, std::vector<double>(N, 0.0));
	  for (int i = 0; i < N; i++)
	  {
		  u[i][0] = sin(M_PI*grid_x[i]);
		  u[i][N-1] = sin(M_PI*grid_x[i])*exp(-M_PI);
		  u[0][i] = 0;
		  u[N-1][i] = 0;
	  }
	  u[N-1][0] = 0;
	  return u;
  }
  

  //Creating the A matrix required in the Linear System of equations
/*  std::vector<std::vector<double> > createAMatrix()
{
	  int sz = N - 2;
	  int N_square = sz*sz;
	  std::vector<std::vector<double> >A(N_square, std::vector<double>(N_square, 0.0));
		for (int i = 0; i < sz; i++)
		{
		    for (int j = 0; j < sz; j++)
		      {
			if (i == j)
			  {
			    for (int k = 0; k < sz; k++)
			      {
				int off = i*sz;
				A[k+off][k+off] = -4;
				if (k+1 < sz)
				  {
				    A[k+off][k+1+off]=1;
				  }
				if (k-1 >=0)
				  {
				    A[k+off][k-1+off]=1;
				  }
			      }
			  }
			else if (j == (i-1) || j == (i+1))
			  {
			    int off_x = i*sz, off_y = j*sz;
			    for (int k = 0; k < sz; k++)
			      {
				A[k+off_x][k+off_y] = 1;
			      }
			  }
		      }
		  }
		return A;
}
*/
  //Creating B vector which is on the RHS of the Linear System of Equations
 std::vector<double> createBVector(int N, int start, int end)
{
	  int N_square = (N-2)*(N-2);
	  std::vector<double> B(N_square, 0.0);

	  //std::vector<int> x_i(N_square, 0);
	  //std::vector<int> y_j(N_square, 0);
	  int xx, yy;
	  double result;
	  for (int i = start; i < end; i++)
	  {
		result = 0.0;
		xx = i%(N-2)+1;
		yy = i/(N-2)+1;
		//x_i[i] = xx;
		//y_j[i] = yy;
		//std::cout << xx << "\t" << yy << "\n";
		double xx_dx = double(xx)/(N-1);
		if(boundary(xx+1, yy))
		{
			result += 0.0;
		}
		if (boundary(xx-1,yy))
		{
			result += 0.0;
		}
		if(boundary(xx,yy+1))
		{
			result += sin(M_PI*xx_dx)*exp(-M_PI);
		}
		if(boundary(xx,yy-1))
		{
			result += sin(M_PI*xx_dx);
		}
		B[i-start] = -result;
	  }

	  return B;
}
  ~grid()
  {
	  std::vector<double>().swap(grid_x);
	  std::vector<double>().swap(grid_y);
	  std::vector<std::vector<double> >().swap(u);
  }


};

bool grid::boundary(int x, int y)
{
	if(x==0 || x== N-1 || y==0 || y==N-1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

#endif
    
