#include <iostream>

#include "PoissonPDE.hpp"
#include "SourceFunction2.hpp"

void printvector(int n, double* vector)
{
  std::cout << "\n";

  for(int i=0; i<n; i++)
  {
    std::cout << vector[i] << "\n";
  }
}



int main(int argc, char* argv[])
{
  int n;
  double length, h, error;
  double* mesh_widths, errors, mesh_nodes, u_approx, u_exact;
  int iterations;

  iterations = 15;
  mesh_widths = new double[iterations];
  errors = new double[iterations];

  SourceFunction2* f2 = new SourceFunction2(); // intialise function object


  // TODO add a loop for n values
  n = 5; // no. of nodes in square mesh along each direction
  length = 1.0; // square domain (0,1) x (0,1)

  PoissonPDE* pde = new PoissonPDE(n, length, *f2); // initialise PDE object

  (*pde).SolvePDE();

  h = (*pde).GetH;
  error = (*pde). GetErrorNorm;



  if (n==5) // mesh nodes at y=0.25, y=0.5 and y=0.75 to plot
  {
    mesh_nodes = new double[n];
    u_approx = new double[pow(n,2)];
    u_exact = new double[pow(n,2)];

    mesh_nodes = (*pde).GetMesh;
    u_approx = (*pde).GetUapprox;
    u_exact = (*pde).GetUexact;
  }

  printvector(n, mesh_nodes);

  printvector(pow(n,2), u_approx);

  printvector(pow(n,2), u_exact);


  delete f2;
  delete pde;
  delete[] mesh_widths;
  delete[] errors;
  delete[] mesh_nodes;
  delete[] u_approx;
  delete[] u_exact;

  return 0;
}
