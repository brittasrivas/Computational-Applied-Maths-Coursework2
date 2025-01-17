#include <iostream>
#include <cmath>
#include <fstream> // Used to create file of data
#include <string>

#include "PoissonPDE.hpp"
#include "SourceFunction2.hpp"



int main(int argc, char* argv[])
{
  int n;
  double length, h, error;
  double *mesh_nodes, *u_approx, *u_exact;
  int iterations;

  iterations = 5;

  SourceFunction2* f2 = new SourceFunction2(); // intialise function object

  //Create error results file
  std::ofstream file;
  file.open("Q2_errors.csv");
  assert(file.is_open());
  file << "h," << "error," << "log(h)," << "log(error)," << "log(h^2.0)" <<  "\n";


  for(int k=0; k<iterations; k++)
  {
    n = 5 * (2*(k+1)-1); // no. of nodes in square mesh along each direction, note: minimum 3 nodes
    length = 1.0; // square domain (0,1) x (0,1)

    PoissonPDE* pde = new PoissonPDE(n, length, *f2); // initialise PDE object

    (*pde).SolvePDE();

    h = (*pde).GetH();
    error = (*pde). GetErrorNorm();

    // save result to file
    file << h << "," << error << "," << log(h) << "," << log(error) << "," << log(pow(h,2.0)) << "\n";

    std::cout<< "n = " << n << "\n";
    std::cout<< " h = " << h << "\n";
    std::cout << " Error = " << error << "\n\n";

    if (n==5) // mesh nodes at y=0.25, y=0.5 and y=0.75 to plot
    {
      mesh_nodes = new double[n];
      u_approx = new double[int(pow(n,2))];
      u_exact = new double[int(pow(n,2))];

      mesh_nodes = (*pde).GetMesh();

      u_approx = (*pde).GetUapprox();
      u_exact = (*pde).GetUexact();

    }


    delete pde;
  }

  // close errors file
  file.close();
  std::string command = "mv Q2_errors.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q2";
  system(command.c_str());



  //Create results file for plots when n=5
  std::ofstream file2;
  file2.open("Q2_function.csv");
  assert(file2.is_open());
  file2 << "mesh,"
    << "approx y=0.25," << "exact y=0.25,"
    << "approx y=0.5," << "exact y=0.5,"
    << "approx y=0.75," << "exact y=0.75"
    << "\n";

  n=5;
  for(int i=0; i<n; i++)
  {
    file2 << mesh_nodes[i]
      << "," << u_approx[n+i] << "," << u_exact[n+i]
      << "," << u_approx[2*n+i] << "," << u_exact[2*n+i]
      << "," << u_approx[3*n+i] << "," << u_exact[3*n+i]
      << "\n";
  }
  file2.close();
  std::string command2 = "mv Q2_function.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q2";
  system(command2.c_str());




  delete f2;
  delete[] mesh_nodes;
  delete[] u_approx;
  delete[] u_exact;

  return 0;
}
