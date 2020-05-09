#include <iostream>
#include <cmath>
#include <fstream> // Used to create file of data

#include "ConvDiffPDE.hpp"
#include "SourceFunction3.hpp"


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
  double nu, beta, gamma;
  double length, h, error;
  double *mesh_nodes, *u_approx, *u_exact;
  int iterations;

  iterations = 5;

  nu = 0.01;
  beta = 1.0;
  gamma = 0.25;

  SourceFunction3* f3 = new SourceFunction3(nu, beta, gamma); // intialise function object

  //Create error results file
  std::ofstream file;
  file.open("Q3_errors.csv");
  assert(file.is_open());
  file << "h," << "error," << "log(h)," << "log(error)" << "\n";


  for(int k=0; k<iterations; k++)
  {
    n = 3 + int(pow(2,k)); // no. of nodes in square mesh along each direction, note: minimum 3 nodes
    length = 1.0; // square domain (0,1) x (0,1)

    ConvDiffPDE* pde = new ConvDiffPDE(n, length, *f3); // initialise PDE object

    (*pde).SolvePDE();

    h = (*pde).GetH();
    error = (*pde). GetErrorNorm();

    // save result to file
    file << h << "," << error << "," << log(h) << "," << log(error) << "\n";

    std::cout<< "n = " << n << "\n";
    std::cout<< " h = " << h << "\n";
    std::cout << " Error = " << error << "\n\n";

    if (n==5) // mesh nodes at y=0.25, y=0.5 and y=0.75 to plot
    {
      mesh_nodes = new double[n];
      u_approx = new double[int(pow(n,2))];
      u_exact = new double[int(pow(n,2))];

      // //TODO REMOVE AND REINSTATE ABOVE
      // u_approx = new double[int(pow(n-2,2))];
      // u_exact = new double[int(pow(n-2,2))];

      mesh_nodes = (*pde).GetMesh();

      u_approx = (*pde).GetUapprox();
      u_exact = (*pde).GetUexact();

      // //TODO REMOVE AND REINSTATE ABOVE
      // u_approx = (*pde).GetUvec_approx();
      // u_exact = (*pde).GetUvec_exact();

      //TODO REMOVE PRINTS
      // std::cout<< "\nMesh:";
      // printvector(n, mesh_nodes);
      // std::cout<< "\nU_approx:";
      // // printvector(int(pow(n,2)), u_approx);
      // printvector(int(pow(n-2,2)), u_approx);
      // std::cout<< "\nU_exact:";
      // // printvector(int(pow(n,2)), u_exact);
      // printvector(int(pow(n-2,2)), u_exact);
    }


    delete pde;
  }

  // close errors file
  file.close();
  std::string command = "mv Q3_errors.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q3";
  system(command.c_str());



  //Create results file for plots when n=5
  std::ofstream file2;
  file2.open("Q3_function.csv");
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
  std::string command2 = "mv Q3_function.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q3";
  system(command2.c_str());




  delete f3;
  delete[] mesh_nodes;
  delete[] u_approx;
  delete[] u_exact;

  return 0;
}
