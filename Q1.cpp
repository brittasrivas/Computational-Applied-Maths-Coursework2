#include <iostream>
#include <cmath>
#include <fstream> // Used to create file of data
#include <string>

//TODO Function Prototypes


double* createMesh(const int n, const double h, const double x_start, const double x_end)
{
  double* mesh;
  mesh = new double[n];


  for (int j=0; j<n; j++)
  {
    mesh[j] = j*h;
  }

  return mesh;
}


double* solveTridiagonal(const int m, double *d, double *u, double *l, double *Fvec)
// Solve Ax=b for A tridiagonal represented by diagonal vectors d, u and l.
// Algorithm from Epperson 2013, Section 2.6.
{
  double* x;
  x = new double[m];

  // Elimination stage
  for(int i=1; i<m; i++)
  {
    d[i] = d[i] - u[i-1]*(l[i-1]/d[i-1]);
    Fvec[i] = Fvec[i] - Fvec[i-1]*(l[i-1]/d[i-1]);
  }

  //Backsolve
  x[m-1] = Fvec[m-1]/d[m-1];
  for(int i=m-2; i>=0; i--)
  {
    x[i] = ( Fvec[i] - u[i]*x[i+1] )/d[i];
  }

  return x;
}


void testCoeffcientMatrix(const int m, const double h, double* A_d, double* A_u, double* A_l)
{
  for (int i=0; i<m-1; i++)
  {
    A_d[i] = 2.0/pow(h,2.0);
    A_u[i] = -1.0/pow(h,2.0);
    A_l[i] = -1.0/pow(h,2.0);
  }

  A_d[m-1] = 2.0/pow(h,2.0);


}


double testExactU(double x)
{
  return x - sin(M_PI*x);
}

void outputExactUvec(const int n, const double* mesh, double* Uvec_exact)
{
  for(int i=0; i<n; i++)
  {
    Uvec_exact[i] = testExactU(mesh[i]);
  }
}


double testExactF(double x)
{
  return (-1 * pow(M_PI,2.0) * sin(M_PI*x));
}


void testFvec(const int m, const double alpha, const double beta, const double h, const double* mesh, double* Fvec)
{
   for (int i=0; i<m; i++)
   {
     Fvec[i] = testExactF(mesh[i+1]);
     //Note use mesh[i+1] since mesh is of length n, whilst Fvec is of length m (i.e. excludes first and last nodes in mesh)
   }

   Fvec[0] += alpha/(pow(h,2.0));
   Fvec[m-1] += beta/(pow(h,2.0));
}


void calculateGridFunctionError(const int m, const double h, const double* Uvec_full, const double* Uvec_exact, double &grid_function_error)
{
  double sum=0.0;

  for(int i=1; i<=m; i++)
  {
    sum += pow((Uvec_exact[i] - Uvec_full[i]),2.0);
  }

  grid_function_error = pow((h * sum), 0.5);
}


double executeQuestion(const int n, const double x_start, const double x_end, double &h, double* Uvec_full, double* Uvec_exact)
{
  int m = n-2;
  h = (x_end - x_start)/double(m+1);
  double alpha = 0.0;
  double beta = 1.0;

  // Mesh
  double *mesh;
  mesh = new double[n]; //mesh of length n
  mesh = createMesh(n, h, x_start, x_end);

  // F vector
  double *Fvec;
  Fvec = new double[m];

  // Construct Question1 F vector
  testFvec(m, alpha, beta, h, mesh, Fvec);

  // Coefficient matrix A
  double *A_d, *A_u, *A_l; //tridiagonal matrix components
  A_d = new double[m]; //diagonal
  A_u = new double[m-1]; //upper diagonal
  A_l = new double[m-1]; //lower diagonal

  // Construct Question1 Coefficient matrix A
  testCoeffcientMatrix(m, h, A_d, A_u, A_l);

  // U approx
  double *Uvec; // holds positions 1 to m of the U vector
  Uvec = new double[m];
  Uvec = solveTridiagonal(m, A_d, A_u, A_l, Fvec);

  // U approx, including boundaries
  Uvec_full[0] = alpha;
  Uvec_full[n-1] = beta;

  for(int j=1; j<=m; j++)
  {
    Uvec_full[j] = Uvec[j-1];
  }


  // U exact, including boundaries
  outputExactUvec(n, mesh, Uvec_exact);


  // Error
  double grid_function_error;
  calculateGridFunctionError(m, h, Uvec_full, Uvec_exact, grid_function_error);

  delete[] mesh;
  delete[] A_d;
  delete[] A_u;
  delete[] A_l;
  delete[] Fvec;
  delete[] Uvec;

  return grid_function_error;
}


int main(int argc, char* argv[])
{

  int n;
  double h;
  double x_start;
  double x_end;
  double grid_function_error;

  int maxIterations = 15;

  double *errors, *mesh_widths;
  errors = new double[maxIterations];
  mesh_widths = new double[maxIterations];

  for(int k=0; k<maxIterations; k++)
  {

    n = 3 + int(pow(2,k)); //no. of nodes, note: minimum 3 nodes

    x_start = 0.0;
    x_end = 1.0;

    double *u_approx, *u_exact;
    u_approx = new double[n];
    u_exact = new double[n];

    grid_function_error = executeQuestion(n, x_start, x_end, h, u_approx, u_exact);

    mesh_widths[k] = h;
    errors[k] = grid_function_error;


    if (n==5)
    {
      std::cout << "\nn=5:\n";
      // print u_approx
      std::cout << "\nU approx:\n";
      for (int i=0; i<n; i++)
      {
        std::cout << u_approx[i] <<'\n';
      }

      // print u_exact
      std::cout << "\nU exact:\n";
      for (int i=0; i<n; i++)
      {
        std::cout << u_exact[i] <<'\n';
      }
    }

    //print grid_function_error
    std::cout << "\nh: " << h;
    std::cout << "\nGridFunctionError: " << grid_function_error << "\n";

    delete[] u_approx;
    delete[] u_exact;

  }

  //Create results file
  std::ofstream file;
  file.open("Q1_errors.csv");
  assert(file.is_open());
  file << "h," << "error" << "\n";
  for(int j=0; j<maxIterations; j++)
  {
    file << mesh_widths[j] << "," << errors[j] << "\n";
  }
  file.close();
  // std::string command = "mv Q1_errors.csv Documents/GitHub/";
  // system(command);



  return 0;

}
