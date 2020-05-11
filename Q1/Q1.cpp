#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

//Function Prototypes
double* createMesh(const int n, const double h, const double x_start,
  const double x_end);
double* solveTridiagonal(const int m, double *d, double *u, double *l, double *Fvec);
void calculateGridFunctionError(const int m, const double h,
  const double* Uvec_full, const double* Uvec_exact, double &grid_function_error);
void testCoeffcientMatrix(const int m, const double h, double* A_d, double* A_u, double* A_l);
double testExactU(double x);
void outputExactUvec(const int n, const double* mesh, double* Uvec_exact);
double testExactF(double x);
void testFvec(const int m, const double alpha, const double beta, const double h,
  const double* mesh, double* Fvec);
double solvePDE(const int n, const double x_start, const double x_end,
  double &h, double* Uvec_full, double* Uvec_exact, double* mesh);



/*----------------------------GENERAL FUNCTIONS-------------------------------*/

double* createMesh(const int n, const double h, const double x_start,
  const double x_end)
/* Outputs mesh of n nodes between x_start and x_end with mesh width h. */
{
  double* mesh;
  mesh = new double[n];


  for (int j=0; j<n; j++)
  {
    mesh[j] = j*h;
  }

  return mesh;
}


double* solveTridiagonal(const int m, double *d, double *u, double *l,
  double *Fvec)
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


void calculateGridFunctionError(const int m, const double h,
  const double* Uvec_full, const double* Uvec_exact,
  double &grid_function_error)
/* Calculates grid function error norm between approx and exact functions. */
{
  double sum=0.0;

  for(int i=1; i<=m; i++)
  {
    sum += pow((Uvec_exact[i] - Uvec_full[i]),2.0);
  }

  grid_function_error = pow((h * sum), 0.5);
}


/*----------------CODE VERIFICATION TEST FUNCTIONS----------------------------*/

void testCoeffcientMatrix(const int m, const double h, double* A_d, double* A_u,
  double* A_l)
/* Constructs Matrix A */
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
/* Exact u(x) function */
{
  return x - sin(M_PI*x);
}


void outputExactUvec(const int n, const double* mesh, double* Uvec_exact)
/* Evaluates exact u function for the interior mesh nodes and
saves in Uvec_exact */
{
  for(int i=0; i<n; i++)
  {
    Uvec_exact[i] = testExactU(mesh[i]);
  }
}


double testExactF(double x)
/* Evaluates f(x) */
{
  return (-1 * pow(M_PI,2.0) * sin(M_PI*x));
}


void testFvec(const int m, const double alpha, const double beta,
  const double h, const double* mesh, double* Fvec)
/* Constructs vector F */
{
   for (int i=0; i<m; i++)
   {
     Fvec[i] = testExactF(mesh[i+1]);
     //Note use mesh[i+1] since mesh is of length n, whilst Fvec is of length m
     //(i.e. excludes first and last nodes in mesh)
   }

   Fvec[0] += alpha/(pow(h,2.0));
   Fvec[m-1] += beta/(pow(h,2.0));
}



/*---------------------------PDE SOLVE FUNCTION-------------------------------*/

double solvePDE(const int n, const double x_start, const double x_end,
  double &h, double* Uvec_full, double* Uvec_exact, double* mesh)
/* Solves AU = F and returns the error norm of the approximation. */
{
  int m = n-2;
  h = (x_end - x_start)/double(m+1);
  double alpha = 0.0;
  double beta = 1.0;

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

  // U approx, including boundaries - used to plot the graph
  Uvec_full[0] = alpha;
  Uvec_full[n-1] = beta;

  for(int j=1; j<=m; j++)
  {
    Uvec_full[j] = Uvec[j-1];
  }


  // U exact, including boundaries - used to plot the graph
  outputExactUvec(n, mesh, Uvec_exact);


  // Error
  double grid_function_error;
  calculateGridFunctionError(m, h, Uvec_full, Uvec_exact, grid_function_error);

  delete[] A_d;
  delete[] A_u;
  delete[] A_l;
  delete[] Fvec;
  delete[] Uvec;

  return grid_function_error;
}


/*----------------------------------MAIN--------------------------------------*/

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

    double *mesh;
    mesh = new double[n]; //mesh of length n

    double *u_approx, *u_exact;
    u_approx = new double[n];
    u_exact = new double[n];

    grid_function_error = solvePDE(n, x_start, x_end, h, u_approx,
      u_exact, mesh);

    mesh_widths[k] = h;
    errors[k] = grid_function_error;


    if (n==5)
    {
      mesh = createMesh(n, h, x_start, x_end);

      //Create results file for plots when n=5
      std::ofstream file2;
      file2.open("Q1_function.csv");
      assert(file2.is_open());
      file2 << "mesh," << "u_approx," << "u_exact" << "\n";

      for(int i=0; i<n; i++)
      {
        file2 << mesh[i] << "," << u_approx[i] << "," << u_exact[i] << "\n";
      }
      file2.close();
      std::string command2 = "mv Q1_function.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q1";
      system(command2.c_str());

    }


    delete[] u_approx;
    delete[] u_exact;
    delete[] mesh;

  }


  //Create results file
  std::ofstream file;
  file.open("Q1_errors.csv");
  assert(file.is_open());
  file << "h," << "error," << "log(h)," << "log(error)," << "log(h^2)" << "\n";
  for(int j=0; j<maxIterations; j++)
  {
    file << mesh_widths[j] << "," << errors[j] << ","
      << log(mesh_widths[j]) << "," << log(errors[j]) << ","
      << log(pow(mesh_widths[j],2.0)) << "\n";
  }
  file.close();
  std::string command = "mv Q1_errors.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q1";
  system(command.c_str());


  return 0;

}
