#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h> //used for rand()
#include <cmath>

//TODO FUNCTION PROTOTYPES
//TODO MOVE MATRIX FUNCTIONS INTO SEPARATE FILE

void printMatrix(const std::vector< std::vector<double> > M, const int r, const int c)
{
  for (int i=0; i<r; i++)
  {
    for (int j=0; j<c; j++)
    {
      std::cout << M[i][j] << "\t";
    }
    std::cout << "\n";
  }
}

double* matvec_mult(const std::vector< std::vector<double> > M, const double* v)
/* Calculates Matrix-Vector multiplication: M*v
Matrix M is size mxn. */
{
  int m = M.size();
  int n = M[0].size();
  double *ans;
  ans = new double[m];

  for(int i=0; i<m; i++)
  {
    double row_sum = 0.0;
    for (int j=0; j<n; j++)
    {
      row_sum += M[i][j] * v[j];
    }
    ans[i] = row_sum;
  }

  return ans;
}


double* vector_minus(const int n, const double* v1, const double* v2)
/* Calculates v1 - v2 where both are of length n. */
{
  double* ans;
  ans = new double[n];

  for (int i=0; i<n; i++)
  {
    ans[i] = v1[i] - v2[i];
  }

  return ans;
}

double* vector_pointwise_mult(const int n, const double* v1, const double* v2)
/* Calculates v1 .* v2 where both are of length n. */
{
  double* ans;
  ans = new double[n];

  for (int i=0; i<n; i++)
  {
    ans[i] = v1[i] * v2[i];
  }

  return ans;
}

std::vector< std::vector<double> > matrix_transpose(std::vector< std::vector<double> > M)
{
  //TODO figure out function return matrix
}

void initialiseData(double* x1, double* x2, double* y1, double* y2)
{
  x1[0] = 0.1;
  x1[1] = 0.3;
  x1[2] = 0.1;
  x1[3] = 0.6;
  x1[4] = 0.4;
  x1[5] = 0.6;
  x1[6] = 0.5;
  x1[7] = 0.9;
  x1[8] = 0.4;
  x1[9] = 0.7;

  x2[0] = 0.1;
  x2[1] = 0.4;
  x2[2] = 0.5;
  x2[3] = 0.9;
  x2[4] = 0.2;
  x2[5] = 0.3;
  x2[6] = 0.6;
  x2[7] = 0.2;
  x2[8] = 0.4;
  x2[9] = 0.6;

  y1[0] = y1[1] = y1[2] = y1[3] = y1[4] = 1.0;
  y1[5] = y1[6] = y1[7] = y1[8] = y1[9] = 0.0;

  y2[0] = y2[1] = y2[2] = y2[3] = y2[4] = 0.0;
  y2[5] = y2[6] = y2[7] = y2[8] = y2[9] = 1.0;

}

void initialiseParams(std::vector< std::vector<double> > &W2,
  std::vector< std::vector<double> > &W3,
  std::vector< std::vector<double> > &W4,
  double *b2, double *b3, double *b4)
{
  srand((unsigned)time(NULL));

  //Initialise W2
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<2; j++)
    {
      W2[i][j] = (double)(rand()) / (double)(RAND_MAX + 1.0);
    }
  }
  //Initialise W3
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<2; j++)
    {
      W3[i][j] = (double)(rand()) / (double)(RAND_MAX + 1.0);
    }
  }
  //Initialise W4
  for(int i=0; i<2; i++)
  {
    for(int j=0; j<3; j++)
    {
      W4[i][j] = (double)(rand()) / (double)(RAND_MAX + 1.0);
    }
  }


  b2[0] = (double)(rand()) / (double)(RAND_MAX + 1.0);
  b2[1] = (double)(rand()) / (double)(RAND_MAX + 1.0);

  b3[0] = (double)(rand()) / (double)(RAND_MAX + 1.0);
  b3[1] = (double)(rand()) / (double)(RAND_MAX + 1.0);
  b3[2] = (double)(rand()) / (double)(RAND_MAX + 1.0);

  b4[0] = (double)(rand()) / (double)(RAND_MAX + 1.0);
  b4[1] = (double)(rand()) / (double)(RAND_MAX + 1.0);

}

void computeZ(double* z, const int m, const double* a_old, const double* b, const std::vector< std::vector<double> > W)
/* calculates z = Wa + b */
{
  z = matvec_mult(W, a_old);

  for (int i=0; i<m; i++)
  {
    z[i] += b[i];
  }
}



double* activate(double* a_old, std::vector< std::vector<double> > W, double* b)
{
  int m = W.size();
  double *z, *a_new;
  z = new double[m];
  a_new = new double[m];

  computeZ(z, m, a_old, b, W);

  for (int i=0; i<m; i++)
  {
    a_new[i] = 1.0 / (1.0 + exp(-1.0*z[i]));
  }

  delete[] z;
  return a_new;
}


double* activate_deriv(const int n, const double* a)
/* Find derivative of a=sigma(z) where a is a vector of length n. */
{
  double deriv;
  deriv = new double[n];

  for(int i=0; i<n; i++)
  {
    deriv[i] = a[i]*(1.0 - a[i]);
  }

  return deriv;
}



int main(int argc, char* argv[])
{
  //INITIALISE DATA
  int N = 10; // No. of data points

  double *x1, *x2, *y1, *y2;
  x1 = new double[N];
  x2 = new double[N];
  y1 = new double[N];
  y2 = new double[N];

  initialiseData(x1, x2, y1, y2);

  //INITIALISE PARAMETERS
  //double W2[2][2], W3[3][2], W4[2][3]; //TODO REMOVE
  std::vector< std::vector<double> > W2(2, std::vector<double>(2));
  std::vector< std::vector<double> > W3(3, std::vector<double>(2));
  std::vector< std::vector<double> > W4(2, std::vector<double>(3));
  double *b2, *b3, *b4;
  b2 = new double[2];
  b3 = new double[3];
  b4 = new double[2];

  initialiseParams(W2, W3, W4, b2, b3, b4);

  //TODO REMOVE PRINTS
  std::cout << "\nW2: \n";
  printMatrix(W2, 2, 2);
  std::cout << "\nW3: \n";
  printMatrix(W3, 3, 2);
  std::cout << "\nW4: \n";
  printMatrix(W4, 2, 3);


  double eta = 0.05;
  int Niter = 1e6;

  //Initialise arrays for auxiliary steps in SGD
  double *costs_x, *x, *y;
  double *a2, *a3, *a4;
  double *delta2, *delta3, *delta4;
  double *delta2_term2 , *delta3_term2 , *delta4_term2;

  costs = new double[Niter];
  x = new double[2];
  y = new double[2];
  a2 = new double[2];
  a3 = new double[3];
  a4 = new double[2];
  delta2 = new double[2];
  delta3 = new double[3];
  delta4 = new double[2];
  delta2_term2 = new double[2];
  delta3_term2 = new double[3];
  delta4_term2 = new double[2];
  std::vector< std::vector<double> > W2_transp(2, std::vector<double>(2));
  std::vector< std::vector<double> > W3_transp(2, std::vector<double>(3));
  std::vector< std::vector<double> > W4_transp(3, std::vector<double>(2));



  //STOCHASTIC GRADIENT DESCENT
  int random;
  srand((unsigned)time(NULL));

  for (int k=0; k<Niter; k++)
  {
    // Pick an (x,y) data pair at random
    random = rand() % N;
    x[0] = x1[random];
    x[1] = x2[random];
    y[0] = y1[random];
    y[1] = y2[random];

    // Forward pass
    a2 = activate(x,W2,b2);
    a3 = activate(a2,W3,b3);
    a4 = activate(a3,W4,b4);

    // Backward pass
    delta4 = activate_deriv(2, a4);
    delta4_term2 = vector_minus(2, a4, y);
    delta4 = vector_pointwise_mult(delta4, delta4_term2);

    delta3 = activate_deriv(3, a3);




  }






  delete[] x1;
  delete[] x2;
  delete[] y1;
  delete[] y2;
  delete[] b2;
  delete[] b3;
  delete[] b4;
  delete[] costs;
  delete[] x;
  delete[] y;
  delete[] a2;
  delete[] a3;
  delete[] a4;
  delete[] delta2;
  delete[] delta3;
  delete[] delta4;
  delete[] delta2_term2;
  delete[] delta3_term2;
  delete[] delta4_term2;
  //Note: do not need to delete matrices - std::vector destructor will clear them.

  return 0;
}
