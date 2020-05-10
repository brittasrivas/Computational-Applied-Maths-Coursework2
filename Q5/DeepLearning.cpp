#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h> //used for rand()

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
  double Niter = 1e6;

  //Initialise arrays for auxiliary steps in SGD
  double *costs_x, *x, *y;
  double *a2, *a3, *a4;
  double *delta1, *delta2, *delta3;

  costs_x = new double[Niter];
  x = new double[2];
  y = new double[2];
  //TODO ALLOCATE MEMORY FOR A AND DELTAS


  //STOCHASTIC GRADIENT Descent
  int random;
  srand((unsigned)time(NULL));

  for (int k=0; k<Niter; k++)
  {
    random = rand() % N;
    x[0] = x1[random];
    x[1] = x2[random];
    y[0] = y1[random];
    y[1] = y2[random];


  }






  delete[] x1;
  delete[] x2;
  delete[] y1;
  delete[] y2;
  delete[] b2;
  delete[] b3;
  delete[] b4;
  delete[] costs_x;
  delete[] x;
  delete[] y;
  //Note: do not need to delete matrices - std::vector destructor will clear them.

  return 0;
}
