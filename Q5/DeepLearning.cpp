#include <iostream>
#include <cmath>
#include <vector>
#include <stdlib.h> //used for rand()

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

void initialiseParams(double W2[2][2], double W3[3][2], double W4[2][3],
  double *b2, double *b3, double *b4)
{
  W2[0][0] = rand() % 5000;
  W2[0][1] = ;
  W2[1][0] = ;
  W2[1][1] = ;


}


int main(int argc, char* argv[])
{
  int N = 10; // No. of data points

  double *x1, *x2, *y1, *y2;
  x1 = new double[N];
  x2 = new double[N];
  y1 = new double[N];
  y2 = new double[N];

  initialiseData(x1, x2, y1, y2);

  double W2[2][2], W3[3][2], W4[2][3];
  double *b2, *b3, *b4;
  b2 = new double[2];
  b3 = new double[3];
  b4 = new double[2];




  delete[] x1;
  delete[] x2;
  delete[] y1;
  delete[] y2;

  return 0;
}
