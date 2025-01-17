#include <iostream>
#include <cmath>
#include <fstream>
#include <stdlib.h> //used for rand()

//FUNCTION PROTOTYPES
double Func(const double x, const double w, const double b);
double cost_x(const double y, const double x, const double w, const double b);
double cost_function(const int n, const double* ys, const double* xs, const double w, const double b);
double deriv_Cx_w(const double y, const double x, const double w, const double b);
double deriv_Cx_b(const double y, const double x, const double w, const double b);
void stochasticGradientDescent(const int n, const double eta, const int MaxIter,
  double* ys, double* xs, double* costs, double* ws, double* bs);



double Func(const double x, const double w, const double b)
{
  return w*x + b;
}

double cost_x(const double y, const double x, const double w, const double b)
{
  double diff = y - Func(x, w, b);
  return 0.5 * pow(diff, 2.0);
}

double cost_function(const int n, const double* ys, const double* xs, const double w, const double b)
{
  double sum = 0.0;

  for (int i=0; i<n; i++)
  {
    sum += cost_x(ys[i], xs[i], w, b);
  }

  return sum/double(n);
}

double deriv_Cx_w(const double y, const double x, const double w, const double b)
{
  return (w*x + b - y) * x;
}


double deriv_Cx_b(const double y, const double x, const double w, const double b)
{
  return w*x + b - y;
}


void stochasticGradientDescent(const int n, const double eta, const int MaxIter,
  double* ys, double* xs, double* costs, double* ws, double* bs)
{
  int random;

  for(int k=0; k<MaxIter; k++)
  {
    random = rand() % n;

    ws[k+1] = ws[k] - eta*deriv_Cx_w(ys[random], xs[random], ws[k], bs[k]);
    bs[k+1] = bs[k] - eta*deriv_Cx_b(ys[random], xs[random], ws[k], bs[k]);

    costs[k+1] = cost_function(n, ys, xs, ws[k+1], bs[k+1]);
  }
}


int main(int argc, char* argv[])
{
  const int n = 3;
  const double w_0 = 0.5;
  const double b_0 = 0.5;


  double *xs, *ys;
  xs = new double[n];
  ys = new double[n];
  xs[0] = 0.0;
  xs[1] = 0.5;
  xs[2] = 1.0;
  ys[0] = 0.0;
  ys[1] = 1.0;
  ys[2] = 1.0;

  double *costs, *ws, *bs;



  double eta = 0.75; //learning rate
  int MaxIter = 64;


  for (int i=0; i<4 ; i++)
  {
    // Max iterations experiment
    if (i==0)
    {
      MaxIter = 64;
    }
    else if (i>0)
    {
      MaxIter = MaxIter*3;
    }

    costs = new double[MaxIter+1];
    ws = new double[MaxIter+1];
    bs = new double[MaxIter+1];
    ws[0] = w_0;
    bs[0] = b_0;
    costs[0] = cost_function(n, ys, xs, ws[0], bs[0]);

    // Learning rate experiment
    for (int j=0; j<4; j++)
    {
      if (j==0)
      {
        eta = 0.75;
      }
      else if (j>0)
      {
        eta = eta/3.0;
      }

      // Stochastic Gradient Descent Algorithm
      stochasticGradientDescent(n, eta, MaxIter, ys, xs, costs, ws, bs);

      std::cout << "\n\nMaxIter: " << MaxIter;
      std::cout << "\neta: " << eta;
      std::cout << "\n                w: " << ws[MaxIter];
      std::cout << "\n                b: " << bs[MaxIter];

    }

  }

  // From the experiments, taking MaxIter: 1728 & eta: 0.0277778
  // converged close to w=1 & b=1/6 as per the LSQ minimizer estimation


  //Create results file
  std::ofstream file;
  file.open("Q4_StochasticGradientDescent.csv");
  assert(file.is_open());
  file << "iteration," << "w," << "b," << "cost" << "\n";

  for (int k=0; k<(MaxIter+1); k++)
  {
    file << k << "," << ws[k] << "," << bs[k] << "," << costs[k] << "\n";
  }

  file.close();
  std::string command = "mv Q4_StochasticGradientDescent.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q4";
  system(command.c_str());


  delete[] xs;
  delete[] ys;
  delete[] costs;
  delete[] ws;
  delete[] bs;

  return 0;
}
