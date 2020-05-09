#include <iostream>
#include <cmath>
#include <fstream>

//FUNCTION PROTOTYPES

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


int main(int argc, char* argv[])
{
  const int n = 3;
  double eta = 0.75;
  int MaxIter = 64;
  double w = 0.5;
  double b = 0.5;


  double *xs, *ys;
  xs = new double[n];
  ys = new double[n];
  xs[0] = 0.0;
  xs[1] = 0.5;
  xs[2] = 1.0;
  ys[0] = 0.0;
  ys[1] = 1.0;
  ys[2] = 1.0;

  double cost;
  cost = cost_function(n, ys, xs, w, b);


  //Create results file
  std::ofstream file;
  file.open("Q4_StochasticGradientDescent.csv");
  assert(file.is_open());
  file << "iteration," << "w," << "b," << "cost" << "\n";
  file << "0," << w << "," << b << "," << cost << "\n";   //Initial values


  // Stochastic Gradient Descent Algorithm
  double w_temp, b_temp;
  int random;

  for(int k=0; k<MaxIter; k++)
  {
    random = rand() % n;

    w_temp = w - eta*deriv_Cx_w(ys[random], xs[random], w, b);
    b_temp = b - eta*deriv_Cx_b(ys[random], xs[random], w, b);

    w = w_temp;
    b = b_temp;
    cost = cost_function(n, ys, xs, w, b);

    file << k+1 << "," << w << "," << b << "," << cost << "\n";
  }


  //Close file
  file.close();
  std::string command = "mv Q4_StochasticGradientDescent.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q4";
  system(command.c_str());


  delete[] xs;
  delete[] ys;

  return 0;
}
