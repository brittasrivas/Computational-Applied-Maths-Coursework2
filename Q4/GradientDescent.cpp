#include <iostream>
#include <cmath>
#include <fstream>

//FUNCTION PROTOTYPES

double cost_function(const int n, const double* y, const double* x, const double w, const double b)
{
  double sum = 0.0;
  double diff;

  for (int i=0; i<n; i++)
  {
    diff = y[i] - (w*x[i] + b);
    sum += pow(diff, 2.0);
  }

  return sum/double(n*2.0);
}


double deriv_w(const int n, const double* y, const double* x, const double w, const double b)
{
  double sum = 0.0;

  for (int i=0; i<n; i++)
  {
    sum += (w*x[i] + b - y[i]) * x[i];
  }

  return sum/double(n);
}


double deriv_b(const int n, const double* y, const double* x, const double w, const double b)
{
  double sum = 0.0;

  for (int i=0; i<n; i++)
  {
    sum += w*x[i] + b - y[i];
  }

  return sum/double(n);
}


int main(int argc, char* argv[])
{
  const int n = 3;
  double eta = 0.75;
  int MaxIter = 64;
  double w = 0.5;
  double b = 0.5;


  double *x, *y;
  x = new double[n];
  y = new double[n];
  x[0] = 0.0;
  x[1] = 0.5;
  x[2] = 1.0;
  y[0] = 0.0;
  y[1] = 1.0;
  y[2] = 1.0;

  double cost;
  cost = cost_function(n, y, x, w, b);


  //Create results file
  std::ofstream file;
  file.open("Q4_GradientDescent.csv");
  assert(file.is_open());
  file << "iteration," << "w," << "b," << "cost" << "\n";
  file << "0," << w << "," << b << "," << cost << "\n";   //Initial values


  // Gradient Descent Algorithm
  double w_temp, b_temp;

  for(int k=0; k<MaxIter; k++)
  {
    w_temp = w - eta*deriv_w(n, y, x, w, b);
    b_temp = b - eta*deriv_b(n, y, x, w, b);

    w = w_temp;
    b = b_temp;
    cost = cost_function(n, y, x, w, b);

    file << k+1 << "," << w << "," << b << "," << cost << "\n";
  }


  //Close file
  file.close();
  std::string command = "mv Q4_GradientDescent.csv Documents/GitHub/Computational-Applied-Maths-Coursework2/Q4";
  system(command.c_str());


  delete[] x;
  delete[] y;

  return 0;
}
