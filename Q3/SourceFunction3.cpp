#include <cmath>
#include <cassert>
#include "SourceFunction3.hpp"


// Constructor
SourceFunction3::SourceFunction3(double nu, double beta, double gamma)
{
  assert(nu > 0);
  assert(beta>=0);
  assert(gamma>=0);

  mNu = nu;
  mBeta = beta;
  mGamma = gamma;
}

// Destructor
SourceFunction3::~SourceFunction3()
{}

// Overriding pure virtual methods
double SourceFunction3::evaluateF(double x, double y)
/* Evaluates f(x,y) */
{
  double u_x, u_xx, u_yy, u;

  u_x = 0.5 * M_PI * cos(M_PI*x) * (1-cos(2*M_PI*y));
  u_xx = -0.5 * pow(M_PI,2.0) * sin(M_PI*x) * (1-cos(2*M_PI*y));
  u_yy = 2.0 * pow(M_PI,2.0) * sin(M_PI*x) * cos(2*M_PI*y);
  u = evaluateUexact(x,y);

  return -1.0*mNu*(u_xx + u_yy) + mBeta*u_x + mGamma*u;
}


double SourceFunction3::evaluateBoundary(double x, double y)
/* Evaluates g(x,y), i.e. boundary function */
{
  return 0.0; //Zero Dirichlet boundaries
}


double SourceFunction3::evaluateUexact(double x, double y)
/* Evaluates u_exact(x,y) */
{
  return 0.5 * sin(M_PI*x) * (1-cos(2*M_PI*y));
}


//Getters
double SourceFunction3::GetNu()
{
  return mNu;
}

double SourceFunction3::GetBeta()
{
  return mBeta;
}

double SourceFunction3::GetGamma()
{
  return mGamma;
}
