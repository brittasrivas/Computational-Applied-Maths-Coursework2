#include "AbstractPDE.hpp"


void AbstractPDE::ConstructMesh()
{
  for (int i=0; i<mN; i++)
  {
    mMesh[i] = i * mH;
  }
}


void AbstractPDE::Gauss()
{
  //TODO
}



void AbstractPDE::SolvePDE()
{
  //TODO
}






// Getters

double AbstractPDE::GetErrorNorm()
{
  return mErrorNorm;
}

double AbstractPDE::GetH()
{
  return mH;
}

double* AbstractPDE::GetUapprox()
{
  return mUapprox;
}

double* AbstractPDE::GetUexact()
{
  return mUexact;
}
