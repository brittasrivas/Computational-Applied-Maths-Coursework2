#include <vector>
#include <cmath>
#include <cassert>

#include "ConvDiffPDE.hpp"

//Specialised Constructor
ConvDiffPDE::ConvDiffPDE(const int n, const double length, Abstract2DFunction& aFunction)
{
  assert(n>=3);
  assert(length>0);

  mN = n;
  mM = mN-2;
  mLength = length;
  mH = mLength/double(mM+1);

  mFunction = &aFunction;

  mMesh = new double[mN];

  mA.resize(pow(mM,2), std::vector<double>(pow(mM,2)));

  mFvec = new double[int(pow(mM,2))];
  mUvec = new double[int(pow(mM,2))];
  mUvec_exact = new double[int(pow(mM,2))];

  mUapprox = new double[int(pow(mN,2))];
  mUexact = new double[int(pow(mN,2))];

}


//Destructor
ConvDiffPDE::~ConvDiffPDE()
{}


void ConvDiffPDE::ConstructA()
/* Creates matrix A for finite difference approximation of the Convection-Diffusion-Reaction equation. */
{
  double nu, beta, gamma;

  nu = (*mFunction).GetNu();
  beta = (*mFunction).GetBeta();
  gamma = (*mFunction).GetGamma();

  for (int i=0; i<pow(mM,2); i++)
  {
    mA[i][i] = ((4.0*nu)/pow(mH,2.0)) + gamma; // leading diagonal

    // Other non-zero entries in each row
    if ( ((i+1) < pow(mM,2)) && (((i+1) % mM) != 0) )
    {
      mA[i][i+1] = (-1.0*nu)/pow(mH,2.0) + (beta/mH); //TODO NOTE MOVED BETA TERM HERE RATHER THAN A[i][i]
    }

    if ((i+mM) < pow(mM,2))
    {
      mA[i][i+mM] = (-1.0*nu)/pow(mH,2.0);
    }

    if ( ((i-1)>=0) && ((i % mM) != 0) )
    {
      mA[i][i-1] = ((-1.0*nu)/pow(mH,2.0)) - (beta/mH);
    }

    if ((i-mM) >= 0)
    {
      mA[i][i-mM] = (-1.0*nu)/pow(mH,2.0);
    }

  }
}


void ConvDiffPDE::ConstructFvec()
/* Creates vector F for finite difference approximation of the Convection-Diffusion-Reaction equation. */
{
  int k=0;

  for (int j=1; j<=mM; j++)
  {
    for (int i=1; i<=mM; i++)
    {
      mFvec[k] = (*mFunction).evaluateF(mMesh[i], mMesh[j]); // f(x,y)

      // Since zero Dirichlet boundaries, no need to adjust for boundary terms

      k++;
    }
  }
}
