#include <vector>
#include <cassert>

#include "PoissonPDE.hpp"

//Specialised Constructor
PoissonPDE::PoissonPDE(const int n, const double length, Abstract2DFunction& aFunction)
{
  assert(n>=2);
  assert(length>0);
  mN = n;
  mM = mN-2;
  mLength = length;
  mH = mLength/double(mM+1);

  mFunction = &aFunction;

  mMesh = new double[mN];

  mA.resize(pow(mM,2), std::vector<double>(pow(mM,2)));

  mFvec = new double[pow(mM,2)];
  mUvec = new double[pow(mM,2)];
  mUvec_exact = new double[pow(mM,2)];

  mUapprox = new double[pow(mN,2)];
  mUexact = new double[pow(mN,2)];

}


//Destructor
PoissonPDE::~PoissonPDE()
{}


PoissonPDE::void ConstructA()
/* Creates matrix A for finite difference approximation of the Poisson equation. */
{
  for (int i=0; i<pow(mM,2); i++)
  {
    mA[i][i] = 4.0/pow(mH,2.0); // leading diagonal

    // Other non-zero entries
    if ( ((i+1) < pow(mM,2)) && (((i+1) % mM) != 0) )
    {
      mA[i][i+1] = -1.0/pow(mH,2.0);
    }

    if ((i+mM) < pow(mM,2))
    {
      mA[i][i+mM] = -1.0/pow(mH,2.0);
    }

    if ( ((i-1)>=0) && ((i % mM) != 0) )
    {
      mA[i][i-1] = -1.0/pow(mH,2.0);
    }

    if ((i-mM) >= 0)
    {
      mA[i][i-mM] = -1.0/pow(mH,2.0);
    }

  }
}


PoissonPDE::void ConstructFvec()
/* Creates vector F for finite difference approximation of the Poisson equation. */
{
  int k=0;

  for (int j=1; j<=mM; j++)
  {
    for (int i=1; i<=mM; i++)
    {
      Fvec[k] = (*mFunction).evaluateF(mMesh[i], mMesh[j]); // f(x,y)

      // Adjust for boundary terms
      
      if (i==1) // add g(0,hj)/(h^2)
      {
        Fvec[k] += (*mFunction).evaluateBoundary(mMesh[i-1], (mH*j)) / pow(mH,2.0);
      }

      if (i==mM) // add g((m+1)h,hj)/(h^2)
      {
        Fvec[k] += (*mFunction).evaluateBoundary(mMesh[i+1], (mH*j)) / pow(mH,2.0);
      }

      if (j==1) // add g(hi,0)/(h^2)
      {
        Fvec[k] += (*mFunction).evaluateBoundary((mH*i), mMesh[j-1]) / pow(mH,2.0);
      }

      if (j==mM) // add g(hi,(m+1)h)/(h^2)
      {
        Fvec += (*mFunction).evaluateBoundary((mH*i), mMesh[j+1]) / pow(mH,2.0);
      }

      k++;
    }
  }
}
