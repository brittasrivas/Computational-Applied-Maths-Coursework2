#include <vector>

#include "EllipticalPDE.hpp"

//Specialised Constructor
EllipticalPDE::EllipticalPDE(const int n, const double length, Abstract2DFunction& aFunction)
{
  mN = n;
  mM = mN-2;
  mLength = length;
  mH = mLength/double(mM+1);

  mFunction = &aFunction;

  mMesh.resize(mM, std::vector<double>(mM));

  mA.resize(pow(mM,2), std::vector<double>(pow(mM,2)));

  mFvec = new double[pow(mM,2)];
  mUvec = new double[pow(mM,2)];

  mUapprox = new double[pow(mN,2)];
  mUexact = new double[pow(mN,2)];


}
