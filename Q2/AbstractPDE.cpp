#include <cmath>
#include <vector>
#include "AbstractPDE.hpp"



void AbstractPDE::ConstructMesh()
/* Creates array of mesh nodes in one direction.
Will be used for both x-direction and y-direction since uniform mesh. */
{
  for (int i=0; i<mN; i++)
  {
    mMesh[i] = i * mH;
  }
}

void AbstractPDE::AugmentA()
/* Makes augmented matrix A|b for use by Gauss. */
{
  mA.resize(pow(mM,2), std::vector<double>(pow(mM,2)+1)); // Add another column to A

  for(int i=0; i<pow(mM,2); i++)
  {
    mA[i][pow(mM,2)] = mFvec[i];
  }
}


void AbstractPDE::Gauss()
/* Gauss function taken from martin-thorma.com */
{
  AugmentA(); // A|b

  int n = mA.size(); // number of rows in matrix A

  for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(mA[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(mA[k][i]) > maxEl) {
                maxEl = abs(mA[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1; k++) {
            double tmp = mA[maxRow][k];
            mA[maxRow][k] = mA[i][k];
            mA[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -mA[k][i]/mA[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    mA[k][j] = 0;
                } else {
                    mA[k][j] += c * mA[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    for (int i=n-1; i>=0; i--) {
        mUvec[i] = mA[i][n]/mA[i][i];
        for (int k=i-1; k>=0; k--) {
            mA[k][n] -= mA[k][i] * mUvec[i];
        }
    }
}


void AbstractPDE::ConstructUvec_exact()
/* Creates array of exact values of u(x,y) at interior mesh nodes. */
{
  int k=0;

  for (int j=1; j<mM; j++)
  {
    for (int i=1; i<mM; i++)
    {
      mUexact[k] = (*mFunction).evaluateUexact(mMesh[i], mMesh[j]);

      k++;
    }
  }
}


void AbstractPDE::CalculateError()
/* Grid function norm. Only calculated on interior nodes by definition. */
{
  double sum=0.0;

  for(int i=0; i<pow(mM,2); i++)
  {
    sum += pow((mUvec_exact[i] - mUvec[i]),2.0);
  }

  mErrorNorm = pow((h * sum), 0.5);
}


void AbstractPDE::ConstructUapprox()
/* Full U approx for whole mesh, including mUvec for the interior nodes. */
{
  int k=0;

  for (int j=0; j<mN; j++)
  {
    for (int i=0; i<mN; i++)
    {
      if (i==0 || i==mM+1 || j==0 || j==mM+1)  //Boundaries
      {
        mUapprox[k] = (*mFunction).evaluateBoundary(mMesh[i], mMesh[j]);
      }
      else // Fill in interior nodes using mUvec
      {
        mUapprox[k] = mUvec[(i-1) + ((j-1)*mM)];
      }

      k++;
    }
  }
}


void AbstractPDE::ConstructUexact()
/* Creates array of exact values of u(x,y) for the whole mesh. */
{
  int k=0;

  for (int j=0; j<mN; j++)
  {
    for (int i=0; i<mN; i++)
    {
      mUexact[k] = (*mFunction).evaluateUexact(mMesh[i], mMesh[j]);

      k++;
    }
  }
}



void AbstractPDE::SolvePDE()
/* Combines auxiliary functions to find Uvec approximation and its error. */
{
  ConstructMesh();
  ConstructA();
  ConstructFvec();
  Gauss(); //includes AugmentA()
  ConstructUvec_exact()
  CalculateError();

  // Creates arrays of U approximated and exact over the whole mesh, used for plotting graphs
  ConstructUapprox();
  ConstructUexact();
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

double* AbstractPDE::GetMesh()
{
  return mMesh;
}

double* AbstractPDE::GetUapprox()
{
  return mUapprox;
}

double* AbstractPDE::GetUexact()
{
  return mUexact;
}
