#ifndef ABSTRACTPDE
#define ABSTRACTPDE

#include <vector>
#include "Abstract2DFunction.hpp"

class AbstractPDE
{

  public:

    void SolvePDE(); // combines auxiliary functions to find Uvec approximation


    // Getters
    double GetErrorNorm();
    double GetH();
    double* GetMesh();
    double* GetUapprox();
    double* GetUexact();



  protected:

    void ConstructMesh();
    virtual void ConstructA() = 0;
    virtual void ConstructFvec() = 0;
    void AugmentA();
    void Gauss();
    void ConstructUvec_exact();
    void CalculateError();
    void ConstructUapprox();
    void ConstructUexact();


    int mN; // number of nodes in mesh along each direction
    int mM; // number of interior nodes in mesh along each direction
    double mLength; // length of mesh in each direction
    double mH; // mesh width

    Abstract2DFunction* mFunction; // function pointer to f(x,y)

    double* mMesh; //mesh nodes in each direction for a square mesh
    std::vector< std::vector<double> > mA; // matrix A
    double* mFvec; // vector F
    double* mUvec; // vector U approximation for interior mesh nodes
    double* mUvec_exact; // vector U exact for interior mesh nodes
    double mErrorNorm; // 2-norm of Uvec approximation error

    double* mUapprox; // vector U approximation for all mesh nodes - used to plot graph
    double* mUexact; // vector U exact for all mesh nodes - used to plot graph




};

#endif
