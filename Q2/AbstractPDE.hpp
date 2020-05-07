#ifndef AbstractPDE
#define AbstractPDE

#include <vector>
#include "Abstract2DFunction.hpp"

class AbstractPDE
{

  public:

    void SolvePDE();


    // Getters
    double GetErrorNorm();
    double GetH();
    double* GetUapprox();
    double* GetUexact();



  protected:

    void ConstructMesh();
    virtual void ConstructA() = 0;
    virtual void ConstructFvec() = 0;
    void Gauss();

    int mN; // number of nodes in mesh along each direction
    int mM; // number of interior nodes in mesh along each direction
    double mLength; // length of mesh in each direction
    double mH; // mesh width

    Abstract2DFunction* mFunction; // function pointer to f(x,y)

    double* mMesh; //mesh in each direction for a square mesh
    std::vector< std::vector<double> > mA; // matrix A
    double* mFvec; // vector F
    double* mUvec; // vector U approximation for interior mesh nodes
    double* mUapprox; // vector U approximation for all mesh nodes
    double* mUexact; // vector U exact

    double mErrorNorm; // 2-norm of Uapprox error


};

#endif
