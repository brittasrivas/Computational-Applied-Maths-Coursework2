#ifndef CONVDIFFPDE
#define CONVDIFFPDE

#include "AbstractPDE.hpp"

class ConvDiffPDE: public AbstractPDE
{
  public:
    // Specialised Constructor
    ConvDiffPDE(const int n, const double length, Abstract2DFunction& aFunction);

    // Destructor
    ~ConvDiffPDE();

  protected:
    // Overiding pure virtual methods
    void ConstructA();
    void ConstructFvec();

};

#endif
