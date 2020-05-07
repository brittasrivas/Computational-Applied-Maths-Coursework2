#ifndef POISSONPDE
#define POISSONPDE

#include "AbstractPDE.hpp"

class PoissonPDE: public AbstractPDE
{
  public:
    // Specialised Constructor
    PoissonPDE(const int n, const double length, Abstract2DFunction& aFunction);

    // Destructor
    ~PoissonPDE();

  protected:
    void ConstructA();
    void ConstructFvec();


};

#endif
