#ifndef SOURCEFUNCTION3
#define SOURCEFUNCTION3

#include "Abstract2DFunction.hpp"

class SourceFunction3: public Abstract2DFunction
{
  public:

    // Constructor
    SourceFunction3(double nu, double beta, double gamma);

    // Destructor
    ~SourceFunction3();

    // Overriding pure virtual methods
    double evaluateF(double x, double y);
    double evaluateBoundary(double x, double y);
    double evaluateUexact(double x, double y);

    //Getters
    double GetNu();
    double GetBeta();
    double GetGamma();

  protected:
    double mNu;
    double mBeta;
    double mGamma;

};

#endif
