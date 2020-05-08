#ifndef ABSTRACT2DFUNCTION
#define ABSTRACT2DFUNCTION

/* Class to specify a function f(x,y) and evaluate the function under set boundary conditions */

class Abstract2DFunction
{
  public:

    // Class methods to be overridden
    virtual double evaluateF(double x, double y) = 0;

    virtual double evaluateBoundary(double x, double y) = 0;

    virtual double evaluateUexact(double x, double y) = 0;

};

#endif
