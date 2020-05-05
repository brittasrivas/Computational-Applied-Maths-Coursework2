#ifndef ABSTRACTFUNCTION
#define ABSTRACTFUNCTION

/* Class to specify a function f(x) and evaluate the function under set boundary conditions */

class AbstractFunction
{
  public:

    // Class methods to be overridden
    virtual double evaluateF() = 0;

    virtual double evaluateBoundary() = 0;

}

#endif
