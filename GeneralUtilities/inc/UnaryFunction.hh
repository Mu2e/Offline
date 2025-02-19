// Ed Callaghan
// Interface for a function of a single variable
// February 2025

#ifndef GeneralUtilities_UnaryFunction_hh
#define GeneralUtilities_UnaryFunction_hh

namespace mu2e{
  class UnaryFunction{
    public:
      UnaryFunction() = default;
     ~UnaryFunction() = default;

      // concrete implementations should override
      virtual double Evaluate(double) = 0;
      double operator() (double x){
        double rv = this->Evaluate(x);
        return rv;
      };
    protected:
      /**/
    private:
      /**/
  };
};

#endif
