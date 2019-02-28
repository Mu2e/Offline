//NOTE THIS WILL BE USED FOR ALIGNMENT PURPOSES LATER ON
#ifndef Mu2eUtilities_PolyCoeff_hh
#define Mu2eUtilities_PolyCoeff_hh

//c++
#include <vector>
//ROOT
#include "Math/VectorUtil.h"
#include "TMatrixD.h"

using namespace std;

namespace mu2e {
  
  class PolyCoeff{
  
	public:
     
      PolyCoeff(std::vector<int> inVariablesByVector,
                                            int outVariable, double coefficient){}
      
      std::vector<int> InVariables() const {return _inVarByVec;}

      int OutVariable() const {return _outVar;}

      
      double Coefficient() const {return _coefficient;}

      
      std::vector<int> InVariables(std::vector<int> inVar ) {
        _inVarByVec  = inVar;
        return _inVarByVec;
      }

      int OutVariable(int  outVar) {
        _outVar      = outVar;
        return _outVar;
      }

      
      double Coefficient(double coeff) {
        _coefficient = coeff;
        return _coefficient;
      }

     
    private:
      std::vector<int> _inVarByVec;
      int              _outVar;
      double           _coefficient;
  };

       
}//end mu2e
#endif

