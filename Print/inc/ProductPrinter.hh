//
//  A base class for classes which print different types of products
//  It serves to simplify the module code a bit
// 
#ifndef Print_inc_ProductPrinter_hh
#define Print_inc_ProductPrinter_hh

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"

#include "CLHEP/Matrix/SymMatrix.h"
#include <boost/io/ios_state.hpp>
#include <ostream>
#include <ios>
#include <iomanip>

namespace mu2e {


  class ProductPrinter {
  public:

    ProductPrinter():_verbose(1) {}
    virtual ~ProductPrinter() {}

    // 0 = none, 1 (default), 2, ..
    void setVerbose(int i) { _verbose = i; }
    int verbose() const {return _verbose; }

    virtual void Print(art::Event const& event,
		       std::ostream& os = std::cout) {}

    void PrintMatrix(const CLHEP::HepSymMatrix& matrix, 
		     std::ostream& os, int mode=0) {
      // when this destructs, it restores the flag state
      boost::io::ios_flags_saver ifs(os);
      // print fixed or scientific
      os.setf(std::ios::floatfield);
      for(int r=1; r<=matrix.num_row(); r++) {
	os << "   ";
	for(int c=1; c<=matrix.num_col(); c++) {
	  if(mode==0) {
	    os << " " << std::setw(13) << std::setprecision(6) << matrix(r,c);
	  } else {
	    double value = sqrt(fabs(matrix(r,r)*matrix(c,c)));
	    if( value!=0.0 ) value = matrix(r,c)/value;
	    os << " " << std::setw(13) << std::setprecision(6) << value;
	  }
	}
	os << "\n";
      }
    }

  private:
    int _verbose;

  };

}
#endif
