#include "Offline/GeneralUtilities/inc/TwoDWeight.hh"
#include "TMath.h"
#include <stdexcept>
namespace mu2e {
  using SVEC = ROOT::Math::SVector<double,2>;
  using SMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;

  TwoDWeight::TwoDWeight(TwoDPoint const& pos,float ivar) {
    wt_ = pos.cov();
    if(ivar>0.0){
      wt_(0,0) += ivar;
      wt_(1,1) += ivar;
    }
    if(!wt_.Invert())throw std::invalid_argument( "Inversion Failure" );
    wtpos_ = wt_*pos.pos();
  }

  TwoDPoint TwoDWeight::point(float ivar) const {
    int ifail(0);
    SMAT cov = wt_.Inverse(ifail);
    if(ifail != 0)throw std::invalid_argument( "Inversion Failure" );
    auto pos = cov*wtpos_;
    if(ivar>0.0){
      cov(0,0) -= ivar;
      cov(1,1) -= ivar;
    }
    return TwoDPoint(pos,cov);
  }

  TwoDWeight& TwoDWeight::operator +=(TwoDWeight const& other) {
    wt_ += other.wt();
    wtpos_ += other.wtPos();
    return *this;
  }

  TwoDWeight& TwoDWeight::operator -=(TwoDWeight const& other) {
    wt_ -= other.wt();
    wtpos_ -= other.wtPos();
    return *this;
  }

  void TwoDWeight::print(std::ostream& os) const {
    os << *this << std::endl;
  }
}

std::ostream& operator << (std::ostream& ost, mu2e::TwoDWeight const& wt) {
  ost << "TwoDWeight position (" << wt.weightPosition()[0] << "," << wt.weightPosition()[1] << "), upper weight matrix ("
    << wt.weight()(0,0) << "," <<  wt.weight()(0,1) << "," << wt.weight()(1,1) << ") ";
  return ost;
}
