#ifndef HELIXTRAJ_HH
#define HELIXTRAJ_HH
#include "CLHEP/Matrix/Vector.h"
#include <CLHEP/Matrix/SymMatrix.h>


class HelixParams;

class HelixTraj{
private:
    CLHEP::HepVector pvec_;
    CLHEP::HepSymMatrix pcov_;
public:
  HelixTraj(const CLHEP::HepVector& pvec,const CLHEP::HepSymMatrix& pcov ) : pvec_(pvec), pcov_(pcov)
  {
    //  Make sure the dimensions of the input matrix and vector are correct
    /*
    if( pvec.num_row() != NHLXPRM || pcov.num_row() != NHLXPRM ){
    ErrMsg(fatal) 
        << "HelixTraj: incorrect constructor vector/matrix dimension" << endmsg;
    }
    */

    if (omega() == 0.0) pvec_[omegaIndex] = 1.e-9;
  }
  HelixTraj( const HelixTraj& h ) : pvec_(h.pvec_), pcov_(h.pcov_) 
  {
  }


  HelixTraj* clone() const
  {
    return new HelixTraj(*this);
  }

  enum ParIndex {d0Index=0, phi0Index, omegaIndex, z0Index, tanDipIndex, NHLXPRM};
  double d0() const {return pvec_[d0Index];}
  double phi0() const;
  double omega() const {return pvec_[omegaIndex]; }
  double z0() const {return pvec_[z0Index]; }
  double tanDip() const {  return pvec_[tanDipIndex]; }


  double z( const double& ) const;
  double zFlight(double zpos) const;
  double dip() const {return atan(tanDip());}
  double cosDip() const {return 1./sqrt(1.+(tanDip()*tanDip())); }
  double sinDip() const {return tanDip()*cosDip(); }
  double translen(const double& f) const {return cosDip()*f;}
  double arc( const double& f) const {return translen(f)*omega();}
  double angle(const double& f) const;
};
#endif

