#ifndef BTrkLegacy_HelixParams_HH
#define BTrkLegacy_HelixParams_HH
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

// Class interface //
namespace mu2e {
  class HelixParams {
    public:
      enum ParIndex {d0Index=0, phi0Index, omegaIndex, z0Index, tanDipIndex, NHLXPRM};

      HelixParams(const CLHEP::HepVector& pvec, const CLHEP::HepSymMatrix& pcov) : parvec(pvec), parcov(pcov) {}
      ~HelixParams(){}
      HelixParams* clone() const { return new HelixParams(*this); }

      double d0() const                              {return parvec[d0Index];}
      double phi0() const                            {return parvec[phi0Index];}
      double omega() const                           {return parvec[omegaIndex];}
      double z0() const                              {return parvec[z0Index];}
      double tanDip() const                          {return parvec[tanDipIndex];}

      const CLHEP::HepVector& params() const                {return parvec;}
      CLHEP::HepVector& params()                            {return parvec;}
      const CLHEP::HepSymMatrix& covariance() const            {return parcov;}
      CLHEP::HepSymMatrix& covariance()                        {return parcov;}

      void setD0(double in)                          {parvec[d0Index] = in;}
      void setPhi0(double in)                        {parvec[phi0Index] = in;}
      void setOmega(double in)                       {parvec[omegaIndex] = in;}
      void setZ0(double in)                          {parvec[z0Index] = in;}
      void setTanDip(double in)                      {parvec[tanDipIndex] = in;}

    private:
      CLHEP::HepVector parvec;
      CLHEP::HepSymMatrix parcov;
  };
  typedef HelixParams HelixTraj;
}

#endif
