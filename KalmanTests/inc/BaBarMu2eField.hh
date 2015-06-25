// $Id: BaBarMu2eField.hh,v 1.1 2012/09/19 20:17:37 brownd Exp $
// Description:	Class Header for |BaBarMu2eField|
//              Provide an arbitary fixed field
// Author List:A. Snyder, Copyright (C) 1998	SLAC
#ifndef BaBarMu2eField_HH
#define BaBarMu2eField_HH

#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BField/BField.hh"
// CLHEP
#include "CLHEP/Vector/ThreeVector.h"

// class interface //
namespace mu2e 
{
  class BaBarMu2eField : public BField {

    public:
      //construct from file; optionally amplify the distortions by the given factor
      BaBarMu2eField(CLHEP::Hep3Vector const& origin=CLHEP::Hep3Vector(0.0,0.0,0.0));
      //destroy
      virtual ~BaBarMu2eField();
      // field vector at a point.
      virtual CLHEP::Hep3Vector bFieldVect
	(const HepPoint & point=HepPoint(0,0,0))const;
      // override the nominal field
      virtual double bFieldNominal()const;
     CLHEP::Hep3Vector const& origin() const { return _origin; }
    private:
      mutable double _bnom;
      CLHEP::Hep3Vector _origin;
  };
}
#endif





