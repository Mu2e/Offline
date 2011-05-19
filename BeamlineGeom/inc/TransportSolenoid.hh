#ifndef BeamlineGeom_TransportSolenoid_hh
#define BeamlineGeom_TransportSolenoid_hh

//
// Class to represent the transport solenoid
//
#include "BeamlineGeom/inc/StraightSection.hh"
#include "BeamlineGeom/inc/Collimator.hh"

namespace mu2e {

  class TransportSolenoid {

  friend class BeamlineMaker;

  public:
    TransportSolenoid() :
      _rTorus(0.0), _rVac(0.0), _rCryo(0.0),
      _ts1(), _ts3(), _ts5(),
      _coll1(), _coll31(), _coll32(), _coll5()
   {;}

    // use compiler-generated copy c'tor, copy assignment, and d'tor

    double torusRadius() const {return _rTorus;}
    double innerRadius() const {return _rVac;}
    double outerRadius() const {return _rCryo;}

    StraightSection const& getTS1() const {return _ts1;}
    StraightSection const& getTS3() const {return _ts3;}
    StraightSection const& getTS5() const {return _ts5;}

    Collimator const& getColl1()  const { return _coll1;  }
    Collimator const& getColl31() const { return _coll31; }
    Collimator const& getColl32() const { return _coll32; }
    Collimator const& getColl5()  const { return _coll5;  }

  protected:

    // All dimensions in mm.

    double _rTorus;
    double _rVac;
    double _rCryo;

    StraightSection _ts1;
    StraightSection _ts3;
    StraightSection _ts5;

    Collimator _coll1;
    Collimator _coll31;
    Collimator _coll32;
    Collimator _coll5;
};

}
#endif /* BeamlineGeom_TransportSolenoid_hh */
