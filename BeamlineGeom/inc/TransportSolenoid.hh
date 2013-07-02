#ifndef BeamlineGeom_TransportSolenoid_hh
#define BeamlineGeom_TransportSolenoid_hh

//
// Class to represent the transport solenoid
//
#include "BeamlineGeom/inc/StraightSection.hh"
#include "BeamlineGeom/inc/Coil.hh"
#include "BeamlineGeom/inc/Collimator_TS1.hh"
#include "BeamlineGeom/inc/Collimator_TS3.hh"
#include "BeamlineGeom/inc/Collimator_TS5.hh"
#include "BeamlineGeom/inc/PbarWindow.hh"
#include "BeamlineGeom/inc/TorusSection.hh"

// C++ includes
#include <algorithm>
#include <map>
#include <memory>
#include <vector>

// cet
#include "cetlib/exception.h"

namespace mu2e {

  class TransportSolenoid {
    
    template <typename T>
    using vector_unique_ptr = std::vector<std::unique_ptr<T>>;
    
    friend class BeamlineMaker;

  public:
    TransportSolenoid() :
      _rTorus(0.), _rVac(0.)
   {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor

    // - only the following enums should be used for the following
    //   accessors
    enum enum_type_ts  { TS1, TS2, TS3, TS4, TS5 };
    enum enum_type_ext { IN, OUT };

    // Cryo-stat dimensions
    double torusRadius() const { return _rTorus; }
    double innerRadius() const { return _rVac; }

    double endWallU1_rIn() const { return _rIn_endWallU1; }
    double endWallU2_rIn() const { return _rIn_endWallU2; }
    double endWallD_rIn() const  { return _rIn_endWallD;  }

    double endWallU1_rOut() const { return _rOut_endWallU1; }
    double endWallU2_rOut() const { return _rOut_endWallU2; }
    double endWallD_rOut() const  { return _rOut_endWallD;  }

    double endWallU1_halfLength() const { return _halfLength_endWallU1; }
    double endWallU2_halfLength() const { return _halfLength_endWallU2; }
    double endWallD_halfLength() const  { return _halfLength_endWallD;  }

    std::string material()       const { return _material; }
    std::string insideMaterial() const { return _insideMaterial; }

    // make non-const to allow type casting
    TSSection* getTSCryo(enum_type_ts i,enum_type_ext j) const { return _cryoMap.find(i)->second.at(j); }

    StraightSection const& getTS1_in()  const { return _ts1in;  }
    StraightSection const& getTS1_out() const { return _ts1out; }
    StraightSection const& getTS3_in()  const { return _ts3in;  }
    StraightSection const& getTS3_out() const { return _ts3out; }
    StraightSection const& getTS5_in()  const { return _ts5in;  }
    StraightSection const& getTS5_out() const { return _ts5out; }

    TorusSection    const& getTS2_in()  const { return _ts2in;  }
    TorusSection    const& getTS2_out() const { return _ts2out; }
    TorusSection    const& getTS4_in()  const { return _ts4in;  }
    TorusSection    const& getTS4_out() const { return _ts4out; }

    // Coils 
    std::string coil_material() const { return _coilMaterial; }
    unsigned getNCoils(enum_type_ts i) const { return _nCoils.at( (unsigned)i); }
    const vector_unique_ptr<Coil>& getTSCoils(enum_type_ts i) const      
    { return _coilMap.find(i)->second; }

    // Collimators
    CollimatorTS1 const& getColl1()  const { return _coll1;  }
    CollimatorTS3 const& getColl31() const { return _coll31; }
    CollimatorTS3 const& getColl32() const { return _coll32; }
    CollimatorTS5 const& getColl5()  const { return _coll5;  }

    PbarWindow const& getPbarWindow() const { return _pbarWindow; }

  protected:

    // Cryostat
    double _rTorus;
    double _rVac;

    double _rIn_endWallU1; 
    double _rIn_endWallU2; 
    double _rIn_endWallD;  

    double _rOut_endWallU1; 
    double _rOut_endWallU2; 
    double _rOut_endWallD;  

    double _halfLength_endWallU1; 
    double _halfLength_endWallU2; 
    double _halfLength_endWallD;  

    std::string _material;
    std::string _insideMaterial;

    std::map<int,std::vector<TSSection*>> _cryoMap;

    StraightSection _ts1in ;
    StraightSection _ts1out;
    StraightSection _ts3in ;
    StraightSection _ts3out;
    StraightSection _ts5in ;
    StraightSection _ts5out;

    TorusSection    _ts2in ;
    TorusSection    _ts2out;
    TorusSection    _ts4in ;
    TorusSection    _ts4out;


    // Coils
    std::string _coilMaterial;
    std::map<int,vector_unique_ptr<Coil>> _coilMap;
    const std::vector<unsigned> _nCoils { 3, 18, 8, 18, 5 };
    
    // Collimators
    CollimatorTS1 _coll1;
    CollimatorTS3 _coll31;
    CollimatorTS3 _coll32;
    CollimatorTS5 _coll5;

    PbarWindow _pbarWindow;

  public:

    // Helper utility
    static void checkSizeOfVectors( const unsigned nCoils, 
                                    const std::initializer_list<std::vector<double>>& tmp_vecs ) {
      
      int sizeMismatches(0);
      std::for_each( tmp_vecs.begin(),
                     tmp_vecs.end(),
                     [nCoils,&sizeMismatches](const std::vector<double>& tmp_vec){
                       if ( tmp_vec.size() != nCoils ) sizeMismatches++;
                     } );
      
      if ( sizeMismatches > 0 ) {
        throw cet::exception("VECTOR") 
          << "Mismatch in specified TS coil parameters and required number of coils \n" 
          << "Should be " << nCoils; 
        
      }
      
    }
    
  };

}
#endif /* BeamlineGeom_TransportSolenoid_hh */
