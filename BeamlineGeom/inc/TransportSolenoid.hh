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
#include "BeamlineGeom/inc/ConeSection.hh"
#include "BeamlineGeom/inc/TSSection.hh"

#include "GeneralUtilities/inc/EnumToStringSparse.hh"

// C++ includes
#include <map>
#include <memory>
#include <vector>


namespace mu2e {

  class TransportSolenoid {
    
    friend class BeamlineMaker;

  public:
    TransportSolenoid() :
      _rTorus(0.), _rVac(0.)
    {
      // Reserve number of coils
      for ( unsigned iTS = TSRegion::TS1 ; iTS <= TSRegion::TS5 ; ++iTS ) 
        _coilMap[ (TSRegion)iTS ].reserve( getNCoils( (TSRegion)iTS ) );
    }
    
    // use compiler-generated copy c'tor, copy assignment, and d'tor

    // - only the following enums should be used for the following
    //   accessors
    class TSRegionDetail {
    public:
      enum enum_type  { unknown, TS1, TS2, TS3, TS4, TS5 };
      static std::string const& typeName() {
        static std::string type("TSRegion"); return type;
      }
      static std::map<enum_type,std::string> const& names() {
        static std::map<enum_type,std::string> nam;
        
        if ( nam.empty() ) {
          nam[unknown] = "unknown";
          nam[TS1]     = "TS1";
          nam[TS2]     = "TS2";
          nam[TS3]     = "TS3";
          nam[TS4]     = "TS4";
          nam[TS5]     = "TS5";
        }

        return nam;
      }
    };

    class TSRadialPartDetail {
    public:
      enum enum_type  { unknown, IN, OUT };
      static std::string const& typeName() {
        static std::string type("TSRadialPart"); return type;
      }
      static std::map<enum_type,std::string> const& names() {
        static std::map<enum_type,std::string> nam;

        if ( nam.empty() ) {
          nam[unknown] = "unknown";
          nam[IN]      = "IN";
          nam[OUT]     = "OUT";
        }
        
        return nam;
      }
    };

    // Coil Assemblies
    class TSCARegionDetail {
    public:
      enum enum_type  { unknown, TS1, TS2, TS3u, TS3uu, TS3dd, TS3d, TS4, TS5 };
      static std::string const& typeName() {
        static std::string type("TSCARegion"); return type;
      }
      static std::map<enum_type,std::string> const& names() {
        static std::map<enum_type,std::string> nam;
        
        if ( nam.empty() ) {
          nam[unknown] = "unknown";
          nam[TS1]     = "TS1";
          nam[TS2]     = "TS2";
          nam[TS3u]    = "TS3u";
	  //          nam[TS3ud]   = "TS3ud";
          nam[TS3uu]   = "TS3uu";
          nam[TS3dd]   = "TS3dd";
          nam[TS3d]    = "TS3d";
          nam[TS4]     = "TS4";
          nam[TS5]     = "TS5";
        }

        return nam;
      }
    };

    typedef EnumToStringSparse<TSRegionDetail> TSRegion;
    typedef EnumToStringSparse<TSRadialPartDetail> TSRadialPart;
    typedef EnumToStringSparse<TSCARegionDetail> TSCARegion;

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

    std::string downstreamVacuumMaterial() const { return _downstreamVacuumMaterial; }
    std::string upstreamVacuumMaterial()   const { return _upstreamVacuumMaterial; }

    template <class T = TSSection>
    T* getTSCryo(TSRegion::enum_type i,TSRadialPart::enum_type j) const {
      return static_cast<T*>( _cryoMap.find(i)->second.find(j)->second.get() );        
    }

    // The coils assemblies are approximated by a torus and cones for now

    const std::vector<double> & caRadii(TSCARegion::enum_type i) const { 
      return _caRadiiMap.at(i);
    }
    // those only make sense for the toruses and cylinders
    double innerCARadius(TSCARegion::enum_type i) const { 
      return _caRadiiMap.at(i)[0]; 
    }
    double outerCARadius(TSCARegion::enum_type i) const { 
      return _caRadiiMap.at(i)[1]; 
    }

    // Rings  (David Norvil Brown, April 2015)
    double rInRingSide() const { return _rInRingSide; }
    double rOutRingSide() const { return _rOutRingSide; }
    double thickRingSide() const { return _thickRingSide; }
    double rInRing() const { return _rInRing; }
    double rOutRing() const { return _rOutRing; }
    double lengthRing() const { return _lengthRing; }
    std::string RingMaterial() const { return _RingMaterial; }
    std::vector<double> xRing() const { return _xRing; }
    std::vector<double> yRing() const { return _yRing; }
    std::vector<double> zRing() const { return _zRing; }
    std::vector<double> thetaRing() const { return _thetaRing; }
    // Coils
    std::string coil_material() const { return _coilMaterial; }
    unsigned getNCoils(TSRegion::enum_type i) const { return _nCoils.at(i); }
    const std::vector<Coil>& getTSCoils(TSRegion::enum_type i) const { return _coilMap.at(i); }

    // Collimators
    CollimatorTS1 const& getColl1()  const { return _coll1;  }
    CollimatorTS3 const& getColl31() const { return _coll31; }
    CollimatorTS3 const& getColl32() const { return _coll32; }
    CollimatorTS5 const& getColl51()  const { return _coll51;  }
    CollimatorTS5 const& getColl52()  const { return _coll52;  }
    CollimatorTS5 const& getColl53()  const { return _coll53;  }

    // Coil Assemblies
    const std::string & caMaterial() const { return _caMaterial; }
    template <class T = TSSection>
    T* getTSCA(TSCARegion::enum_type i) const {
      return static_cast<T*>(_caMap.at(i).get() );
    }

    // Vacua
    template <class T = TSSection>
    T* getTSVacuum(TSRegion::enum_type i) const {
      return static_cast<T*>(_vacuumMap.find(i)->second.get() );
    }

    // Poly-lining
    const TorusSection* getTSPolyLining( TSRegion::enum_type i) const {
      return _polyLiningMap.find(i)->second.get();
    }

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
    std::string _downstreamVacuumMaterial;
    std::string _upstreamVacuumMaterial;

    // cryostat map
    typedef std::map<TSRadialPart::enum_type,std::unique_ptr<TSSection>> map_unique_ptrs_TSSection;
    std::map<TSRegion::enum_type,map_unique_ptrs_TSSection> _cryoMap;

    // Rings 
    double _rInRingSide, _rOutRingSide, _thickRingSide;
    double _rInRing, _rOutRing, _lengthRing;
    std::string _RingMaterial;
    std::vector<double> _xRing;
    std::vector<double> _yRing;
    std::vector<double> _zRing;
    std::vector<double> _thetaRing;

    // Coils
    std::string _coilMaterial;
    std::map<TSRegion::enum_type,std::vector<Coil>> _coilMap;
    const std::map<TSRegion::enum_type,unsigned> _nCoils = {
      {TSRegion::TS1, 4},
      {TSRegion::TS2,18},
      {TSRegion::TS3,10},
      {TSRegion::TS4,18},
      {TSRegion::TS5, 6}
    };
    
    // Collimators
    CollimatorTS1 _coll1;
    CollimatorTS3 _coll31;
    CollimatorTS3 _coll32;
    CollimatorTS5 _coll51;
    CollimatorTS5 _coll52;
    CollimatorTS5 _coll53;

    // Vacuum map
    std::map<TSRegion::enum_type,std::unique_ptr<TSSection>> _vacuumMap;

    // TS Coil Assemblies (CA)

    std::string _caMaterial;

    // if a section is a torus or a cylinder the size is 2; 4 for cone; in out; upstream, downstream
    std::map<TSCARegion::enum_type, std::vector<double>> _caRadiiMap;

    // CA Coil Assemblies map
    std::map<TSCARegion::enum_type, std::unique_ptr<TSSection>> _caMap;

    // Poly-lining map
    std::map<TSRegion::enum_type,std::unique_ptr<TorusSection>> _polyLiningMap;

    PbarWindow _pbarWindow;   

  };

}
#endif /* BeamlineGeom_TransportSolenoid_hh */
