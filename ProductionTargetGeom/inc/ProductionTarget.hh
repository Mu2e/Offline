// Geometry of the production target. This also defines the proton beam direction.
//
// Andrei Gaponenko, 2011
// Update to v2 in 2018 by David Norvil Brown

#ifndef PRODUCTIONTARGET_HH
#define PRODUCTIONTARGET_HH

#include <map>

#include "canvas/Persistency/Common/Wrapper.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Offline/GeomPrimitives/inc/Polycone.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class ProductionTargetMaker;

  class ProductionTarget : virtual public Detector {
  public:

    // cylinder parameters
    std::string tier1TargetType() const {return _tier1TargetType;}
    int    version() const { return _version; }
    double rOut() const { return _rOut; }
    double halfLength() const { return _halfLength; }
    double envHalfLength() const { return _envelHalfLength; }

    // Parameters for version 2 for stiffening parts
    double nFins       () const   { return _nFins;    }
    double finHeight   () const   { return _finHeight;    }
    double finThickness() const   { return _finThickness; }
    double hubDistUS() const      { return _hubDistUS;    }
    double hubDistDS() const      { return _hubDistDS;    }
    double hubAngleUS() const     { return _hubAngleUS;   }
    double hubAngleDS() const     { return _hubAngleDS;   }
    double hubOverhangUS() const  { return _hubOverhangUS;}
    double hubOverhangDS() const  { return _hubOverhangDS;}
    double hubLenUS() const       { return _hubLenUS;     }
    double hubLenDS() const       { return _hubLenDS;     }
    void   setHubLenUS(double& aLen ) { _hubLenUS = aLen; }
    void   setHubLenDS(double& aLen ) { _hubLenDS = aLen; }


    // in mu2e coordinates
    const CLHEP::Hep3Vector& position() const { return _prodTargetPosition; }

    // this is used to transorm particle momentum and position from
    // the PrimaryProtonGun frame to the Mu2e frame
    const CLHEP::HepRotation& protonBeamRotation() const { return _protonBeamRotation; }

    const CLHEP::Hep3Vector& haymanPosition() const { return _haymanProdTargetPosition; }



    // "passive" rotation, used for placing the production target.
    // This is the inverse of protonBeamRotation.
    const CLHEP::HepRotation& productionTargetRotation() const { return _protonBeamInverseRotation; }

    const Polycone * getHubsRgtPtr() const {return _pHubsRgtParams/*.get()*/;}
    const Polycone * getHubsLftPtr() const {return _pHubsLftParams/*.get()*/;}
    const std::map<double, CLHEP::Hep3Vector> & anchoringPntsRgt() const { return _anchoringPntsRgt; }
    const std::map<double, CLHEP::Hep3Vector> & anchoringPntsLft() const { return _anchoringPntsLft; }

    ~ProductionTarget() override {

      if (!_tier1TargetType.empty()){
        delete _pHubsRgtParams;
        delete _pHubsLftParams;
      }
    }
    ProductionTarget( ProductionTarget const&  ) = default;
    ProductionTarget( ProductionTarget&&       ) = default;
    ProductionTarget& operator=(ProductionTarget const&  ) = default;
    ProductionTarget& operator=(ProductionTarget &&      ) = default;


    //
    // accessors for hayman v2
    std::string haymanTargetType()           const  {return _haymanTargetType;}
    double halfHaymanLength()                const  {return _halfHaymanLength;}
    std::string targetCoreMaterial()         const  {return _targetCoreMaterial;}
    std::string targetFinMaterial()          const  {return _targetFinMaterial;}
    std::string targetVacuumMaterial()       const  {return _targetVacuumMaterial;}
    std::string supportRingMaterial()        const  {return _supportRingMaterial;}
    std::string spokeMaterial()              const  {return _spokeMaterial;}

    double rotHaymanX()                      const  {return _rotHaymanX;}
    double rotHaymanY()                      const  {return _rotHaymanY;}
    double rotHaymanZ()                      const  {return _rotHaymanZ;}
    CLHEP::Hep3Vector haymanProdTargetPosition()           const {return _haymanProdTargetPosition; }
    int numberOfTargetSections()                           const {return _numberOfTargetSections;}
    std::vector<double> startingSectionThickness()         const {return _startingSectionThickness;}
    std::vector<int> numberOfSegmentsPerSection()          const {return _numberOfSegmentsPerSection;}
    std::vector<double> thicknessOfSegmentPerSection()     const {return _thicknessOfSegmentPerSection;}
    std::vector<double> heightOfRectangularGapPerSection() const {return _heightOfRectangularGapPerSection;}
    std::vector<double> thicknessOfGapPerSection()         const {return _thicknessOfGapPerSection;}
    int nHaymanFins()                       const {return _nHaymanFins;}
    std::vector<double> finAngles()         const {return _finAngles;}
    double haymanFinThickness()                   const {return _haymanFinThickness;}
    double finOuterRadius()                 const {return _finOuterRadius;}
    double supportRingLength()              const {return _supportRingLength;}
    double supportRingInnerRadius()         const {return _supportRingInnerRadius;}
    double supportRingOuterRadius()         const {return _supportRingOuterRadius;}
    double supportRingCutoutThickness()     const {return _supportRingCutoutThickness;}
    double supportRingCutoutLength()        const {return _supportRingCutoutLength;}

    //define the supports for the Hayman_v2
    bool                       supportsBuild       () const {return _supportsBuild       ;}
    double                     supportWheelRIn     () const {return _supportWheelRIn     ;}
    double                     supportWheelROut    () const {return _supportWheelROut    ;}
    double                     supportWheelHL      () const {return _supportWheelHL      ;}
    std::string                supportWheelMaterial() const {return _supportWheelMaterial;}
    int                        nSpokesPerSide      () const {return _nSpokesPerSide      ;}
    //features on the bicycle wheel
    const std::vector<double>& supportWheelFeatureAngles() const {return _supportWheelFeatureAngles;}
    const std::vector<double>& supportWheelFeatureArcs  () const {return _supportWheelFeatureArcs  ;}
    const std::vector<double>& supportWheelFeatureRIns  () const {return _supportWheelFeatureRIns  ;}
    //rods in the support wheel
    const std::vector<double>& supportWheelRodHL          () const {return _supportWheelRodHL          ;}
    const std::vector<double>& supportWheelRodOffset      () const {return _supportWheelRodOffset      ;}
    const std::vector<double>& supportWheelRodRadius      () const {return _supportWheelRodRadius      ;}
    const std::vector<double>& supportWheelRodRadialOffset() const {return _supportWheelRodRadialOffset;}
    const std::vector<double>& supportWheelRodWireOffsetD () const {return _supportWheelRodWireOffsetD ;}
    const std::vector<double>& supportWheelRodWireOffsetU () const {return _supportWheelRodWireOffsetU ;}
    const std::vector<double>& supportWheelRodAngles      () const {return _supportWheelRodAngles      ;}

    //spoke parameters
    const std::vector<double>& spokeTargetAnglesD() const {return _spokeTargetAnglesD;}
    const std::vector<double>& spokeTargetAnglesU() const {return _spokeTargetAnglesU;}
    double                     spokeRadius       () const {return _spokeRadius       ;}

    double productionTargetMotherOuterRadius() const {return _productionTargetMotherOuterRadius;}
    double productionTargetMotherHalfLength()  const {return _productionTargetMotherHalfLength;}

    std::string hayman_v_2_0 = "Hayman_v_2_0";
    std::string tier1 = "MDC2018";

    CLHEP::Hep3Vector targetPositionByVersion() const {
      if (_haymanTargetType == hayman_v_2_0){
        return _haymanProdTargetPosition;}
      else if  (_tier1TargetType == "MDC2018"){
        return _prodTargetPosition;}
      else throw cet::exception("BADCONFIG")
             << "in ProductionTarget.hh, no valid target specified"<< std::endl;
    }
    double targetHalfLengthByVersion() const {
     if (_haymanTargetType == hayman_v_2_0){
        return _halfHaymanLength;}
     else if  (_tier1TargetType == "MDC2018"){
       return _halfLength;}
     else throw cet::exception("BADCONFIG")
            << "in ProductionTarget.hh, no valid target specified"<< std::endl;
    }



    //----------------------------------------------------------------

  private:
    friend class ProductionTargetMaker;

    // Private ctr: the class should be only obtained via ProductionTargetFNAL::ProductionTargetMaker.
    ProductionTarget(std::string tier1TargetType, int version, double rOut, double halfLength, double rotX,
                     double rotY, const CLHEP::Hep3Vector& position,
                     int    nFins,
                     double finHeight, double finThickness,
                     double hubDistUS, double hubDistDS,
                     double hubAngleUS, double hubAngleDS,
                     double hubOverhangUS, double hubOverhangDS );

    ProductionTarget(
                     std::string haymanTargetType, int version
                     ,double productionTargetMotherOuterRadius
                     ,double productionTargetMotherHalfLength
                     ,double rOut
                     ,double halfHaymanLength
                     ,double rotHaymanX
                     ,double rotHaymanY
                     ,double rotHaymanZ
                     ,const CLHEP::Hep3Vector& haymanProdTargetPosition
                     ,std::string targetCoreMaterial
                     ,std::string targetFinMaterial
                     ,std::string targetVacuumMaterial
                     ,std::string supportRingMaterial
                     ,std::string spokeMaterial
                     ,int numberOfTargetSections
                     ,std::vector<double> startingSectionThickness
                     ,std::vector<int> numberOfSegmentsPerSection
                     ,std::vector<double> thicknessOfSegmentPerSection
                     ,std::vector<double> heightOfRectangularGapPerSection
                     ,std::vector<double> thicknessOfGapPerSection
                     ,int nHaymanFins
                     ,std::vector<double> finAngles
                     ,double haymanFinThickness
                     ,double finOuterRadius
                     ,double supportRingLength
                     ,double supportRingInnerRadius
                     ,double supportRingOuterRadius
                     ,double supportRingCutoutThickness
                     ,double supportRingCutoutLength
                     );

    CLHEP::HepRotation _protonBeamRotation;

    // can't return by const ref if invert on the fly so need to store redundant data
    CLHEP::HepRotation _protonBeamInverseRotation;// FIXME: should be transient

//    std::unique_ptr<Polycone> _pHubsRgtParams;
//    std::unique_ptr<Polycone> _pHubsLftParams;
    // Fixme: replace with unique_ptr.
    Polycone * _pHubsRgtParams;
    Polycone * _pHubsLftParams;
    std::map<double,CLHEP::Hep3Vector> _anchoringPntsRgt;
    std::map<double,CLHEP::Hep3Vector> _anchoringPntsLft;

    CLHEP::Hep3Vector _prodTargetPosition;

    std::string _tier1TargetType;
    std::string _haymanTargetType;

    int    _version;
    double _productionTargetMotherOuterRadius;
    double _productionTargetMotherHalfLength;
    double _rOut;
    double _halfLength;
    double _envelHalfLength;

    // version 1+ parameters
    int    _nFins;
    double _finHeight;
    double _finThickness;
    double _hubDistUS;
    double _hubDistDS;
    double _hubAngleUS;
    double _hubAngleDS;
    double _hubOverhangUS;
    double _hubOverhangDS;
    double _hubLenUS;
    double _hubLenDS;


    // hayman parameters


    double _halfHaymanLength;
    double _rotHaymanX;
    double _rotHaymanY;
    double _rotHaymanZ;
    CLHEP::Hep3Vector _haymanProdTargetPosition;

    std::string _targetCoreMaterial;
    std::string _targetFinMaterial;
    std::string _targetVacuumMaterial;
    std::string _supportRingMaterial;
    std::string _spokeMaterial;
    int _numberOfTargetSections;
    std::vector<double> _startingSectionThickness;
    std::vector<int> _numberOfSegmentsPerSection;
    std::vector<double> _thicknessOfSegmentPerSection;
    std::vector<double> _heightOfRectangularGapPerSection;
    std::vector<double> _thicknessOfGapPerSection;
    int _nHaymanFins;
    std::vector<double> _finAngles;
    double _haymanFinThickness;
    double _finOuterRadius;
    double _supportRingLength;
    double _supportRingInnerRadius;
    double _supportRingOuterRadius;
    double _supportRingCutoutThickness;
    double _supportRingCutoutLength;

    //parameters for the support wheel
    bool                _supportsBuild       ; //whether or not to build the supports
    double              _supportWheelRIn     ;
    double              _supportWheelROut    ;
    double              _supportWheelHL      ; //half thickness in z
    std::string         _supportWheelMaterial;

    //parameters for the non-wheel features near the support rods that are on the wheel
    std::vector<double> _supportWheelFeatureAngles;
    std::vector<double> _supportWheelFeatureArcs  ;
    std::vector<double> _supportWheelFeatureRIns  ;

    //parameters for rods in the support wheel
    std::vector<double> _supportWheelRodHL          ; //includes length through the wheel
    std::vector<double> _supportWheelRodOffset      ; //z offset with respect to the wheel
    std::vector<double> _supportWheelRodRadius      ; //radius of the rod
    std::vector<double> _supportWheelRodRadialOffset; //radius from the wheel center the rod is centered at
    std::vector<double> _supportWheelRodWireOffsetD ; //z offset from the end of the rod the wire connects (downstream)
    std::vector<double> _supportWheelRodWireOffsetU ; //z offset from the end of the rod the wire connects (upstream)
    std::vector<double> _supportWheelRodAngles      ; //angle about the wheel the rod is

    //parameters for the wires (spokes) connecting the support wheel and the target
    // std::string         _spokeMaterial; (defined above)
    int                 _nSpokesPerSide    ; //also constrains the number of rods in the support wheel
    std::vector<double> _spokeTargetAnglesD; //angle about the target the wire connects to (downstream)
    std::vector<double> _spokeTargetAnglesU; //angle about the target the wire connects to (upstream)
    double              _spokeRadius       ; //radius of the wire

    // Needed for persistency
    template<class T> friend class art::Wrapper;
    ProductionTarget():_pHubsRgtParams(NULL), _pHubsLftParams(NULL) {}
  };
}

#endif/*PRODUCTIONTARGET_HH*/
