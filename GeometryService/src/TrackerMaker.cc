//
// Construct and return a Tracker.
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/GeometryService/inc/TrackerMaker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeneralUtilities/inc/HepTransform.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <bitset>
#include <string>
#include <sstream>

using cet::square;
using cet::diff_of_squares;
using cet::sum_of_squares;

using namespace std;

namespace mu2e {
  using StrawCollection = std::array<Straw,StrawId::_nustraws>;

  // Constructor that gets information from the config file instead of
  // from arguments.
  TrackerMaker::TrackerMaker( SimpleConfig const& config){
    parseConfig(config);

    buildIt( );

    // print straw layout for debbugging pourposes

    if (_verbosityLevel>2) {

      int ipln = -1;
      int ipnl = -1;
      int ilay = -1;

      cout << __func__ << " (_tt->_allStraws).size(), StrawId::_nustraws "
           << fixed << setw(6) << _tt->nStraws()
           << fixed << setw(6) << StrawId::_nustraws
           << endl;

      for ( const Straw& straw : _tt->getStraws() ){

        int cpln = straw.id().getPlane();
        int cpnl = straw.id().getPanel();
        int clay = straw.id().getLayer();

        const Plane& plane = _tt->getPlane(cpln);
        const Panel& panel = plane.getPanel(cpnl);

        size_t nStrawsPerPanel = panel.nStraws();
        size_t nStrawsPerPlane = plane.nPanels() * nStrawsPerPanel;


        size_t ipnlf = nStrawsPerPanel*cpnl + nStrawsPerPlane*cpln;

        cout << __func__ << " Straw "
             << fixed << " "
             << " plnfloor " << setw(6) << ipnlf << " "
             << straw.id() << " "
             << " mid point " << straw.getMidPoint()
             << " r " << sqrt(straw.getMidPoint()[0]*straw.getMidPoint()[0]+
                              straw.getMidPoint()[1]*straw.getMidPoint()[1])
             << " direction " << straw.getDirection()
             << " straw0MidPoint  " << panel.straw0MidPoint()
             << " straw0Direction " << panel.straw0Direction()
             << " origin " << plane.origin();

        if (ipnl>cpnl && ipln==cpln) cout << " <--S";
        if (ilay>clay && ipnl==cpnl) cout << " <--L";
        if (ipln!=cpln) ipln=cpln;
        if (ipnl!=cpnl) ipnl=cpnl;
        if (ilay!=clay) ilay=clay;

        cout << endl;

      } // end straw loop

    } // end if verbosity

  }

  void TrackerMaker::parseConfig( const SimpleConfig& config ){

    _verbosityLevel     = config.getInt("tracker.verbosityLevel",0);
    _ttVersion          = config.getInt("TrackerVersion",3);

    _motherRIn        = config.getDouble("tracker.mother.rIn"        )*CLHEP::mm;
    _motherROut       = config.getDouble("tracker.mother.rOut"       )*CLHEP::mm;
    _motherHalfLength = config.getDouble("tracker.mother.halfLength" )*CLHEP::mm;
    _motherZ0         = config.getDouble("tracker.mother.z0"         )*CLHEP::mm;
    _manifoldsPerEnd    = config.getInt("tracker.manifoldsPerEnd");
    _strawsPerManifold  = config.getInt("tracker.strawsPerManifold");
    _rotationPattern    = config.getInt("tracker.rotationPattern");
    _panelZPattern      = config.getInt("tracker.panelZPattern");
    _layerZPattern      = config.getInt("tracker.layerZPattern");
    _spacingPattern     = config.getInt("tracker.spacingPattern");
    _innermostLayer     = config.getInt("tracker.innermostLayer");

    _planePadding      = config.getDouble("tracker.planePadding",0.5)*CLHEP::mm;
    _panelPadding      = config.getDouble("tracker.panelPadding",0.25)*CLHEP::mm;

    _oddStationRotation   =  config.getDouble("tracker.oddStationRotation")*CLHEP::degree;
    _zCenter              =  config.getDouble("tracker.z0")*CLHEP::mm;
    _xCenter              = -config.getDouble("mu2e.solenoidOffset")*CLHEP::mm;
    _envelopeInnerRadius  =  config.getDouble("tracker.envelopeInnerRadius")*CLHEP::mm;
    _rInnermostWire       =  config.getDouble("tracker.rInnermostWire")*CLHEP::mm;
    _strawOuterRadius     =  config.getDouble("tracker.strawOuterRadius")*CLHEP::mm;
    _strawWallThickness   =  config.getDouble("tracker.strawWallThickness")*CLHEP::mm;
    _strawGap             =  config.getDouble("tracker.strawGap")*CLHEP::mm;
    _planeSpacing        =  config.getDouble("tracker.planeSpacing")*CLHEP::mm; // this is really the station spacing
    _planeHalfSeparation =  config.getDouble("tracker.planeHalfSeparation")*CLHEP::mm;

    _outerSupportRadius   =  config.getDouble("tracker.outerSupportRadius")*CLHEP::mm;
    _innerSupportRadius   =  config.getDouble("tracker.innerSupportRadius")*CLHEP::mm;
    _supportHalfThickness =  config.getDouble("tracker.supportHalfThickness")*CLHEP::mm;
    _wireRadius           =  config.getDouble("tracker.wireRadius")*CLHEP::mm;
    _virtualDetectorHalfLength = config.getDouble("vd.halfLength")*CLHEP::mm;

    config.getVectorInt("tracker.nonExistingPlanes", _nonExistingPlanes,  vector<int>() );

    _verbosityLevel > 0 && _nonExistingPlanes.size()>0 &&
      cout << __func__ << " inactive planes : f/l   "
           << _nonExistingPlanes.front() << " / "
           << _nonExistingPlanes.back()
           << endl;

    config.getVectorDouble("tracker.manifoldHalfLengths", _manifoldHalfLengths, 3);
    for ( size_t i=0; i<_manifoldHalfLengths.size(); ++i ){
      _manifoldHalfLengths.at(i) *= CLHEP::mm;
    }

    config.getVectorString("tracker.strawMaterials", _strawMaterials, 3);

    _envelopeMaterial = config.getString("tracker.mat.vacuum");
    _supportMaterial = config.getString("tracker.mat.support");

    _passivationMargin        =  config.getDouble("tracker.passivationMargin")*CLHEP::mm;

    // For the detailv0 support model there are lots of extra parameters.
    _supportModel = SupportModel( config.getString("trackerSupport.model","simple"));
    if ( _supportModel == SupportModel::detailedv0 ) {
      _endRingOuterRadius      = config.getDouble( "trackerSupport.endRing.outerRadius" );
      _endRingInnerRadius      = config.getDouble( "trackerSupport.endRing.innerRadius" );
      _endRingHalfLength       = config.getDouble( "trackerSupport.endRing.halfLength"  );
      _endRingZOffset          = config.getDouble( "trackerSupport.endRing.zOffset"     );
      _endRingMaterial         = config.getString( "trackerSupport.endRing.material"    );

      _hasDownRing    = config.getBool( "trackerSupport.downRing.build",false);
      if ( _hasDownRing ) {
        _downRingOuterRadius      = config.getDouble( "trackerSupport.downRing.outerRadius" );
        _downRingInnerRadius      = config.getDouble( "trackerSupport.downRing.innerRadius" );
        _downRingHalfLength       = config.getDouble( "trackerSupport.downRing.halfLength"  );
        _downRingZOffset          = config.getDouble( "trackerSupport.downRing.zOffset"     );
        _downRingMaterial         = config.getString( "trackerSupport.downRing.material"    );
      }

      config.getVectorInt( "trackerSupport.midRing.slot", _midRingSlot );
      _midRingHalfLength       = config.getDouble(    "trackerSupport.midRing.halfLength" );
      _midRingPhi0 = config.getDouble( "trackerSupport.midRing.Phi0",180.0)*CLHEP::degree;
      _midRingdPhi = config.getDouble( "trackerSupport.midRing.dPhi",180.0)*CLHEP::degree;
      _midRingMaterial = config.getString( "trackerSupport.midRing.material", "StainlessSteel316");
      // support beams;
      // fixme use vectors to contain them all (e.g. vector<SupportBeamParams>)

      config.getVectorDouble( "trackerSupport.beam0.phiRange", _beam0_phiRange );
      _beam0_innerRadius     = config.getDouble( "trackerSupport.beam0.innerRadius" );
      _beam0_outerRadius     = config.getDouble( "trackerSupport.beam0.outerRadius" );
      _beam0_material        = config.getString( "trackerSupport.beam0.material" );

      config.getVectorDouble( "trackerSupport.beam1.phiRange", _beam1_phiRange );
      config.getVectorDouble( "trackerSupport.beam1.phiSpans", _beam1_phiSpans );
      config.getVectorDouble( "trackerSupport.beam1.servicePhi0s", _beam1_servicePhi0s );
      config.getVectorDouble( "trackerSupport.beam1.servicePhiEnds", _beam1_servicePhiEnds );
      _beam1_innerRadius = config.getDouble( "trackerSupport.beam1.innerRadius" );
      _beam1_midRadius1  = config.getDouble( "trackerSupport.beam1.midRadius1" );
      _beam1_midRadius2  = config.getDouble( "trackerSupport.beam1.midRadius2" );
      _beam1_outerRadius = config.getDouble( "trackerSupport.beam1.outerRadius" );
      _beam1_material    = config.getString( "trackerSupport.beam1.material" );
      config.getVectorDouble( "trackerSupport.beam1.serviceOuterRadii", _beam1_serviceOuterRadii );
      config.getVectorString( "trackerSupport.beam1.serviceMaterials", _beam1_serviceMaterials );
      config.getVectorDouble( "trackerSupport.beam1.serviceCovRelThickness", _beam1_serviceCovRelThickness );
      config.getVectorString( "trackerSupport.beam1.serviceMaterialsCov", _beam1_serviceMaterialsCov );


      _panelPhi  = config.getDouble("trackerSupport.phiCoverage",120.0)*CLHEP::degree;
      _dphiRibs = config.getDouble("trackerSupport.dphiRibs",27.0)*CLHEP::degree;
      _ribHalfAngle = config.getDouble("trackerSupport.ribHalfAngle",1.0)*CLHEP::degree;

      _innerRingInnerRadius    = config.getDouble( "trackerSupport.innerRing.innerRadius" );
      _innerRingOuterRadius    = config.getDouble( "trackerSupport.innerRing.outerRadius" );
      _innerRingHalfLength     = config.getDouble( "trackerSupport.innerRing.halfLength"  );
      _innerRingMaterial       = config.getString( "trackerSupport.innerRing.material"    );

      _centerPlateHalfLength   = config.getDouble( "trackerSupport.centerPlate.halfLength" );
      _centerPlateMaterial     = config.getString( "trackerSupport.centerPlate.material"   );

      _outerRingInnerRadius    = config.getDouble( "trackerSupport.outerRing.innerRadius" );
      _outerRingOuterRadius    = config.getDouble( "trackerSupport.outerRing.outerRadius" );
      _outerRingMaterial       = config.getString( "trackerSupport.outerRing.material"    );

      _coverHalfLength          = config.getDouble( "trackerSupport.cover.halfLength"           );
      _coverMaterial            = config.getString( "trackerSupport.cover.material"             );
      _electronicsG10HalfLength = config.getDouble( "trackerSupport.electronics.g10.halfLength" );
      _electronicsG10Material   = config.getString( "trackerSupport.electronics.g10.material"   );
      _electronicsCuHhalfLength = config.getDouble( "trackerSupport.electronics.cu.halfLength"  );
      _electronicsCuMaterial    = config.getString( "trackerSupport.electronics.cu.material"    );
      _channelZOffset           = config.getDouble( "trackerSupport.channel.zOffset"            );
      _panelZOffset             = config.getDouble( "trackerSupport.panel.zOffset", 0.0 );
      _channelDepth             = config.getDouble( "trackerSupport.channel.depth"              );
      _channelMaterial          = config.getString( "trackerSupport.channel.material"           );
      _electronicsSpaceMaterial = config.getString( "trackerSupport.electronicsSpace.material"  );

      const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
      geomOptions->loadEntry( config, "trackerElectronicsKey",   "trackerSupport.electronics.key");
      geomOptions->loadEntry( config, "trackerElectronicsShield","trackerSupport.electronics.key.shield");
      _EBKeyVisible          = geomOptions->isVisible("trackerElectronicsKey");
      _EBKeySolid            = geomOptions->isSolid("trackerElectronicsKey");
      _EBKeyShieldVisible    = geomOptions->isVisible("trackerElectronicsShield");
      _EBKeyShieldSolid      = geomOptions->isSolid("trackerElectronicsShield");

      _EBKeyHalfLength         = config.getDouble("trackerSupport.electronics.key.halfLength");
      _EBKeyShieldHalfLength   = config.getDouble("trackerSupport.electronics.key.shieldHalfLength");
      _EBKeyInnerRadius        = config.getDouble("trackerSupport.electronics.key.innerRadius");
      _EBKeyOuterRadius        = config.getDouble("trackerSupport.electronics.key.outerRadius");
      _EBKeyShiftFromPanelFace = config.getDouble("trackerSupport.electronics.key.shiftFromPanelFace");
      //_EBKeyVisible            = config.getBool(  "trackerSupport.electronics.key.visible");
      //_EBKeySolid              = config.getBool(  "trackerSupport.electronics.key.solid");
      //_EBKeyShieldVisible      = config.getBool(  "trackerSupport.electronics.key.shieldVisible");
      //_EBKeyShieldSolid        = config.getBool(  "trackerSupport.electronics.key.shieldSolid");
      _EBKeyMaterial           = config.getString("trackerSupport.electronics.key.material");
      _EBKeyShieldMaterial     = config.getString("trackerSupport.electronics.key.shieldMaterial");
      _EBKeyPhiRange           = config.getDouble("trackerSupport.electronics.key.phiRange")*CLHEP::degree;
      _EBKeyPhiExtraRotation   = config.getDouble("trackerSupport.electronics.key.phiExtraRotation")*CLHEP::degree;

      _wallOuterMetalThickness  = config.getDouble("tracker.straw.wallOuterMetal.thickness")*CLHEP::mm;
      _wallInnerMetal1Thickness = config.getDouble("tracker.straw.wallInnerMetal1.thickness")*CLHEP::mm;
      _wallInnerMetal2Thickness = config.getDouble("tracker.straw.wallInnerMetal2.thickness")*CLHEP::mm;
      _wirePlateThickness       = config.getDouble("tracker.straw.wirePlate.thickness")*CLHEP::mm;
      _wallOuterMetalMaterial   = config.getString("tracker.straw.wallOuterMetal.material");
      _wallInnerMetal1Material  = config.getString("tracker.straw.wallInnerMetal1.material");
      _wallInnerMetal2Material  = config.getString("tracker.straw.wallInnerMetal2.material");
      _wirePlateMaterial        = config.getString("tracker.straw.wirePlate.material");

    }

    //string tracker.mat.manifold  = "G4_Al";  // Placeholder.

    // Also define some parameters that may become variable some day.
    _panelBaseRotations.clear();
    _panelZSide.clear();

    if (StrawId::_npanels == 6 ){
  // the Z pattern here is forced by the 'alternating panel id' convention
  // that permeates the rest of the code.
        _panelZSide.push_back(-1.);
        _panelZSide.push_back(+1.);
        _panelZSide.push_back(-1.);
        _panelZSide.push_back(+1.);
        _panelZSide.push_back(-1.);
        _panelZSide.push_back(+1.);
      if(_rotationPattern == 1){
// cdr geometry, taken from DOC 888, also alternatives 1 and 3 from doc 2799
        // faces overlap by 60 degrees
// Implicitly define the rotations for the even and odd (sequentially) panels.
        _panelBaseRotations.push_back(  45.*CLHEP::degree);
        _panelBaseRotations.push_back( 105.*CLHEP::degree);
        _panelBaseRotations.push_back( 165.*CLHEP::degree);
        _panelBaseRotations.push_back( 225.*CLHEP::degree);
        _panelBaseRotations.push_back( 285.*CLHEP::degree);
        _panelBaseRotations.push_back( 345.*CLHEP::degree);
        _panelBaseRotations.push_back(  75.*CLHEP::degree);
        _panelBaseRotations.push_back(  15.*CLHEP::degree);
        _panelBaseRotations.push_back( 195.*CLHEP::degree);
        _panelBaseRotations.push_back( 135.*CLHEP::degree);
        _panelBaseRotations.push_back( 315.*CLHEP::degree);
        _panelBaseRotations.push_back( 255.*CLHEP::degree);

      } else if(_rotationPattern==2){
        // alternative 2 from DOC 2799
        // faces overlap by 60 degrees
        _panelBaseRotations.push_back(  45.*CLHEP::degree);
        _panelBaseRotations.push_back(  75.*CLHEP::degree);
        _panelBaseRotations.push_back( 165.*CLHEP::degree);
        _panelBaseRotations.push_back( 195.*CLHEP::degree);
        _panelBaseRotations.push_back( 285.*CLHEP::degree);
        _panelBaseRotations.push_back( 315.*CLHEP::degree);
        _panelBaseRotations.push_back( 105.*CLHEP::degree);
        _panelBaseRotations.push_back(  15.*CLHEP::degree);
        _panelBaseRotations.push_back( 225.*CLHEP::degree);
        _panelBaseRotations.push_back( 135.*CLHEP::degree);
        _panelBaseRotations.push_back( 345.*CLHEP::degree);
        _panelBaseRotations.push_back( 255.*CLHEP::degree);

      } else if(_rotationPattern==3){
      // faces overlap by 60 degrees, second plane 'flipped'
        _panelBaseRotations.push_back(   0.*CLHEP::degree);
        _panelBaseRotations.push_back(  90.*CLHEP::degree);
        _panelBaseRotations.push_back(  120.*CLHEP::degree);
        _panelBaseRotations.push_back(  210.*CLHEP::degree);
        _panelBaseRotations.push_back(  240.*CLHEP::degree);
        _panelBaseRotations.push_back(  330.*CLHEP::degree);
        _panelBaseRotations.push_back(  30.*CLHEP::degree);
        _panelBaseRotations.push_back(  60.*CLHEP::degree);
        _panelBaseRotations.push_back(  150.*CLHEP::degree);
        _panelBaseRotations.push_back(  180.*CLHEP::degree);
        _panelBaseRotations.push_back(  270.*CLHEP::degree);
        _panelBaseRotations.push_back(  300.*CLHEP::degree);
      } else if(_rotationPattern==4){
//-----------------------------------------------------------------------------
// Mu2e-2 studies: 2 faces within one plane have parallel straws within
//                 each 120 deg panel
//-----------------------------------------------------------------------------
        _panelBaseRotations.push_back(   0.*CLHEP::degree);
        _panelBaseRotations.push_back(   0.*CLHEP::degree);
        _panelBaseRotations.push_back(  120.*CLHEP::degree);
        _panelBaseRotations.push_back(  120.*CLHEP::degree);
        _panelBaseRotations.push_back(  240.*CLHEP::degree);
        _panelBaseRotations.push_back(  240.*CLHEP::degree);
        _panelBaseRotations.push_back(   60.*CLHEP::degree);
        _panelBaseRotations.push_back(   60.*CLHEP::degree);
        _panelBaseRotations.push_back(  210.*CLHEP::degree);
        _panelBaseRotations.push_back(  210.*CLHEP::degree);
        _panelBaseRotations.push_back(  270.*CLHEP::degree);
        _panelBaseRotations.push_back(  270.*CLHEP::degree);
      } else if(_rotationPattern==5){
        _panelBaseRotations.push_back(   15.*CLHEP::degree);
        _panelBaseRotations.push_back(  105.*CLHEP::degree);
        _panelBaseRotations.push_back(  135.*CLHEP::degree);
        _panelBaseRotations.push_back(  225.*CLHEP::degree);
        _panelBaseRotations.push_back(  255.*CLHEP::degree);
        _panelBaseRotations.push_back(  345.*CLHEP::degree);
        _panelBaseRotations.push_back(  165.*CLHEP::degree);
        _panelBaseRotations.push_back(   75.*CLHEP::degree);
        _panelBaseRotations.push_back(   45.*CLHEP::degree);
        _panelBaseRotations.push_back(  315.*CLHEP::degree);
        _panelBaseRotations.push_back(  285.*CLHEP::degree);
        _panelBaseRotations.push_back(  195.*CLHEP::degree);
        _panelBaseRotations.push_back(  165.*CLHEP::degree);
        _panelBaseRotations.push_back(   75.*CLHEP::degree);
        _panelBaseRotations.push_back(   45.*CLHEP::degree);
        _panelBaseRotations.push_back(  315.*CLHEP::degree);
        _panelBaseRotations.push_back(  285.*CLHEP::degree);
        _panelBaseRotations.push_back(  195.*CLHEP::degree);
        _panelBaseRotations.push_back(   15.*CLHEP::degree);
        _panelBaseRotations.push_back(  105.*CLHEP::degree);
        _panelBaseRotations.push_back(  135.*CLHEP::degree);
        _panelBaseRotations.push_back(  225.*CLHEP::degree);
        _panelBaseRotations.push_back(  255.*CLHEP::degree);
        _panelBaseRotations.push_back(  345.*CLHEP::degree);
       } else if(_rotationPattern==6){ // as-built planes, doc 888v13
        _panelBaseRotations.push_back(  105.*CLHEP::degree);
        _panelBaseRotations.push_back(  195.*CLHEP::degree);
        _panelBaseRotations.push_back(  225.*CLHEP::degree);
        _panelBaseRotations.push_back(  315.*CLHEP::degree);
        _panelBaseRotations.push_back(  345.*CLHEP::degree);
        _panelBaseRotations.push_back(  435.*CLHEP::degree);
        _panelBaseRotations.push_back(  255.*CLHEP::degree);
        _panelBaseRotations.push_back(  165.*CLHEP::degree);
        _panelBaseRotations.push_back(  135.*CLHEP::degree);
        _panelBaseRotations.push_back(  405.*CLHEP::degree);
        _panelBaseRotations.push_back(  375.*CLHEP::degree);
        _panelBaseRotations.push_back(  285.*CLHEP::degree);
        _panelBaseRotations.push_back(  255.*CLHEP::degree);
        _panelBaseRotations.push_back(  165.*CLHEP::degree);
        _panelBaseRotations.push_back(  135.*CLHEP::degree);
        _panelBaseRotations.push_back(  405.*CLHEP::degree);
        _panelBaseRotations.push_back(  375.*CLHEP::degree);
        _panelBaseRotations.push_back(  285.*CLHEP::degree);
        _panelBaseRotations.push_back(  105.*CLHEP::degree);
        _panelBaseRotations.push_back(  195.*CLHEP::degree);
        _panelBaseRotations.push_back(  225.*CLHEP::degree);
        _panelBaseRotations.push_back(  315.*CLHEP::degree);
        _panelBaseRotations.push_back(  345.*CLHEP::degree);
        _panelBaseRotations.push_back(  435.*CLHEP::degree);
      } else {
        throw cet::exception("GEOM")
          << "Unrecognized rotation pattern in TrackerMaker. \n";
      }
    } else {
      throw cet::exception("GEOM")
        << "Unrecognized rotation pattern in TrackerMaker. \n";
    }

    _layerHalfSpacing = 0.0;
    _layerHalfShift = 0.0;
    _manifoldXEdgeExcessSpace = 0.0;
    _manifoldZEdgeExcessSpace = 0.0;

  } // end TrackerMaker::parseConfig


  void plntest( const Plane& plane){
    cout << "Plane: "
         << plane.id() << " "
         << plane.origin()
         << endl;
  }


  void TrackerMaker::buildIt(){
  // old straw bookkeeping: still being used
    _strawTrckrConstrCount = -1; // first straw will be at 0

    // Straw Properties
    StrawProperties strawprops;
    strawprops._strawInnerRadius     = _strawOuterRadius - _strawWallThickness;
    strawprops._strawOuterRadius     = _strawOuterRadius;
    strawprops._strawWallThickness   = _strawWallThickness;
    strawprops._outerMetalThickness  = _wallOuterMetalThickness;
    strawprops._innerMetal1Thickness = _wallInnerMetal1Thickness;
    strawprops._innerMetal2Thickness = _wallInnerMetal2Thickness;
    strawprops._wireRadius           = _wireRadius;
    strawprops._wirePlateThickness   = _wirePlateThickness;
// straw
    StrawCollection allStraws;
    // Z location of the first plane.  This is used below
    _z0 = -findFirstPlaneZ0();
    computeLayerSpacingAndShift();
    computeManifoldEdgeExcessSpace();
    // fill their content
    makeStraws(allStraws);
    // create an empty TrackerG4Info: this gets filled further down
    auto g4trackerptr = shared_ptr<TrackerG4Info>(new TrackerG4Info);
    // see which planes actually exist.  This is deprecated, TrackerStatus should be used instead  FIXME!
    std::array<bool,StrawId::_nplanes> planeExists;
    for ( int ipln=0; ipln<StrawId::_nplanes; ++ipln ){
      planeExists[ipln] =  ( find ( _nonExistingPlanes.begin(), _nonExistingPlanes.end(), ipln) == _nonExistingPlanes.end() );
    }
    // build tracker with these
    _tt = unique_ptr<Tracker>(new Tracker(allStraws,strawprops,g4trackerptr,planeExists));
    auto g4tt = _tt->g4Tracker();

// now G4 stuff
    for ( int ipln=0; ipln<StrawId::_nplanes; ++ipln ){
      makePlane(StrawId(ipln,0,0));
    }

    makeMother();


    g4tt->_supportModel = _supportModel;

    // Fill information about the new style, fully detailed support structure.
    if ( _supportModel == SupportModel::detailedv0 ) {
      makeSupportStructure();
    }

    // Fill the information about the, old style minimal supports.
    g4tt->_supportParams = Support( _innerSupportRadius,
                                   _outerSupportRadius,
                                   _supportHalfThickness,
                                   _supportMaterial);

    g4tt->_z0                  = _zCenter;
    g4tt->_envelopeMaterial    = _envelopeMaterial;

    g4tt->_wallMaterialName    = _strawMaterials[0];
    g4tt->_outerMetalMaterial  = _wallOuterMetalMaterial;
    g4tt->_innerMetal1Material = _wallInnerMetal1Material;
    g4tt->_innerMetal2Material = _wallInnerMetal2Material;
    g4tt->_gasMaterialName     = _strawMaterials[1];
    g4tt->_wireMaterialName    = _strawMaterials[2];
    g4tt->_wirePlateMaterial   = _wirePlateMaterial;
// the following is needed by the Detail constructor
    g4tt->_panelZOffset        = _panelZOffset;


    //computeStrawHalfLengths();

    // Order is important here.
    computePlaneEnvelope();
    computeTrackerEnvelope();
    // This uses information from the planes
    makeThinSupportRings();


    finalCheck();

    if ( _verbosityLevel > 0 ) {
      cout << "Tracker Support Structure: \n" << g4tt->_supportStructure << endl;
    }

  } //end TrackerMaker::buildIt.

  void TrackerMaker::makeMother(){
    auto g4tt = _tt->g4Tracker();

    g4tt->_mother = PlacedTubs ( "TrackerMother",
                                TubsParams( _motherRIn, _motherROut, _motherHalfLength),
                                CLHEP::Hep3Vector( _xCenter, 0., _motherZ0),
                                _envelopeMaterial );

  } //end TrackerMaker::makeMother

  // In the present code the straw positions are computed using the manifold information.
  // The channel position is computed using the SupportStructure information.  These two
  // must be self consistent.  Relatively soon the manifold information will be removed
  // and they will be replaced with the support structure information.
  // For now, check for self-consistency.
  void TrackerMaker::finalCheck( ){

    // Only do this test for the new model(s).

    StrawId s0(0,1,0);
    StrawId s1(0,1,1);
    double ztest = 0.5*(_tt->getStraw ( s0 ).getMidPoint().z()+_tt->getStraw ( s1 ).getMidPoint().z())
      -  _tt->getPlane( StrawId(0) ).origin().z();
    ztest *= panelZSide(1,0);

    double tolerance = 1.e-6;
    if ( std::abs(ztest-_channelZOffset) > tolerance ){
      throw cet::exception("GEOM")  << "Inconsistent channel center and wire z location. \n"
                                    << "channel Z offset: " <<  _channelZOffset << "\n"
                                    << "plane z center:  " << _tt->getPlane( StrawId(0) ).origin().z() << "\n "
                                    << "Straw z layer 0:  " << _tt->getStraw ( s0 ).getMidPoint().z() << "\n"
                                    << "Straw z layer 1:  " << _tt->getStraw ( s1 ).getMidPoint().z() << "\n"
                                    << "z test:           " << ztest << " delta " << std::abs(ztest-_channelZOffset) << "\n";
    }

  }

  void TrackerMaker::makeStraws(StrawCollection& straws) {
    // loop over straws
    for ( int ipln=0; ipln<StrawId::_nplanes; ++ipln ){
      double planeDeltaZ = choosePlaneSpacing(ipln);
      CLHEP::Hep3Vector planeorigin( 0.0, 0.0, _z0+planeDeltaZ);
      for ( int ipnl=0; ipnl<StrawId::_npanels; ++ipnl ){
        _strawPanelConstrCount = -1; // these counts are useless and should go away FIXME!
        for ( int ilay=0; ilay<StrawId::_nlayers; ++ilay ){
          StrawId sid(ipln,ipnl,ilay);
          // we use the straw field to indicate the layer.  Layer structure is deprecated FIXME
          makeLayer(sid,planeorigin,straws);
        }
      }
    }
  }

  void TrackerMaker::makePlane( const StrawId& planeId ){
    //std::cout << "->->-> makePlane\n";
    int ipln = planeId.getPlane();
    auto const& planes = _tt->planes();
    auto const& plane = planes.at(ipln);

    if (_verbosityLevel>2) {
      cout << __func__ << " making plane " <<  ipln;
    }
    if (_verbosityLevel>2) {
      cout << ", exists " <<  _tt->planeExists(planeId) << endl;
    }
    for ( int ipnl=0; ipnl<StrawId::_npanels; ++ipnl ){
      makePanel ( StrawId(ipln,ipnl,0));
      if (_verbosityLevel>2) {
        size_t istr = -1; // local index in the panel
        Panel const& panel = *( plane.panels().at(ipnl) );
        std::ostringstream pnlid("",std::ios_base::ate); // to write at the end
        pnlid << panel.id();
        cout << __func__ << " straws in panel "
             << setw(7) << pnlid.str();
        // test of Plane, Panel functions
        cout << setw(4) << panel.id().getPlane();
        cout << setw(2) << panel.id().getPanel() << endl;

        for (const auto istr_p : panel.straws()) {
          if ( istr_p == nullptr ) continue;
          StrawId const & lsid =  (*istr_p).id();
          std::ostringstream nsid("",std::ios_base::ate); // to write at the end
          nsid << lsid;
          cout << setw(3) << ++istr
               << setw(6) << lsid.asUint16()
               << setw(17) << std::bitset<16>(lsid.asUint16())
               << " "
               << setw(6) << std::showbase << std::hex << lsid.asUint16()
               << " " << std::dec << std::noshowbase << setw(7) << nsid.str();
          nsid.str("");
          nsid << panel.getStraw(lsid).id();
          cout << setw(8) << nsid.str();
          nsid.str("");
          StrawId sid  = (*istr_p).id();
          nsid << panel.getStraw(sid).id(); // old id of the straw in the old container
          cout << setw(10) << nsid.str()
               << endl;
        } // straw loop
      } // verbosity check
    } // loop over panels

  }

  void TrackerMaker::makePanel( const PanelId& pnlId ){

    auto const& panel = _tt->getPlane(pnlId).getPanel(pnlId);
    auto g4tt = _tt->g4Tracker();

    // check if the opposite panels do not overlap
    static double const tolerance = 1.e-6; // this should be in a config file

    if ((2.*_manifoldHalfLengths.at(2)+_supportHalfThickness)>_planeHalfSeparation + tolerance) {
      cout << "(2.*_manifoldHalfLengths.at(2)+_supportHalfThickness), _planeHalfSeparation " <<
        (2.*_manifoldHalfLengths.at(2)+_supportHalfThickness) << ", " <<_planeHalfSeparation << endl;
      throw cet::exception("GEOM")  << "Planes are too close \n";
    }

    double strawSpacing = _strawGap+2.*_strawOuterRadius;

    if (_verbosityLevel>2) {
      cout << __func__ << " Expecting the straw spacing to be: " << strawSpacing << endl;
    }

    for ( int ilay=0; ilay<StrawId::_nlayers; ++ilay ){

      // checking spacing of the individual layers
      // are the manifolds sized correctly for the straws?

      // print out the Layer containers (or thier equivalence)

      if (_verbosityLevel>2) {
        uint16_t is = -1;
        for (const auto straw_p : panel.getStrawPointers() ) {
          const Straw& straw(*straw_p);
          ++is;
          if ( ( straw.id().getStraw())%2 != ilay ) continue;
          StrawId lid = straw.id().getLayerId();
          cout << __func__ << " Printing Layer _straws info: " << lid
               << setw(3) << is
               << " " << straw.id()
               << endl;
        } // loop over straws in the panel

      } // end verbosity

      for (  auto is = panel.getStrawPointers().cbegin();
             is < (panel.getStrawPointers().cend()-2); ++is) {
        const Straw& straw0(**is);
        uint16_t sn = straw0.id().getStraw();
        if ( sn%2 != ilay ) continue;
        StrawId lid = straw0.id().getLayerId();
        if (_verbosityLevel>2) {
          cout << __func__ << " Checking spacing"
               << " for layer " << lid << " straw " << straw0.id()  << endl;
        }
        const Straw& straw1(**(is+2));
        double layerDeltaMag =
          (straw1.getMidPoint() - straw0.getMidPoint()).mag();
        if ( abs(layerDeltaMag-strawSpacing)> tolerance ) {
          cout << "Layer straw spacing is (mm)   : " << layerDeltaMag
               << " for layer " << lid << " straw " << straw0.id()  << endl;
          cout << "It should be                  : " << strawSpacing << " diff: "
               << (layerDeltaMag-strawSpacing) << endl;

          throw cet::exception("GEOM")  << "Incorrect intralayer straw spacing, check manifold sizes etc..\n";

        } // end if tolerance check
      } // loop over straws in the panel
    } // loop over layers

    // check spacing between layers/straws

    if (StrawId::_nlayers>1) {

      for (  auto is = panel.getStrawPointers().cbegin();
             is < (panel.getStrawPointers().cend()-1); ++is) {
        const Straw& straw0(**is);
        uint16_t sn = straw0.id().getStraw();
        if ( sn%2 != 0 ) continue;
        const Straw& straw1(**(is+1));
        StrawId lid0 = straw0.id().getLayerId();
        StrawId lid1 = straw1.id().getLayerId();

        double xLayerDeltaMag =
          (straw0.getMidPoint() - straw1.getMidPoint()).mag();

        if (_verbosityLevel>2) {
          cout << __func__ << " Checking spacing"
               << " for layer " << lid0 << " straw " << straw0.id()
               << " and for layer " << lid1 << " straw " << straw1.id()  << endl;
        }

        if ( abs(xLayerDeltaMag-strawSpacing)> tolerance ) {
          cout << "xLayer straw spacing is (mm)   : "
               << xLayerDeltaMag
               << " for straws: "
               << straw0.id() << ", " << straw1.id()
               << endl;
          cout << "It should be                   : "
               << strawSpacing << " diff: "
               << (xLayerDeltaMag-strawSpacing) << endl;

          throw cet::exception("GEOM")  << "Incorrect interlayer straw spacing \n";

        }
      } // loop over straws in panel

      for (  auto is = panel.getStrawPointers().cbegin()+2;
             is < (panel.getStrawPointers().cend()); ++is) {
        const Straw& straw(**is);
        uint16_t sn = straw.id().getStraw();
        if ( sn%2 != 0 ) continue;
        auto i0 = is;
        auto i1 = is-1;
        if ( _innermostLayer == 1 ){
          i0 = is-1;
          i1 = is;
        }
        const Straw& straw0(**i0);
        const Straw& straw1(**i1);
        double xLayerDeltaMag =
          (straw0.getMidPoint() - straw1.getMidPoint()).mag();
        StrawId lid0 = straw0.id().getLayerId();
        StrawId lid1 = straw1.id().getLayerId();
        if (_verbosityLevel>2) {
          cout << __func__ << " Checking spacing"
               << " for layer " << lid0 << " straw " << straw0.id()
               << " and for layer " << lid1 << " straw " << straw1.id()  << endl;
        }
        if ( abs(xLayerDeltaMag-strawSpacing)> tolerance ) {
          cout << "xLayer straw spacing is (mm)   : "
               << xLayerDeltaMag
               << " for straws: "
               << straw0.id() << ", " << straw1.id()
               << endl;
          cout << "It should be                   : "
               << strawSpacing << " diff: "
               << (xLayerDeltaMag-strawSpacing) << endl;

          throw cet::exception("GEOM")  << "Incorrect interlayer straw spacing \n";

        }
      }
    }

    // Alignment of angle of full plane


    g4tt->_panelEB._EBKey       = TubsParams(_EBKeyInnerRadius,
                                    _EBKeyOuterRadius,
                                    _EBKeyHalfLength,
                                    0.,
                                    _EBKeyPhiRange);

    g4tt->_panelEB._EBKeyShield = TubsParams(_EBKeyInnerRadius,
                                    _EBKeyOuterRadius,
                                    _EBKeyShieldHalfLength,
                                    0.,
                                    _EBKeyPhiRange);

    g4tt->_panelEB._EBKeyMaterial         = _EBKeyMaterial;
    g4tt->_panelEB._EBKeyShieldMaterial   = _EBKeyShieldMaterial;
    g4tt->_panelEB._EBKeyPhiExtraRotation = _EBKeyPhiExtraRotation;

  }  // end makePanel

  void TrackerMaker::makeLayer ( const StrawId& layId, CLHEP::Hep3Vector const& planeorigin, StrawCollection& allStraws ){

    // straws per panel
    constexpr int spp = StrawId::_nstraws;

    int ilay = layId.getLayer();
    int ipnl = layId.getPanel();
    int iplane = layId.getPlane();

    // Is layer zero closest to the straw or farthest from it.
    int factor = ( _layerZPattern == 0 ) ? ilay : (StrawId::_nlayers - ilay - 1);

    //    cout << "Debugging TrackerMaker ilay: " << ilay << endl;

    // Start to populate the layer.
    // |z| of straw center, relative to the center of the plane.
    // Sign is taken care of elsewhere.

    // see similar calc in computePanelBoxParams
    // the above commented out calculation places the straws at the edge of the manifold in Z

    double zOffset = _supportHalfThickness + _manifoldZEdgeExcessSpace +
      _strawOuterRadius + factor*2.*_layerHalfSpacing;

    // Rotation that puts wire direction and wire mid-point into their
    // correct orientations.
    CLHEP::HepRotationZ RZ(panelRotation(ipnl,layId.getPlane()));

    // Unit vector in the wire (U) direction. This depends on Station, Plane and Panel  Hardware convention is U points from HV to Cal (Duke convention)
    // See doc 888 table 4 and figure 9 and doc 5703 figure 1 for details
    CLHEP::Hep3Vector Udir;
    if((ipnl +iplane + iplane/2)%2)
      Udir = RZ*CLHEP::Hep3Vector(0.,1.,0.); // upstream-facing panel
    else
      Udir = RZ*CLHEP::Hep3Vector(0.,-1.,0.); // downstream-facing panel

    // Straw number within the layer; does not reset to zero at each manifold.
    // we number the straws starting from the most inner one across the two layers in the panel/panel
    // it will be 0 for layer0 and 1 for layer1
    int listraw(ilay-2);

    // we increase the number by 2, not 1

    // Add all of the straws
    for ( int iman=0; iman<_manifoldsPerEnd; ++iman ){

      // Inner edge of the innermost wire connected to this manifold.
      // this is layer dependent, take care of it at xstraw below

      double xA = _rInnermostWire + 2.*_manifoldHalfLengths.at(0)*iman;

      // Add all of the straws connected to this manifold.
      for ( int istr=0; istr<_strawsPerManifold; ++istr ){
        listraw +=2;

        // layers with fewer straws would complicate StrawSD, constructTrackerv, TrackerMaker

        // Construct straw midpoint in its base position in the
        // coord system of the plane envelope.
        // we will shift the "second" layer from the manifold edge

        double xstraw = (ilay%2==_innermostLayer) ?
          xA + (2*istr)*_strawOuterRadius + istr*_strawGap :
          xA + (2*istr)*_strawOuterRadius + istr*_strawGap + 2.0*_layerHalfShift;

        CLHEP::Hep3Vector mid( xstraw, 0., zOffset*panelZSide(ipnl, iplane ) );

        mid += planeorigin;

        // Rotate straw midpoint to its actual location.
        CLHEP::Hep3Vector offset = RZ*mid;

        // compute the half-length; this is just the distance to the manifold
        double halflen = sqrt( _innerSupportRadius*_innerSupportRadius - xstraw*xstraw);

        ++_strawTrckrConstrCount;
        ++_strawPanelConstrCount;

        // _tt->_strawExists[index.asInt()] = plane.exists();

        StrawId lsid(iplane, ipnl, listraw);
        // in the new tracker model the straws are placed in the
        // panels, not layers, so we have to reshuffle them here
        // before we place them in allStraws

        // number of panels placed
        int npp = _strawTrckrConstrCount/spp;
        // current straw in the panel is listraw (same as StrawId::straw())
        // we use int in case we "overcount" and rely on at() to tell us that
        // counter used to *place* the straws in an order 0..95, not 0,2..93,95
        int strawCountReCounted = npp*spp+listraw;

        allStraws.at(strawCountReCounted) =
          Straw( lsid,
                 offset,
                 Udir,
                 halflen
                 );

//        allStraws_p.at(lsid.asUint16()) = &allStraws.at(strawCountReCounted);
        // straw pointers are always stored by StrawId order
//        panelStraws2_p.at(lsid.straw()) = &allStraws.at(strawCountReCounted);
//        strawExists2.at(lsid.asUint16()) = plane.exists();

        if (_verbosityLevel>3) {
          std::ostringstream nsid("",std::ios_base::ate); // to write at the end
          nsid << lsid;

          cout << __func__ << " strCnt, strCntRc, npp, iplane, ipnl, opsc, listraw, _sid, StrawId:"
               << setw(6) << _strawTrckrConstrCount
               << setw(6) << strawCountReCounted
               << setw(6) << npp
               << setw(3) << iplane
               << setw(2) << ipnl
               << setw(3) << _strawPanelConstrCount
               << setw(3) << listraw
               << setw(6) << lsid.asUint16()
               << setw(17) << std::bitset<16>(lsid.asUint16())
               << " "
               << setw(6) << std::showbase << std::hex << lsid.asUint16()
               << " " << std::dec << std::noshowbase << setw(7) << nsid.str()
               // << " " << std::showbase << std::hex << setw(10) << &allStraws.at(_strawTrckrConstrCount)
               // << std::dec << std::noshowbase
               << endl;
        }

        // layer._straws.push_back(&allStraws.at(_strawTrckrConstrCount));
      }
    }

//std::cout << "<-<-<- makeLayer\n";
  } // end TrackerMaker::makeLayer



  void TrackerMaker::computeLayerSpacingAndShift(){
    double rspacing = 2*_strawOuterRadius + _strawGap;
    _layerHalfSpacing = (StrawId::_nlayers<=1) ? 0.0 : 0.25*sqrt(3.0)*rspacing; // sqrt(3.0) = 2*sin(60)
    _layerHalfShift   = (StrawId::_nlayers<=1) ? 0.0 : 0.25*rspacing;
  }

  void TrackerMaker::computeManifoldEdgeExcessSpace(){

    // Computes space between first/last straw and edge of manifold

    _manifoldXEdgeExcessSpace = _manifoldHalfLengths.at(0) -
      _strawOuterRadius*_strawsPerManifold -
      _strawGap*(_strawsPerManifold-1)*0.5;

    _manifoldZEdgeExcessSpace = _manifoldHalfLengths.at(2) - _strawOuterRadius -
      (StrawId::_nlayers-1)*_layerHalfSpacing;

    //cout << "Debugging,  _manifoldXEdgeExcessSpace, _manifoldZEdgeExcessSpace: " <<
    //     _manifoldXEdgeExcessSpace << ", " << _manifoldZEdgeExcessSpace << endl;

    if ( _manifoldXEdgeExcessSpace < 0.0 || _manifoldZEdgeExcessSpace < 0.0){
      throw cet::exception("GEOM")
        << "Manifolds are too small to hold straws!\n";
    }

  }

  // Compute the spacing for the given plane.
  double TrackerMaker::choosePlaneSpacing( int ipln ) const {


    if ( _spacingPattern == 0 ) {
      return ipln * _planeSpacing;
    }

    else if ( _spacingPattern == 1 ) {
      int station = ipln/2;
      int k = ipln%2;
      if (k == 0 ) {
        return  station*_planeSpacing - _planeHalfSeparation;
      } else if (k == 1 ) {
        return  station*_planeSpacing + _planeHalfSeparation;
      }
    }

    throw cet::exception("GEOM")
      << "Unrecognized separation pattern in TrackerMaker. \n";

  }

  double TrackerMaker::findFirstPlaneZ0() const{

    if ( _spacingPattern == 0 ) {
      return _planeSpacing*double(StrawId::_nplanes-1)/2.0;
    }

    else if ( _spacingPattern == 1 ) {
      int nStations = StrawId::_nstations;
      return double(nStations-1)/2.0 * _planeSpacing;
    }
    else {
      throw cet::exception("GEOM")
        << "Unrecognized separation pattern in TrackerMaker. \n";
    }
    return 0.0;
  }

  // Identify the neighbour straws for all straws in the tracker
  void TrackerMaker::makeSupportStructure(){
    auto g4tt = _tt->g4Tracker();

    SupportStructure& sup  = g4tt->_supportStructure;

    // Positions for the next few objects in Mu2e coordinates.
    TubsParams endRingTubs( _endRingInnerRadius, _endRingOuterRadius, _endRingHalfLength);

    //    TubsParams midRingTubs ( _endRingInnerRadius, _endRingOuterRadius, _midRingHalfLength);
    sup._stiffRings.push_back(PlacedTubs ( "TrackerEndRingUpstream",   endRingTubs, CLHEP::Hep3Vector( _xCenter, 0., _zCenter-_endRingZOffset), _endRingMaterial ));

    if ( _hasDownRing ) {
      TubsParams downRingTubs( _downRingInnerRadius, _downRingOuterRadius, _downRingHalfLength);
      sup._stiffRings.push_back(PlacedTubs ( "TrackerEndRingDownstream",
                                             downRingTubs, CLHEP::Hep3Vector( _xCenter, 0., _zCenter + _downRingZOffset),
                                             _downRingMaterial ));
    }

    {
      if ( StrawId::_nplanes%2 !=0 ){
        throw cet::exception("GEOM")
          << "TrackerMaker::makeSupportStructure expected an even number of planes. Saw " << StrawId::_nplanes << " planes.\n";
      }

      // From upstream end of most upstream station to the
      // downstream end of the most downstream station.
      // Including all materials.
      double overallLength = (StrawId::_nstations-1)*_planeSpacing + 2.*_planeHalfSeparation + 2.* _innerRingHalfLength;
      if ( _ttVersion > 3 ) overallLength = (StrawId::_nstations-1)*_planeSpacing +
                              2.*_planeHalfSeparation + 4.*_innerRingHalfLength
                              + _planePadding + 2.*_panelPadding;

      // we make support beams here (they used to be called staves)

      // Staves touch the big ring at the upstream end; they end flush
      // with the downstream edge of the most downstream station.

      double z1 = -(_endRingZOffset-_endRingHalfLength);
      double z2 = overallLength/2.;
      if ( _hasDownRing ) z2 = _downRingZOffset - _downRingHalfLength;
      double zoff = (z1+z2)/2.;
      double zHalf = ( z2-z1)/2.;

      // hold the beamTubsParams in a map

      std::unordered_map<std::string,TubsParams> supportBeamParams;
      std::unordered_map<std::string,TubsParams> supportServiceParams;

      // the top beam is different

      // size_t ibeam(0);

      // top beam is no longer there

      std::ostringstream bos("TrackerSupportBeam_",std::ios_base::ate); // to write at the end
      bos << std::setfill('0') << std::setw(2);


      size_t nsServices =  _beam1_servicePhi0s.size();

      if ( nsServices !=  _beam1_servicePhiEnds.size()
           || nsServices != _beam1_serviceMaterials.size()
           || nsServices != _beam1_serviceOuterRadii.size()
           || nsServices != _beam1_serviceCovRelThickness.size()
           || nsServices != _beam1_serviceMaterialsCov.size() ) {
        throw cet::exception("GEOM")
          << __func__ << "the number of all the beam service paramters has to be the same"
          << endl;
      }

      const size_t nssbeams = _beam1_phiSpans.size() - 1;

      for (size_t ibeam = 1; ibeam!=3; ++ibeam) {

        double phi00 = (ibeam == 1)
          ? _beam1_phiRange[0] -_beam1_phiSpans[0]
          : 180.0 - _beam1_phiRange[0] + _beam1_phiSpans[0]; // effectively 180 or 0;

        for (size_t ssbeam = 0; ssbeam != nssbeams; ++ssbeam) {

          bos.str("TrackerSupportBeam_");
          bos << std::setw(1) << ibeam << ssbeam;

          if ( _verbosityLevel > 0 ) {
            cout << __func__ << " making " <<  bos.str() << endl;
          }

          double deltaPhi = _beam1_phiSpans[ssbeam+1] - _beam1_phiSpans[ssbeam];
          double phi0     = (ibeam == 1)
            ? phi00 + _beam1_phiSpans[ssbeam]
            : phi00 - _beam1_phiSpans[ssbeam] - deltaPhi;

          double outerRadius = (ssbeam != 1)
            ? _beam1_outerRadius
            : _beam1_midRadius1; // the service section is different

          supportBeamParams.insert(std::pair<std::string,
                                   TubsParams>(bos.str(),
                                               TubsParams(_beam1_innerRadius,
                                                          outerRadius,
                                                          zHalf,
                                                          phi0*CLHEP::degree,
                                                          deltaPhi*CLHEP::degree)));

          sup._beamBody.push_back(PlacedTubs( bos.str(),
                                              supportBeamParams.at(bos.str()),
                                              CLHEP::Hep3Vector(_xCenter, 0., _zCenter+zoff),
                                              _beam1_material) );

        }

        // a service envelope

        {

          size_t ssbeam = 1;

          bos.str("TrackerSupportServiceEnvelope_");
          bos << std::setw(1) << ibeam << ssbeam;

          if ( _verbosityLevel > 0 ) {
            cout << __func__ << " making " <<  bos.str() << endl;
          }

          double deltaPhi = _beam1_phiSpans[ssbeam+1]-_beam1_phiSpans[ssbeam];
          double phi0     =  (ibeam == 1)
            ? phi00 + _beam1_phiSpans[ssbeam]
            : phi00 - _beam1_phiSpans[ssbeam] - deltaPhi;

          supportBeamParams.insert(std::pair<std::string,
                                   TubsParams>(bos.str(),
                                               TubsParams(_beam1_midRadius1,
                                                          _beam1_midRadius2,
                                                          zHalf,
                                                          phi0*CLHEP::degree,
                                                          deltaPhi*CLHEP::degree)));

          sup._beamBody.push_back( PlacedTubs( bos.str(),
                                               supportBeamParams.at(bos.str()),
                                               CLHEP::Hep3Vector(_xCenter, 0., _zCenter+zoff),
                                               _envelopeMaterial) );

        }

        // services

        // subdivide the services into a number of groups of same lengts for a given phi span

        // adjust the lengts and place all tubes touching the last station
        // adjust the covered arc accordingly

        for (size_t sservice = 0; sservice!=nsServices; ++sservice) {

          // the span of this service section
          double deltaPhi0 = _beam1_servicePhiEnds[sservice] - _beam1_servicePhi0s[sservice];

          bos.str("TrackerSupportServiceSectionEnvelope_");
          bos << std::setw(1) << ibeam << sservice;

          if ( _verbosityLevel > 0 ) {
            cout << __func__ << " making " <<  bos.str() << endl;
          }

          std::string boses =  bos.str();

          bos.str("TrackerSupportService_");
          bos << std::setw(1) << ibeam << sservice << "_";

          std::string boss =  bos.str();

          if ( _verbosityLevel > 0 ) {
            cout << __func__ << " making " <<  boss << endl;
          }

          for ( int ssservice = 0; ssservice!=StrawId::_nstations; ++ssservice) {

            bos.str(boss);
            bos << std::setw(2) << ssservice;

            if ( _verbosityLevel > 0 ) {
              cout << __func__ << " making " <<  bos.str() << endl;
            }

            double sHLength = zHalf*(StrawId::_nstations-ssservice)/StrawId::_nstations;
            double sOffset  = zoff + zHalf - sHLength;

            if ( _verbosityLevel > 0 ) {
              cout << __func__ << " sHLength, sOffset "
                   << sHLength << ", " << sOffset << endl;
            }

            // the span of one "sub" service of this section
            double deltaPhi  = deltaPhi0/StrawId::_nstations;
            // the starting position of the "sub" service
            double phi0      = (ibeam == 1)
              ? phi00 + _beam1_servicePhi0s[sservice] + deltaPhi*ssservice
              : phi00 - _beam1_servicePhi0s[sservice] - deltaPhi*(1+ssservice);

            if ( _verbosityLevel > 0 ) {
              cout << __func__ << " deltaPhi0, phi0, deltaPhi "
                   << deltaPhi0 << ", " << phi0 << ", " << deltaPhi << endl;
            }

            // approximate the service by the main part an a top cover/envelope with different materials

            if ( _beam1_serviceCovRelThickness[sservice] > 1.
                 || _beam1_serviceCovRelThickness[sservice] < 0. ) {
              throw cet::exception("GEOM")
                << __func__ << " beam1_serviceCovRelThickness out of 0...1 range "
                << _beam1_serviceCovRelThickness[sservice]
                << endl;
            }

            double cRadius = _beam1_serviceOuterRadii[sservice] +
              _beam1_serviceCovRelThickness[sservice] *
              ( _beam1_midRadius1 - _beam1_serviceOuterRadii[sservice] );

            // an envelope for this service section

            if ( ssservice == 0) {

              double phi0      = (ibeam == 1)
                ? phi00 + _beam1_servicePhi0s[sservice]
                : phi00 - _beam1_servicePhi0s[sservice] - deltaPhi0;

              supportServiceParams.insert(std::pair<std::string,
                                          TubsParams>(boses,
                                                      TubsParams(_beam1_midRadius1,
                                                                 _beam1_serviceOuterRadii[sservice],
                                                                 sHLength,
                                                                 phi0*CLHEP::degree,
                                                                 deltaPhi0*CLHEP::degree)));

              sup._beamServices.push_back(PlacedTubs( boses,
                                                      supportServiceParams.at(boses),
                                                      CLHEP::Hep3Vector(_xCenter, 0., _zCenter+sOffset),
                                                       _envelopeMaterial));


            }

            supportServiceParams.insert(std::pair<std::string,
                                        TubsParams>(bos.str(),
                                                    TubsParams(_beam1_midRadius1,
                                                               cRadius,
                                                               sHLength,
                                                               phi0*CLHEP::degree,
                                                               deltaPhi*CLHEP::degree)));

            sup._beamServices.push_back(PlacedTubs( bos.str(),
                                                    supportServiceParams.at(bos.str()),
                                                    CLHEP::Hep3Vector(_xCenter, 0., _zCenter+sOffset),
                                                    _beam1_serviceMaterials[sservice]));

            if (_beam1_serviceCovRelThickness[sservice]>0.) {

              bos << std::setw(1) << "_c";

              supportServiceParams.insert(std::pair<std::string,
                                          TubsParams>(bos.str(),
                                                      TubsParams(cRadius,
                                                                  _beam1_serviceOuterRadii[sservice],
                                                                 sHLength,
                                                                 phi0*CLHEP::degree,
                                                                 deltaPhi*CLHEP::degree)));

              sup._beamServices.push_back(PlacedTubs( bos.str(),
                                                      supportServiceParams.at(bos.str()),
                                                      CLHEP::Hep3Vector(_xCenter, 0., _zCenter+sOffset),
                                                      _beam1_serviceMaterialsCov[sservice]));

            }

          }

        }

      }

    }

    // Special for three items added to sup for v5
    {
      sup._panelPhiRange = _panelPhi;
      sup._panelPhiRibs  = _dphiRibs;
      sup._ribHalfAngle = _ribHalfAngle;
    }

    // Positions from here onward are in the coordinates of the plane envelope.
    // For versions above 3, coordinates of the panel envelope
    { // Build and place the "center" plate
      TubsParams centerPlateTubs( _innerRingOuterRadius, _outerRingOuterRadius, _centerPlateHalfLength,0., _panelPhi);
      CLHEP::Hep3Vector centerSpot(0.,0.,0.);
      if ( _ttVersion > 3) centerSpot.set(0.,0.,_innerRingHalfLength - _centerPlateHalfLength );
      sup._centerPlate = PlacedTubs( "TrackerSupportCenterPlate", centerPlateTubs, centerSpot, _centerPlateMaterial);
    }

    { // This is the inner ring
      TubsParams innerRingTubs( _innerRingInnerRadius, _innerRingOuterRadius, _innerRingHalfLength,0., _panelPhi);
      sup._innerRing = PlacedTubs( "TrackerSupportInnerRing", innerRingTubs, CLHEP::Hep3Vector(0.,0.,0), _innerRingMaterial );
    }

    { // This is the channel info for the inner ring
      double outerRadius        = _innerRingInnerRadius + _channelDepth;
      double channelHalfLength  = (_layerHalfSpacing-_strawOuterRadius) + 2.*_strawOuterRadius;
      TubsParams innerChannelTubs( _innerRingInnerRadius, outerRadius, channelHalfLength, 0., _panelPhi );
      sup._innerChannelUpstream   = PlacedTubs( "TrackerSupportInnerChannelUpstream",   innerChannelTubs, CLHEP::Hep3Vector(0.,0.,-_channelZOffset+_panelZOffset), _envelopeMaterial );
      sup._innerChannelDownstream = PlacedTubs( "TrackerSupportInnerChannelDownstream", innerChannelTubs, CLHEP::Hep3Vector(0.,0., _channelZOffset-_panelZOffset), _envelopeMaterial );
    }

    { // This is the Outer ring
      double halfLength = (_innerRingHalfLength - _centerPlateHalfLength)/2.;
      double dz         = _centerPlateHalfLength + halfLength;
      if ( _ttVersion > 3 ) {
        halfLength *= 2.0;
        dz = -_centerPlateHalfLength;
      }
      TubsParams outerRingTubs( _outerRingInnerRadius, _outerRingOuterRadius, halfLength,0.0, _panelPhi);
      sup._outerRingUpstream   = PlacedTubs ( "TrackerSupportOuterRingUpstream",   outerRingTubs, CLHEP::Hep3Vector(0.,0.,-dz), _outerRingMaterial );
      sup._outerRingDownstream = PlacedTubs ( "TrackerSupportOuterRingDownstream", outerRingTubs, CLHEP::Hep3Vector(0.,0., dz), _outerRingMaterial );
    }

    { // Cover plate
      double dz = _innerRingHalfLength-_coverHalfLength;
      if ( _ttVersion > 3 ) dz = -dz;
      TubsParams coverTubs( _innerRingOuterRadius, _outerRingInnerRadius, _coverHalfLength,0., _panelPhi);
      sup._coverUpstream   = PlacedTubs( "TrackerSupportCoverUpstream",   coverTubs, CLHEP::Hep3Vector(0.,0.,-dz), _coverMaterial );
      sup._coverDownstream = PlacedTubs( "TrackerSupportCoverDownstream", coverTubs, CLHEP::Hep3Vector(0.,0., dz), _coverMaterial );
    }

    { // Gas volume
      double halfLength = (_innerRingHalfLength - _centerPlateHalfLength - 2.*_coverHalfLength)/2.;
      double dz         = _centerPlateHalfLength + halfLength;
      if ( _ttVersion > 3 ) {
        halfLength = _innerRingHalfLength - _centerPlateHalfLength - 2.0 * _coverHalfLength;
        dz = 2.0 * _coverHalfLength - _centerPlateHalfLength;
      }

      TubsParams gasTubs( _innerRingOuterRadius, _outerRingInnerRadius, halfLength, 0., _panelPhi );
      sup._gasUpstream   = PlacedTubs ( "TrackerSupportGasUpstream",  gasTubs, CLHEP::Hep3Vector(0.,0.,-dz), _electronicsSpaceMaterial );
      sup._gasDownstream = PlacedTubs ( "TrackerSupportGasDownstream", gasTubs, CLHEP::Hep3Vector(0.,0., dz), _electronicsSpaceMaterial );
    }

    // Positions for the next two are in the coordinates of the electronics space.
    { // G10 (now just "Electronics" )
      TubsParams g10Tubs( _innerRingOuterRadius, _outerRingInnerRadius, _electronicsG10HalfLength, 0., _panelPhi );
      sup._g10Upstream   = PlacedTubs ( "TrackerSupportElecG10Upstream",   g10Tubs, CLHEP::Hep3Vector(0.,0.,-_electronicsG10HalfLength), _electronicsG10Material );
      sup._g10Downstream = PlacedTubs ( "TrackerSupportElecG10Downstream", g10Tubs, CLHEP::Hep3Vector(0.,0.,-_electronicsG10HalfLength), _electronicsG10Material );
    }

    {
      TubsParams cuTubs( _innerRingOuterRadius, _outerRingInnerRadius, _electronicsCuHhalfLength, 0., _panelPhi);
      sup._cuUpstream   = PlacedTubs ( "TrackerSupportElecCuUpstream",   cuTubs, CLHEP::Hep3Vector(0.,0.,_electronicsCuHhalfLength), _electronicsCuMaterial);
      sup._cuDownstream = PlacedTubs ( "TrackerSupportElecCuDownstream", cuTubs, CLHEP::Hep3Vector(0.,0.,_electronicsCuHhalfLength), _electronicsCuMaterial);
    }

  } // end of makeSupportStructure()
  // *******************************


  // This needs to know the z positions of the support rings
  void TrackerMaker::makeThinSupportRings(){
    auto g4tt = _tt->g4Tracker();
    SupportStructure& sup  = g4tt->_supportStructure;

    TubsParams thinRingTubs ( _endRingInnerRadius, _outerRingOuterRadius, _midRingHalfLength,
                              _midRingPhi0, _midRingdPhi); // by default, half rings, on the bottom part, but user-configurable

    for ( size_t i=0; i< _midRingSlot.size(); ++i){
      std::ostringstream name;
      int station = _midRingSlot.at(i);
      int ipln1 = station*2+1;
      int ipln2 = ipln1+1;
      if ( ipln2 >= (int)_tt->nPlanes() ){
        throw cet::exception("GEOM")
          << "Requested a thin support after station: "
          << station
          << " This is between planes: "
          << ipln1 << " and "
          << ipln2 << "\n"
          << "But there are only "
          << _tt->nPlanes()
          << " planes\n";
      }
      name << "ThinSupportRing_" << i;

      // Center the support in the gap between two stations.
      double z = 0.5*( _tt->getPlane(ipln1).origin().z() +  _tt->getPlane(ipln2).origin().z());
      sup._stiffRings.push_back(PlacedTubs ( name.str(),  thinRingTubs, CLHEP::Hep3Vector( _xCenter, 0., _zCenter+z), _midRingMaterial ));
    }

  }


  // Envelope that holds one plane ("TrackerPlaneEnvelope")
  void TrackerMaker::computePlaneEnvelope(){
    auto g4tt = _tt->g4Tracker();

    if ( _supportModel == SupportModel::simple ){
      double halfThick = g4tt->_supportParams.halfThickness() +
        2.*_manifoldHalfLengths[2];
      g4tt->_planeEnvelopeParams = TubsParams( _envelopeInnerRadius,
                                              g4tt->_supportParams.outerRadius(),
                                              halfThick);
    } else if ( _supportModel == SupportModel::detailedv0 ){
      // This is new for version 5, its existence doesn't affect earlier
      // versions.
      g4tt->_panelEnvelopeParams = TubsParams( _envelopeInnerRadius,
                                               _outerRingOuterRadius,
                                               _innerRingHalfLength
                                              + _panelPadding,
                                              0., _panelPhi);
      if ( _ttVersion > 3 ) {
        g4tt->_planeEnvelopeParams = TubsParams( _envelopeInnerRadius,
                                                _outerRingOuterRadius,
                                                2.0 * (_innerRingHalfLength
                                                       + _panelPadding )
                                                + _planePadding );
      } else {
        g4tt->_planeEnvelopeParams = TubsParams( _envelopeInnerRadius,
                                                _outerRingOuterRadius,
                                                _innerRingHalfLength);
      } // end of if on version in assigning plane envelope params
    }else{
      throw cet::exception("GEOM")
        << "Unknown value of _supportModel in TrackerMaker::computePlaneEnvelopeParams "
        << _supportModel
        << "\n";
    }
  }

  // Envelope that holds the full Tracker ("TrackerMother")
  void TrackerMaker::computeTrackerEnvelope(){
    auto g4tt = _tt->g4Tracker();
    if ( _supportModel == SupportModel::simple ){

      // Envelope of a single plane.
      TubsParams planeEnvelope = _tt->g4Tracker()->getPlaneEnvelopeParams();

      // Full length from center to center of the first and last planes.
      double fullLength = _tt->planes().back().origin().z()-_tt->planes().front().origin().z();

      // Remember the thickness of the planes.
      double halfLength = fullLength/2. + planeEnvelope.zHalfLength();

      g4tt->_innerTrackerEnvelopeParams = TubsParams( planeEnvelope.innerRadius(),
                                                     planeEnvelope.outerRadius(),
                                                     halfLength);

    } else if ( _supportModel == SupportModel::detailedv0 ){

      double fullLength = _tt->planes().back().origin().z()-_tt->planes().front().origin().z();
      double halfLength = fullLength/2. + _tt->g4Tracker()->getPlaneEnvelopeParams().zHalfLength();

      TubsParams val( _envelopeInnerRadius, _outerRingOuterRadius, halfLength);
      g4tt->_innerTrackerEnvelopeParams = val;

    } else{

      throw cet::exception("GEOM")
        << "Unknown value of _supportModel in TrackerMaker::computeTrackerEnvelopeParams "
        << _supportModel
        << "\n";
    }

  } //end TrackerMaker::computeTrackerEnvelopeParams



  double
  TrackerMaker::panelRotation(int ipnl,int ipln) const {
    if ( _rotationPattern >= 5 ){
      int jplane = ipln%4;
      int jpln   = ipnl + jplane*StrawId::_npanels;
      double phi = _panelBaseRotations.at(jpln);
      return phi;
    }
    int jplane = ipln%2;
    int ista = (ipln/2)%2;
    int jpln = ipnl + jplane*StrawId::_npanels;
    double phi = _panelBaseRotations.at(jpln);
    if(ista==1)phi += _oddStationRotation;
    return phi;
  }

  double
  TrackerMaker::panelZSide(int ipanel, int iplane) const {

    // Even panels are upstream.
    if ( _panelZPattern == 0 ){
      return _panelZSide.at(ipanel);
    }
    if ( _panelZPattern != 1){
      // throw
    }

    // Mu2e-doc-888-v9, Figure 9;
    int jplane = iplane%4;
    double sign = ( jplane == 0 || jplane == 3) ? -1. : 1.;
    return sign*_panelZSide.at(ipanel);
  }

  // made inside the panel
  // void
  // TrackerMaker::makePanelEBKey(int ipanel, int iplane) {
  // }

} // namespace mu2e
