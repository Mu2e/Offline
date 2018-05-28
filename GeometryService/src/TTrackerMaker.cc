//
// Construct and return a TTracker.
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "Alignment/inc/AlignmentMap.hh"
#include "Alignment/inc/AlignmentObj.hh"
#include "Alignment/inc/AlignmentService.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "GeometryService/inc/TTrackerMaker.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "cetlib/pow.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

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

  // Constructor that gets information from the config file instead of
  // from arguments.
  TTrackerMaker::TTrackerMaker( SimpleConfig const& config){
    parseConfig(config);

    // Determine if Alignment is being used and set up
    myAlignMap = NULL;
    if( useAlignment) myAlignMap = art::ServiceHandle<AlignmentService>()->alignmentMap();

    buildIt( );

    // print straw layout for debbugging pourposes

    if (_verbosityLevel>2) {

      int ipln = -1;
      int ipnl = -1;
      int ilay = -1;
      double iang = -36000;

      cout << __func__ << " (_tt->_allStraws2).size(), TTracker::_nttstraws, _tt->_nStraws "
           << fixed << setw(6) << _tt->_allStraws2.size()
           << fixed << setw(6) << TTracker::_nttstraws
           << fixed << setw(6) << _tt->_nStraws
           << endl;

      size_t istr = -1;
      for (const auto& istr_p : _tt->_allStraws2_p) {
        cout << __func__ << setw(10) << ++istr
          // << setw(20) << istr_p;
             << setw(20) << " ";
        if (istr_p != nullptr ) {
          StrawId const & lsid =  (*istr_p).id();
          std::ostringstream nsid("",std::ios_base::ate); // to write at the end
          nsid << lsid;
          cout << setw(6) << lsid.asUint16()
               << setw(17) << std::bitset<16>(lsid.asUint16())
               << " "
               << setw(10) << std::showbase << std::hex << lsid.asUint16()
               << " " << std::dec << std::noshowbase << setw(7) << nsid.str()
               << " " << _tt->_strawExists2.at(lsid.asUint16())
               << endl;
        } else {
          cout << endl;
        }
      }

      for ( const Straw& straw : _tt->getAllStraws() ){

        size_t istr= straw.index().asInt();

        int cpln = straw.id().getPlane();
        int cpnl = straw.id().getPanel();
        int clay = straw.id().getLayer();

        const Plane& plane = _tt->getPlane(cpln);
        const Panel& panel = plane.getPanel(cpnl);

        size_t nStrawsPerPanel = panel.nStraws();
        size_t nStrawsPerPlane = plane.nPanels() * nStrawsPerPanel;

        double cang = panel.boxRzAngle()/M_PI*180.;
        double dang = plane.rotation()/M_PI*180.;

        size_t ipnlf = nStrawsPerPanel*cpnl + nStrawsPerPlane*cpln;

        cout << __func__ << " Straw "
             << fixed << setw(6) << istr << " "
             << " plnfloor " << setw(6) << ipnlf << " "
             << straw.id() << " "
             << " mid point " << straw.getMidPoint()
             << " r " << sqrt(straw.getMidPoint()[0]*straw.getMidPoint()[0]+
                              straw.getMidPoint()[1]*straw.getMidPoint()[1])
             << " direction " << straw.getDirection()
             << " panel rotation: " << cang
             << " straw0MidPoint  " << panel.straw0MidPoint()
             << " straw0Direction " << panel.straw0Direction()
             << " plane rotation: " << dang
             << " origin " << plane.origin()
             << " straw exists " << _tt->strawExists(straw.id())
             << " plane exists " << plane.exists();

        if (ipnl>cpnl && ipln==cpln) cout << " <--S";
        if (iang>cang && ipln==cpln) cout << " <--A";
        if (ilay>clay && ipnl==cpnl) cout << " <--L";
        if (ipln!=cpln) ipln=cpln;
        if (ipnl!=cpnl) ipnl=cpnl;
        if (ilay!=clay) ilay=clay;

        cout << endl;

      }

    }

  }

  void TTrackerMaker::parseConfig( const SimpleConfig& config ){

    useAlignment        = config.getBool("hasAlignment",false);
    _verbosityLevel     = config.getInt("ttracker.verbosityLevel",0);
    _ttVersion          = config.getInt("TTrackerVersion",3);

    _motherRIn        = config.getDouble("ttracker.mother.rIn"        )*CLHEP::mm;
    _motherROut       = config.getDouble("ttracker.mother.rOut"       )*CLHEP::mm;
    _motherHalfLength = config.getDouble("ttracker.mother.halfLength" )*CLHEP::mm;
    _motherZ0         = config.getDouble("ttracker.mother.z0"         )*CLHEP::mm;
    _numPlanes         = config.getInt("ttracker.numPlanes");
    if ( _numPlanes != StrawId::_nplanes ){
      if ( StrawId::_nplanes == 40 && _numPlanes == 36 ) {
        _verbosityLevel > -1 &&
          cout << __func__ << " Running with ttracker.numPlanes !=  StrawId::_nplanes "
               << _numPlanes << " != " << StrawId::_nplanes
               << " please make sure it is intended " << endl;
      } else {
        cout << __func__
             << " Inconsistent TTracker parameters: "
             << " ttracker.numPlanes !=  StrawId::_nplanes "
             << _numPlanes << " != " << StrawId::_nplanes
             << " please double check and act accordingly " << endl;
        throw cet::exception("GEOM") << "See above for the message from " << __func__;
      }
    }
    _panelsPerPlane    = config.getInt("ttracker.panelsPerPlane");
    if ( _panelsPerPlane != StrawId::_npanels ){
      cout << __func__
           << " Inconsistent TTracker parameters: "
           << " ttracker.numPlanes !=  StrawId::_npanels "
           << _panelsPerPlane << " != " << StrawId::_npanels
           << " please double check and act accordingly " << endl;
      throw cet::exception("GEOM") << "See above for the message from " << __func__;
    }
    _layersPerPanel    = config.getInt("ttracker.layersPerPanel");
    if ( _layersPerPanel != StrawId::_nlayers ){
      cout << __func__
           << " Inconsistent TTracker parameters: "
           << " ttracker.numPlanes !=  StrawId::_nlayers "
           << _layersPerPanel << " != " << StrawId::_nlayers
           << " please double check and act accordingly " << endl;
      throw cet::exception("GEOM") << "See above for the message from " << __func__;
    }
    _manifoldsPerEnd    = config.getInt("ttracker.manifoldsPerEnd");
    _strawsPerManifold  = config.getInt("ttracker.strawsPerManifold");
    if ( _manifoldsPerEnd*_layersPerPanel*_strawsPerManifold != StrawId::_nstraws ){
      cout << __func__
           << " Inconsistent TTracker parameters: "
           << " ttracker.strawsPerManifold*ttracker.manifoldsPerEnd*ttracker.layersPerPanel !=  StrawId::_nstraws "
           << _manifoldsPerEnd*_layersPerPanel*_strawsPerManifold
           << " != " << StrawId::_nstraws
           << " please double check and act accordingly " << endl;
      throw cet::exception("GEOM") << "See above for the message from " << __func__;
   }
    _rotationPattern    = config.getInt("ttracker.rotationPattern");
    _panelZPattern      = config.getInt("ttracker.panelZPattern");
    _layerZPattern      = config.getInt("ttracker.layerZPattern");
    _spacingPattern     = config.getInt("ttracker.spacingPattern");
    _innermostLayer     = config.getInt("ttracker.innermostLayer");

    _planePadding      = config.getDouble("ttracker.planePadding",0.5)*CLHEP::mm;
    _panelPadding      = config.getDouble("ttracker.panelPadding",0.25)*CLHEP::mm;

    _oddStationRotation   =  config.getDouble("ttracker.oddStationRotation")*CLHEP::degree;
    _zCenter              =  config.getDouble("ttracker.z0")*CLHEP::mm;
    _xCenter              = -config.getDouble("mu2e.solenoidOffset")*CLHEP::mm;
    _envelopeInnerRadius  =  config.getDouble("ttracker.envelopeInnerRadius")*CLHEP::mm;
    _rInnermostWire       =  config.getDouble("ttracker.rInnermostWire")*CLHEP::mm;
    _strawOuterRadius     =  config.getDouble("ttracker.strawOuterRadius")*CLHEP::mm;
    _strawWallThickness   =  config.getDouble("ttracker.strawWallThickness")*CLHEP::mm;
    _strawGap             =  config.getDouble("ttracker.strawGap")*CLHEP::mm;
    _planeSpacing        =  config.getDouble("ttracker.planeSpacing")*CLHEP::mm;
    _planeHalfSeparation =  config.getDouble("ttracker.planeHalfSeparation")*CLHEP::mm;

    _outerSupportRadius   =  config.getDouble("ttracker.outerSupportRadius")*CLHEP::mm;
    _innerSupportRadius   =  config.getDouble("ttracker.innerSupportRadius")*CLHEP::mm;
    _supportHalfThickness =  config.getDouble("ttracker.supportHalfThickness")*CLHEP::mm;
    _wireRadius           =  config.getDouble("ttracker.wireRadius")*CLHEP::mm;
    _manifoldYOffset      =  config.getDouble("ttracker.manifoldYOffset")*CLHEP::mm;
    _virtualDetectorHalfLength = config.getDouble("vd.halfLength")*CLHEP::mm;

    config.getVectorInt("ttracker.nonExistingPlanes", _nonExistingPlanes,  vector<int>() );

    _verbosityLevel > 0 && _nonExistingPlanes.size()>0 &&
      cout << __func__ << " inactive planes : f/l   "
           << _nonExistingPlanes.front() << " / "
           << _nonExistingPlanes.back()
           << endl;

    config.getVectorDouble("ttracker.manifoldHalfLengths", _manifoldHalfLengths, 3);
    for ( size_t i=0; i<_manifoldHalfLengths.size(); ++i ){
      _manifoldHalfLengths.at(i) *= CLHEP::mm;
    }

    config.getVectorString("ttracker.strawMaterials", _strawMaterials, 3);

    _envelopeMaterial = config.getString("ttracker.mat.vacuum");
    _supportMaterial = config.getString("ttracker.mat.support");

    _passivationMargin        =  config.getDouble("ttracker.passivationMargin")*CLHEP::mm;

    // For the detailv0 support model there are lots of extra parameters.
    _supportModel = SupportModel( config.getString("ttrackerSupport.model","simple"));
    if ( _supportModel == SupportModel::detailedv0 ) {
      _endRingOuterRadius      = config.getDouble( "ttrackerSupport.endRing.outerRadius" );
      _endRingInnerRadius      = config.getDouble( "ttrackerSupport.endRing.innerRadius" );
      _endRingHalfLength       = config.getDouble( "ttrackerSupport.endRing.halfLength"  );
      _endRingZOffset          = config.getDouble( "ttrackerSupport.endRing.zOffset"     );
      _endRingMaterial         = config.getString( "ttrackerSupport.endRing.material"    );

      _hasDownRing    = config.getBool( "ttrackerSupport.downRing.build",false);
      if ( _hasDownRing ) {
	_downRingOuterRadius      = config.getDouble( "ttrackerSupport.downRing.outerRadius" );
	_downRingInnerRadius      = config.getDouble( "ttrackerSupport.downRing.innerRadius" );
	_downRingHalfLength       = config.getDouble( "ttrackerSupport.downRing.halfLength"  );
	_downRingZOffset          = config.getDouble( "ttrackerSupport.downRing.zOffset"     );
	_downRingMaterial         = config.getString( "ttrackerSupport.downRing.material"    );
      }

      config.getVectorInt( "ttrackerSupport.midRing.slot", _midRingSlot );
      _midRingHalfLength       = config.getDouble(    "ttrackerSupport.midRing.halfLength" );
      _midRingPhi0 = config.getDouble( "ttrackerSupport.midRing.Phi0",180.0)*CLHEP::degree;
      _midRingdPhi = config.getDouble( "ttrackerSupport.midRing.dPhi",180.0)*CLHEP::degree;
      _midRingMaterial = config.getString( "ttrackerSupport.midRing.material", "StainlessSteel316");
      // support beams;
      // fixme use vectors to contain them all (e.g. vector<SupportBeamParams>)

      config.getVectorDouble( "ttrackerSupport.beam0.phiRange", _beam0_phiRange );
      _beam0_innerRadius     = config.getDouble( "ttrackerSupport.beam0.innerRadius" );
      _beam0_outerRadius     = config.getDouble( "ttrackerSupport.beam0.outerRadius" );
      _beam0_material        = config.getString( "ttrackerSupport.beam0.material" );

      config.getVectorDouble( "ttrackerSupport.beam1.phiRange", _beam1_phiRange );
      config.getVectorDouble( "ttrackerSupport.beam1.phiSpans", _beam1_phiSpans );
      config.getVectorDouble( "ttrackerSupport.beam1.servicePhi0s", _beam1_servicePhi0s );
      config.getVectorDouble( "ttrackerSupport.beam1.servicePhiEnds", _beam1_servicePhiEnds );
      _beam1_innerRadius = config.getDouble( "ttrackerSupport.beam1.innerRadius" );
      _beam1_midRadius1  = config.getDouble( "ttrackerSupport.beam1.midRadius1" );
      _beam1_midRadius2  = config.getDouble( "ttrackerSupport.beam1.midRadius2" );
      _beam1_outerRadius = config.getDouble( "ttrackerSupport.beam1.outerRadius" );
      _beam1_material    = config.getString( "ttrackerSupport.beam1.material" );
      config.getVectorDouble( "ttrackerSupport.beam1.serviceOuterRadii", _beam1_serviceOuterRadii );
      config.getVectorString( "ttrackerSupport.beam1.serviceMaterials", _beam1_serviceMaterials );
      config.getVectorDouble( "ttrackerSupport.beam1.serviceCovRelThickness", _beam1_serviceCovRelThickness );
      config.getVectorString( "ttrackerSupport.beam1.serviceMaterialsCov", _beam1_serviceMaterialsCov );


      _panelPhi  = config.getDouble("ttrackerSupport.phiCoverage",120.0)*CLHEP::degree;
      _dphiRibs = config.getDouble("ttrackerSupport.dphiRibs",27.0)*CLHEP::degree;
      _ribHalfAngle = config.getDouble("ttrackerSupport.ribHalfAngle",1.0)*CLHEP::degree;

      _innerRingInnerRadius    = config.getDouble( "ttrackerSupport.innerRing.innerRadius" );
      _innerRingOuterRadius    = config.getDouble( "ttrackerSupport.innerRing.outerRadius" );
      _innerRingHalfLength     = config.getDouble( "ttrackerSupport.innerRing.halfLength"  );
      _innerRingMaterial       = config.getString( "ttrackerSupport.innerRing.material"    );

      _centerPlateHalfLength   = config.getDouble( "ttrackerSupport.centerPlate.halfLength" );
      _centerPlateMaterial     = config.getString( "ttrackerSupport.centerPlate.material"   );

      _outerRingInnerRadius    = config.getDouble( "ttrackerSupport.outerRing.innerRadius" );
      _outerRingOuterRadius    = config.getDouble( "ttrackerSupport.outerRing.outerRadius" );
      _outerRingMaterial       = config.getString( "ttrackerSupport.outerRing.material"    );

      _coverHalfLength          = config.getDouble( "ttrackerSupport.cover.halfLength"           );
      _coverMaterial            = config.getString( "ttrackerSupport.cover.material"             );
      _electronicsG10HalfLength = config.getDouble( "ttrackerSupport.electronics.g10.halfLength" );
      _electronicsG10Material   = config.getString( "ttrackerSupport.electronics.g10.material"   );
      _electronicsCuHhalfLength = config.getDouble( "ttrackerSupport.electronics.cu.halfLength"  );
      _electronicsCuMaterial    = config.getString( "ttrackerSupport.electronics.cu.material"    );
      _channelZOffset           = config.getDouble( "ttrackerSupport.channel.zOffset"            );
      _panelZOffset             = config.getDouble( "ttrackerSupport.panel.zOffset", 0.0 );
      _channelDepth             = config.getDouble( "ttrackerSupport.channel.depth"              );
      _channelMaterial          = config.getString( "ttrackerSupport.channel.material"           );
      _electronicsSpaceMaterial = config.getString( "ttrackerSupport.electronicsSpace.material"  );

      _EBKeyHalfLength         = config.getDouble("ttrackerSupport.electronics.key.halfLength");
      _EBKeyShieldHalfLength   = config.getDouble("ttrackerSupport.electronics.key.shieldHalfLength");
      _EBKeyInnerRadius        = config.getDouble("ttrackerSupport.electronics.key.innerRadius");
      _EBKeyOuterRadius        = config.getDouble("ttrackerSupport.electronics.key.outerRadius");
      _EBKeyShiftFromPanelFace = config.getDouble("ttrackerSupport.electronics.key.shiftFromPanelFace");
      _EBKeyVisible            = config.getBool(  "ttrackerSupport.electronics.key.visible");
      _EBKeySolid              = config.getBool(  "ttrackerSupport.electronics.key.solid");
      _EBKeyShieldVisible      = config.getBool(  "ttrackerSupport.electronics.key.shieldVisible");
      _EBKeyShieldSolid        = config.getBool(  "ttrackerSupport.electronics.key.shieldSolid");
      _EBKeyMaterial           = config.getString("ttrackerSupport.electronics.key.material");
      _EBKeyShieldMaterial     = config.getString("ttrackerSupport.electronics.key.shieldMaterial");
      _EBKeyPhiRange           = config.getDouble("ttrackerSupport.electronics.key.phiRange")*CLHEP::degree;
      _EBKeyPhiExtraRotation   = config.getDouble("ttrackerSupport.electronics.key.phiExtraRotation")*CLHEP::degree;

      _wallOuterMetalThickness  = config.getDouble("ttracker.straw.wallOuterMetal.thickness")*CLHEP::mm;
      _wallInnerMetal1Thickness = config.getDouble("ttracker.straw.wallInnerMetal1.thickness")*CLHEP::mm;
      _wallInnerMetal2Thickness = config.getDouble("ttracker.straw.wallInnerMetal2.thickness")*CLHEP::mm;
      _wirePlateThickness       = config.getDouble("ttracker.straw.wirePlate.thickness")*CLHEP::mm;
      _wallOuterMetalMaterial   = config.getString("ttracker.straw.wallOuterMetal.material");
      _wallInnerMetal1Material  = config.getString("ttracker.straw.wallInnerMetal1.material");
      _wallInnerMetal2Material  = config.getString("ttracker.straw.wallInnerMetal2.material");
      _wirePlateMaterial        = config.getString("ttracker.straw.wirePlate.material");

    }

    if ( _numPlanes%2 != 0 ) {
      throw cet::exception("GEOM")  << "_numPlanes = " << _numPlanes
                                    << ": Current TTracker geometry assumes even number of planes  \n";
    }
    _numStations = _numPlanes/_planesPerStation;

    //string ttracker.mat.manifold  = "G4_Al";  // Placeholder.

    // Also define some parameters that may become variable some day.
    _panelBaseRotations.clear();
    _panelZSide.clear();

    if (_panelsPerPlane == 6 ){
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
      } else {
        throw cet::exception("GEOM")
          << "Unrecognized rotation pattern in TTrackerMaker. \n";
      }
    } else {
      throw cet::exception("GEOM")
        << "Unrecognized rotation pattern in TTrackerMaker. \n";
    }

    // Parts of the algorithm require that the Aseet style tracker is built
    // of stations with 2 planes each.
    if ( _spacingPattern == 1 ){
      if ( _numPlanes%2 == 1 ){
        throw cet::exception("GEOM")
          << "Aseet style tracker requires 2 planes per station.\n"
          << "So ttracker.numPlanes must be even.  It was: "
          << _numPlanes
          << "\n";
      }
    }

    _layerHalfSpacing = 0.0;
    _layerHalfShift = 0.0;
    _manifoldXEdgeExcessSpace = 0.0;
    _manifoldZEdgeExcessSpace = 0.0;

  } // end TTrackerMaker::parseConfig

  // void lptest( const Layer& lay){
  //   cout << lay.id() << " |  "
  //        << lay.nStraws()  <<  " "
  //        << endl;
  // }

  void plntest( const Plane& plane){
    cout << "Plane: "
         << plane.id() << " "
         << plane.origin() << " "
         << plane.rotation()
         << endl;
  }

  // void positionTest( const Layer& lay){
  //   const Straw& straw((lay.getStraw(0)));
  //   cout << "Layer: "
  //        << lay.id() << " "
  //        << straw.getMidPoint().z() <<  " "
  //        << straw.getDirection().z() << " "
  //        << endl;
  // }

  void TTrackerMaker::buildIt(){

    // as the array has constant size, we need a straw counter during construction
    _strawTrckrConstrCount = -1; // first straw will be at 0

    // Make an empty TTracker.
    _tt = unique_ptr<TTracker>(new TTracker());
    _tt->_nPlanes = _numPlanes;
    _tt->_nStraws = _numPlanes *
      StrawId::_npanels *
      StrawId::_nstraws;

    makeMother();

    computeLayerSpacingAndShift();
    computeManifoldEdgeExcessSpace();

    _tt->_supportModel = _supportModel;

    // Fill information about the new style, fully detailed support structure.
    if ( _supportModel == SupportModel::detailedv0 ) {
      makeSupportStructure();
    }

    // Fill the information about the, old style minimal supports.
    _tt->_supportParams = Support( _innerSupportRadius,
                                   _outerSupportRadius,
                                   _supportHalfThickness,
                                   _supportMaterial);

    _tt->_z0                  = _zCenter;
    _tt->_envelopeInnerRadius = _envelopeInnerRadius;
    _tt->_manifoldHalfLengths = _manifoldHalfLengths;
    _tt->_envelopeMaterial    = _envelopeMaterial;

    // Z location of the first plane.
    _z0 = -findFirstPlaneZ0();

    makeDetails();

    // Reserve space for straws so that pointers are valid.
    _nStrawsToReserve = _numPlanes * _panelsPerPlane * _layersPerPanel *
      _manifoldsPerEnd * _strawsPerManifold;

    if (_verbosityLevel>2) {
      cout << __func__ << " _nStrawsToReserve "
           << fixed << setw(6) << _nStrawsToReserve << endl;
    }

    // Construct the planes and their internals.
    for ( int ipln=0; ipln<_numPlanes; ++ipln ){
      makePlane( StrawId(ipln,0,0) );
    }

    // Fill all of the non-persistent information.
    _tt->fillPointers();

    identifyNeighbourStraws();

    // Stations
    // Construct the stations and their internals based on planes internals

    if (_verbosityLevel>2) {
      cout << __func__ << " _numStations: " << _numStations << endl;
    }
    _tt->_stations.reserve(_numStations);

    for ( int istation=0; istation<_numStations; ++istation ){
      makeStation( StationId(istation) );
    }

    // Order is important here.
    computePlaneEnvelope();
    computeTrackerEnvelope();
    recomputeHalfLengths();
    makeStrawTubs();

    // This uses information from the planes
    makeThinSupportRings();


    finalCheck();

    if ( _verbosityLevel > 0 ) {
      cout << "TTracker Support Structure: \n" << _tt->_supportStructure << endl;
    }

    // Test the forAll methods.
    //_tt->forAllLayers( lptest);
    //_tt->forAllPlanes( plntest);
    //_tt->forAllLayers( positionTest);

  } //end TTrackerMaker::buildIt.

  void TTrackerMaker::makeMother(){

    _tt->_mother = PlacedTubs ( "TrackerMother",
                                TubsParams( _motherRIn, _motherROut, _motherHalfLength),
                                CLHEP::Hep3Vector( _xCenter, 0., _motherZ0),
                                _envelopeMaterial );

  } //end TTrackerMaker::makeMother

  // In the present code the straw positions are computed using the manifold information.
  // The channel position is computed using the SupportStructure information.  These two
  // must be self consistent.  Relatively soon the manifold information will be removed
  // and they will be replaced with the support structure information.
  // For now, check for self-consistency.
  void TTrackerMaker::finalCheck( ){

    // Only do this test for the new model(s).

    StrawId s0(0,1,0);
    StrawId s1(0,1,1);
    double ztest = 0.5*(_tt->getStraw ( s0 ).getMidPoint().z()+_tt->getStraw ( s1 ).getMidPoint().z())
      -  _tt->getPlane( PlaneId(0) ).origin().z();
    ztest *= panelZSide(1,0);

    double tolerance = 1.e-6;
    if ( std::abs(ztest-_channelZOffset) > tolerance ){
      throw cet::exception("GEOM")  << "Inconsistent channel center and wire z location. \n"
                                    << "channel Z offset: " <<  _channelZOffset << "\n"
                                    << "plane z center:  " << _tt->getPlane( PlaneId(0) ).origin().z() << "\n "
                                    << "Straw z layer 0:  " << _tt->getStraw ( s0 ).getMidPoint().z() << "\n"
                                    << "Straw z layer 1:  " << _tt->getStraw ( s1 ).getMidPoint().z() << "\n"
                                    << "z test:           " << ztest << " delta " << std::abs(ztest-_channelZOffset) << "\n";
    }

  }


  void TTrackerMaker::makePlane( StrawId planeId ){

    //std::cout << "->->-> makePlane\n";
    int ipln = planeId.getPlane();

    double planeDeltaZ = choosePlaneSpacing(ipln);

    // Handle Alignment
    std::ostringstream tmpName;
    tmpName << "ttPlane" << ipln;
    std::string tmpNameS(tmpName.str());
    AlignmentObj tmpObj;
    if ( NULL != myAlignMap ) tmpObj = myAlignMap->find(tmpNameS);
    CLHEP::Hep3Vector alignTranslate(0,0,0);
    if ( tmpObj.isValid() ) {
      alignTranslate = tmpObj.displacement();
    }

    CLHEP::Hep3Vector origin( alignTranslate.x(), alignTranslate.y(), _z0+planeDeltaZ+alignTranslate.z());

    auto& planes = _tt->_planes;

    // plane rotation is no longer used.
    double phi = 0.0;
    planes.at(ipln) = Plane(planeId,origin, phi);
    Plane& plane = planes.at(ipln);
    if (_verbosityLevel>2) {
      cout << __func__ << " making plane " <<  plane.id();
    }
    plane._exists = ( find ( _nonExistingPlanes.begin(), _nonExistingPlanes.end(), ipln) ==
                      _nonExistingPlanes.end() );
    if (_verbosityLevel>2) {
      cout << ", exists " << plane._exists  << endl;
    }
    for ( int ipnl=0; ipnl<_panelsPerPlane; ++ipnl ){
      makePanel ( StrawId(ipln,ipnl,0), plane );
      if (_verbosityLevel>2) {
        size_t istr = -1; // local index in the panel
        Panel& panel = plane._panels.at(ipnl);
        std::ostringstream pnlid("",std::ios_base::ate); // to write at the end
        pnlid << panel.id();
        cout << __func__ << " straws in panel "
             << setw(7) << pnlid.str();
        // test of Plane, Panel functions
        cout << setw(4) << panel.id().getPlane();
        cout << setw(2) << panel.id().getPanel() << endl;

        for (const auto istr_p : panel._straws2_p) {
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
          // now the straw by the Panel::getStraw(const StrawId& strid2) which uses local index
          // panel.getStraw(StrawId(0,0,0)); // test
          nsid.str("");
          nsid << panel.getStraw(lsid).id();
          cout << setw(8) << nsid.str();
          nsid.str("");
          StrawId sid  = (*istr_p).id();
          nsid << panel.getStraw(sid).id(); // old id of the straw in the old container
          cout << setw(10) << nsid.str()
               << endl;
        }
      }
    }

//std::cout << "<-<-<- makePlane\n";
  }

  void TTrackerMaker::makePanel( const PanelId& pnlId, Plane& plane ){
//std::cout << "->->-> makePanel\n";

    plane._panels.at(pnlId.getPanel()) = Panel(pnlId);
    Panel& panel = plane._panels.at(pnlId.getPanel());

    _strawPanelConstrCount = -1;
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

    for ( int ilay=0; ilay<_layersPerPanel; ++ilay ){
      // we use the straw field to indicate the layer
      makeLayer(StrawId(pnlId.getPlane(),pnlId.getPanel(),ilay), panel);
    }

    for ( int ilay=0; ilay<_layersPerPanel; ++ilay ){

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
               << " " << straw.index()
               << endl;
        }

      }

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

        }
      }
    }

    // check spacing between layers/straws

    if (_layersPerPanel>1) {

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
      }

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

    panel._boxRzAngle = panelRotation( panel.id().getPanel(),plane.id().getPlane() ); //  is it really used? needed?

    // make EBkey

    // need to decide how to apply the rotation; it seems it is deferred to ConstructTTrackerTDR
    // let's do it the same way it is done for the panels

    // panel._EBKeys = PlacedTubs("EBKey",
    //                            TubsParams(_EBKeyInnerRadius,
    //                                       _EBKeyOuterRadius,
    //                                       _EBKeyHalfLength,
    //                                       0.,
    //                                       _EBKeyPhiRange),
    //                            EBKeyPosition,
    //                            EBKeyRotation,
    //                            _EBKeyMaterial);

    panel._EBKey       = TubsParams(_EBKeyInnerRadius,
                                    _EBKeyOuterRadius,
                                    _EBKeyHalfLength,
                                    0.,
                                    _EBKeyPhiRange);

    panel._EBKeyShield = TubsParams(_EBKeyInnerRadius,
                                    _EBKeyOuterRadius,
                                    _EBKeyShieldHalfLength,
                                    0.,
                                    _EBKeyPhiRange);

    panel._EBKeyMaterial         = _EBKeyMaterial;
    panel._EBKeyShieldMaterial   = _EBKeyShieldMaterial;
    panel._EBKeyPhiExtraRotation = _EBKeyPhiExtraRotation;

//std::cout << "<-<-<- makePanel\n";
  }  // makePanel

  void TTrackerMaker::makeLayer ( const StrawId& layId, Panel& panel ){
//std::cout << "->->-> makeLayer\n";

    // Make an empty layer object.
    // panel._layers.push_back( Layer(layId) );
    // Layer& layer = panel._layers.back(); // fixme: try to avoid this construction

    // Get additional bookkeeping info.

    // array type containers of straws and pointers, ttracker ones
    array<Straw,TTracker::_nttstraws>& allStraws2  = _tt->_allStraws2;
    array<Straw const*,TTracker::_maxRedirect>& allStraws2_p  = _tt->_allStraws2_p;
    // panel ones
    array<Straw const*, StrawId::_nstraws>& panelStraws2_p = panel._straws2_p;
    array<bool,TTracker::_maxRedirect>& strawExists2 = _tt->_strawExists2;
    // straws per panel
    constexpr int spp = StrawId::_nstraws;

    const Plane& plane = _tt->getPlane( layId );
    int ilay = layId.getLayer();
    int ipnl = layId.getPanel();

    int iplane = layId.getPlane();

    // Is layer zero closest to the straw or farthest from it.
    int factor = ( _layerZPattern == 0 ) ? ilay : (_layersPerPanel - ilay - 1);

    //    cout << "Debugging TTrackerMaker ilay: " << ilay << endl;

    // Start to populate the layer.
    // layer._nStraws      = _manifoldsPerEnd*_strawsPerManifold; // not really used
    // layer._straws.reserve(_manifoldsPerEnd*_strawsPerManifold);

    // |z| of straw center, relative to the center of the plane.
    // Sign is taken care of elsewhere.

    // see similar calc in computePanelBoxParams
    //    double zOffset = _supportHalfThickness + _strawOuterRadius + ilay*2.*_layerHalfSpacing;
    // the above commented out calculation places the straws at the edge of the manifold in Z

    double zOffset = _supportHalfThickness + _manifoldZEdgeExcessSpace +
      _strawOuterRadius + factor*2.*_layerHalfSpacing;

    // Rotation that puts wire direction and wire mid-point into their
    // correct orientations.
    // CLHEP::HepRotationZ RZ(_panelBaseRotations.at(ipnl));
    CLHEP::HepRotationZ RZ(panelRotation(ipnl,layId.getPlane()));

    // Unit vector in the wire direction. (nominal is the panel 0 to the right?)
    CLHEP::Hep3Vector unit = RZ*CLHEP::Hep3Vector(0.,1.,0.);

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

        // layers with fewer straws would complicate StrawSD, constructTTrackerv, TTrackerMaker

        // Construct straw midpoint in its base position in the
        // coord system of the plane envelope.
        // we will shift the "second" layer from the manifold edge

        double xstraw = (ilay%2==_innermostLayer) ?
          xA + (2*istr)*_strawOuterRadius + istr*_strawGap :
          xA + (2*istr)*_strawOuterRadius + istr*_strawGap + 2.0*_layerHalfShift;

        CLHEP::Hep3Vector mid( xstraw, 0., zOffset*panelZSide(ipnl, plane.id().getPlane()) );
        mid += plane.origin();

        // Rotate straw midpoint to its actual location.
        CLHEP::Hep3Vector offset = RZ*mid;

        ++_strawTrckrConstrCount;
        ++_strawPanelConstrCount;

        // _tt->_strawExists[index.asInt()] = plane.exists();

        StrawId lsid(iplane, ipnl, listraw);
        // in the new tracker model the straws are placed in the
        // panels, not layers, so we have to reshuffle them here
        // before we place them in allStraws2

        // number of panels placed
        int npp = _strawTrckrConstrCount/spp;
        // current straw in the panel is listraw (same as StrawId::straw())
        // we use int in case we "overcount" and rely on at() to tell us that
        // counter used to *place* the straws in an order 0..95, not 0,2..93,95
        int strawCountReCounted = npp*spp+listraw;

        allStraws2.at(strawCountReCounted) =
          Straw( lsid,
                 StrawIndex(strawCountReCounted),
                 offset,
                 &_tt->_strawDetails.at(iman*2+ilay%2),
                 iman*2+ilay%2,
                 unit
                 );

        allStraws2_p.at(lsid.asUint16()) = &allStraws2.at(strawCountReCounted);
        // straw pointers are always stored by StrawId order
        panelStraws2_p.at(lsid.straw()) = &allStraws2.at(strawCountReCounted);
        strawExists2.at(lsid.asUint16()) = plane.exists();

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
               // << " " << std::showbase << std::hex << setw(10) << &allStraws2.at(_strawTrckrConstrCount)
               // << std::dec << std::noshowbase
               << endl;
        }

        // layer._straws.push_back(&allStraws2.at(_strawTrckrConstrCount));
      }
    }

//std::cout << "<-<-<- makeLayer\n";
  } // end TTrackerMaker::makeLayer


// ======= Station view makers ============

  void TTrackerMaker::makeStation( StationId stationId ){

    int ist = stationId;
    int ipln1 = _planesPerStation*ist; // _planesPerStation is/has to be 2
    int ipln2 = ipln1 + 1;

    if (_verbosityLevel>2) {
      cout << __func__ << " StationId, plane1, plane2 :"
           << stationId << ", "
           << ipln1 << ", "
           << ipln2 << ", "
           << endl;
    }

    double stationZ = 0.5 *
        ( _tt->_planes.at(ipln1).origin().z() +
          _tt->_planes.at(ipln2).origin().z() );
    _tt->_stations.push_back(Station(stationId, stationZ));

    Station & st = _tt->_stations.back();
    st._planes.reserve (_planesPerStation);
    st._planes.push_back(_tt->_planes.at(ipln1));
    st._planes.push_back(_tt->_planes.at(ipln2));

  }

  // Assumes all planes and all panels are the same.
  // The straw length depends only on the manifold number.
  // See Mu2e-doc-??? for the algorithm.

  void TTrackerMaker::computeStrawHalfLengths(){

    // we use resize as we will set specific straw parameters using a priori known straw numbers

    _strawHalfLengths.clear();
    _strawHalfLengths.resize(_manifoldsPerEnd*2);
    _strawActiveHalfLengths.clear();
    _strawActiveHalfLengths.resize(_manifoldsPerEnd*2);

    for ( int i=0; i<_manifoldsPerEnd; ++i ){
      double xA = _rInnermostWire  + 2.*_manifoldHalfLengths.at(0)*i;

      //    double xA = (_layersPerPanel==1) ?
      //      _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*i + _manifoldXEdgeExcessSpace :
      //      _envelopeInnerRadius + 2.*_manifoldHalfLengths.at(0)*i + _manifoldXEdgeExcessSpace +
      //       _strawOuterRadius;

      // we ignore the further laying straws in the multi layer case,
      // as this would make the straws shorter than they need to be
      // the wire positioning is not affected by this though

      double yA = sqrt( diff_of_squares(_innerSupportRadius, xA) );
      double yB = yA + _manifoldYOffset;

      // this needs to be inserted at a specific position

      _strawHalfLengths[i*2]=yB;

      // This variable is only used if SupportModel==simple.
      double activeHLen=sqrt( diff_of_squares(_innerSupportRadius,xA+2.5))-_passivationMargin;
      activeHLen = std::max( activeHLen, 1.0);

      // this needs to be an insert at a specific position

      _strawActiveHalfLengths[i*2]=activeHLen;

    }

  } // end TTrackerMaker::computeStrawHalfLengths

  void TTrackerMaker::makeDetails(){

    computeStrawHalfLengths();

    if ( _supportModel == SupportModel::simple ){
      _tt->_strawDetails.reserve(_manifoldsPerEnd);
    } else{
      // This will be extended in recomputeHalfLengths - reserve enough space for the extention.
      // as we need to insert sparsly will resize now
      _tt->_strawDetails.resize(_manifoldsPerEnd*_layersPerPanel);
    }

    for ( int i=0; i<_manifoldsPerEnd; ++i ){
      // inserting at a specific place
      _tt->_strawDetails[i*2]=
        StrawDetail
        ( i*2,
          _strawMaterials,
          _strawOuterRadius,
          _strawWallThickness,
          _strawHalfLengths.at(i*2),
          _strawActiveHalfLengths.at(i*2),
          _wireRadius
          );
    }

  } // end TTrackerMaker::makeDetails

  void TTrackerMaker::computeLayerSpacingAndShift(){

    _layerHalfSpacing = (_layersPerPanel<=1) ? 0.0 :
      sqrt(3.0*(square(_strawOuterRadius)+_strawOuterRadius*_strawGap+0.25*square(_strawGap)))*0.5;
    _layerHalfShift   = (_layersPerPanel<=1) ? 0.0 : _strawGap*0.25 + _strawOuterRadius*0.5;

  }

  void TTrackerMaker::computeManifoldEdgeExcessSpace(){

    // Computes space between first/last straw and edge of manifold

    _manifoldXEdgeExcessSpace = _manifoldHalfLengths.at(0) -
      _strawOuterRadius*_strawsPerManifold -
      _strawGap*(_strawsPerManifold-1)*0.5;

    _manifoldZEdgeExcessSpace = _manifoldHalfLengths.at(2) - _strawOuterRadius -
      (_layersPerPanel-1)*_layerHalfSpacing;

    //cout << "Debugging,  _manifoldXEdgeExcessSpace, _manifoldZEdgeExcessSpace: " <<
    //     _manifoldXEdgeExcessSpace << ", " << _manifoldZEdgeExcessSpace << endl;

    if ( _manifoldXEdgeExcessSpace < 0.0 || _manifoldZEdgeExcessSpace < 0.0){
      throw cet::exception("GEOM")
        << "Manifolds are too small to hold straws!\n";
    }

  }

  // Compute the spacing for the given plane.
  double TTrackerMaker::choosePlaneSpacing( int ipln ) const {


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
      << "Unrecognized separation pattern in TTrackerMaker. \n";

  }

  double TTrackerMaker::findFirstPlaneZ0() const{

    if ( _spacingPattern == 0 ) {
      return _planeSpacing*double(_numPlanes-1)/2.0;
    }

    else if ( _spacingPattern == 1 ) {
      int nStations = _numPlanes/2;
      return double(nStations-1)/2.0 * _planeSpacing;
    }
    else {
      throw cet::exception("GEOM")
        << "Unrecognized separation pattern in TTrackerMaker. \n";
    }
    return 0.0;
  }


  // Identify the neighbour straws for all straws in the tracker
  void TTrackerMaker::identifyNeighbourStraws() {

    for (auto& straw : _tt->_allStraws2) {

      if (_verbosityLevel>2) {
        cout << __func__ << " "
             << straw.id() << ", index "
             << straw.index()
             << " Straw " << straw.id().getStraw()
             << " plane: "
             << _tt->getPlane(straw.id()).id()
             << endl;
      }

      // throw exception if more than 2 layers per panel
      if (_tt->getPlane(straw.id()).getPanel(straw.id()).nLayers() != 2 ) {
        throw cet::exception("GEOM")
          << "The code works with 2 layers per panel. \n";
      } // fixme: rewrite using panels

      //      LayerId lId = straw.id().getLayerId();
      int layer = straw.id().getLayer();
      // int nStrawLayer = _tt->getLayer(lId).nStraws();
      int nStrawLayer = StrawId::_nstraws/StrawId::_nlayers;

      if ( _verbosityLevel>2 ) {
        cout << __func__ << " layer " << layer
             << " of panel "  << straw.id().getPanelId()
             << " has " << nStrawLayer << " straws" << endl;
        cout << __func__ << " Analyzed straw: " << straw.id() << '\t' << straw.index() << endl;
      }

      // add the "same layer" n-2 neighbours straw (if exist)
      // in the new model straw numbers increase by 2 in a given layer

      if ( straw.id().getStraw() > 1 ) {
        // const StrawId nsId(straw.id().getPlane(), straw.id().getPanel(), straw.id().getStraw() - 2 );
        const StrawId nsId( straw.id().asUint16() - 2 );
        straw._nearestById.push_back( nsId );
        if ( _verbosityLevel>2 ) {
          const Straw& temp = _tt->getStraw( nsId );
          cout << __func__ << setw(34) << " Neighbour left straw: " << temp.id() << '\t' << temp.index() << endl;
        }
        straw._nearestByIndex.push_back( _tt->getStraw(nsId).index() );
      }

      // add the "same layer" n+2 neighbours straw (if exist)
      // in the new model straw numbers increase by 2 in a given layer

      if ( straw.id().getStraw() < (2*nStrawLayer-2) ) {
        const StrawId nsId( straw.id().asUint16() + 2 );
        straw._nearestById.push_back( nsId );
        if ( _verbosityLevel>2 ) {
          const Straw& temp = _tt->getStraw( nsId );
          cout << __func__ << setw(34) << " Neighbour right straw: " << temp.id() << '\t' << temp.index() << endl;
        }
        straw._nearestByIndex.push_back( _tt->getStraw(nsId).index() );
      }

      // add the "opposite layer" n neighbours straw (if more than 1 layer)

      if (_layersPerPanel == 2) {

        // throw exception if the two layer of the same panel have different
        // number of straws
        // fixme rewrite the check using panels only
        // if (_tt->getLayer(lId).nStraws() != nStrawLayer) { // fixme: always? the same lId for both?
        //   throw cet::exception("GEOM")
        //     << "The code works only with the same number of straws "
        //     << "per layer in the same panel. \n";
        // }

        // add all straws sharing a preamp
        // assumes current two channel preamps

        if ( straw.id().getStraw() % 2 == 0){
          const StrawId nsId( straw.id().asUint16() + 1 );
          straw._preampById.push_back( nsId );
          straw._preampByIndex.push_back( _tt->getStraw(nsId).index());
        }else{
          const StrawId nsId( straw.id().asUint16() - 1 );
          straw._preampById.push_back( nsId );
          straw._preampByIndex.push_back( _tt->getStraw(nsId).index());
        }

        // add neighbors

        if (layer==0 && straw.id().getStraw()<2*nStrawLayer) {
          const StrawId nsId( straw.id().asUint16() + 1 );
          straw._nearestById.push_back( nsId );
          straw._nearestByIndex.push_back( _tt->getStraw( nsId ).index() );
          if ( _verbosityLevel>2 ) {
            cout << __func__ << setw(34) << " Neighbour opposite up straw: "
                 << straw._nearestById.back() << '\t' <<  straw._nearestByIndex.back() << endl;
          }
        }

        if (layer==1 && straw.id().getStraw()<2*nStrawLayer-1) {
          const StrawId nsId( straw.id().asUint16() + 1 );
          straw._nearestById.push_back( nsId );
          straw._nearestByIndex.push_back( _tt->getStraw( nsId ).index() );
          if ( _verbosityLevel>2 ) {
            cout << __func__ << setw(34) << " Neighbour opposite up straw: "
                 << straw._nearestById.back() << '\t' <<  straw._nearestByIndex.back() << endl;
          }
        }

        if (layer==0 && straw.id().getStraw()>0) {
          const StrawId nsId( straw.id().asUint16() - 1 );
          straw._nearestById.push_back( nsId );
          straw._nearestByIndex.push_back( _tt->getStraw( nsId ).index() );
          if ( _verbosityLevel>2 ) {
            cout << __func__ << setw(34) << " Neighbour opposite down straw: "
                 << straw._nearestById.back() << '\t' <<  straw._nearestByIndex.back() << endl;
          }
        }

        if (layer==1 && straw.id().getStraw()>0) { // layer 1 straw 1 is ok
          const StrawId nsId( straw.id().asUint16() - 1 );
          straw._nearestById.push_back( nsId );
          straw._nearestByIndex.push_back( _tt->getStraw( nsId ).index() );
          if ( _verbosityLevel>2 ) {
            cout << __func__ << setw(34) << " Neighbour opposite down straw: "
                 << straw._nearestById.back() << '\t' <<  straw._nearestByIndex.back() << endl;
          }
        }

      }
    }

  } // identifyNeighborStraws

  void TTrackerMaker::makeSupportStructure(){

    SupportStructure& sup  = _tt->_supportStructure;
    _tt->_panelZOffset = _panelZOffset;

    // Positions for the next few objects in Mu2e coordinates.
    TubsParams endRingTubs( _endRingInnerRadius, _endRingOuterRadius, _endRingHalfLength);

    //    TubsParams midRingTubs ( _endRingInnerRadius, _endRingOuterRadius, _midRingHalfLength);
    sup._stiffRings.push_back(PlacedTubs ( "TTrackerEndRingUpstream",   endRingTubs, CLHEP::Hep3Vector( _xCenter, 0., _zCenter-_endRingZOffset), _endRingMaterial ));

    if ( _hasDownRing ) {
      TubsParams downRingTubs( _downRingInnerRadius, _downRingOuterRadius, _downRingHalfLength);
      sup._stiffRings.push_back(PlacedTubs ( "TTrackerEndRingDownstream",
					     downRingTubs, CLHEP::Hep3Vector( _xCenter, 0., _zCenter + _downRingZOffset),
					     _downRingMaterial ));
    }

    {
      if ( _numPlanes%2 !=0 ){
        throw cet::exception("GEOM")
          << "TTrackerMaker::makeSupportStructure expected an even number of planes. Saw " << _numPlanes << " planes.\n";
      }

      // From upstream end of most upstream station to the downstream end of the most downstream station.
      // Including all materials.
      double overallLength = (_numPlanes/2-1)*_planeSpacing + 2.*_planeHalfSeparation + 2.* _innerRingHalfLength;
      if ( _ttVersion > 3 ) overallLength = (_numPlanes/2-1)*_planeSpacing +
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

      std::ostringstream bos("TTrackerSupportBeam_",std::ios_base::ate); // to write at the end
      bos << std::setfill('0') << std::setw(2);

      // supportBeamParams.insert(std::pair<std::string,
      //                          TubsParams>(bos.str(),
      //                                      TubsParams(_beam0_innerRadius,
      //                                                 _beam0_outerRadius,
      //                                                 zHalf,
      //                                                 _beam0_phiRange[0]*CLHEP::degree,
      //                                                 (_beam0_phiRange[1]- _beam0_phiRange[0])*CLHEP::degree)));

      // sup._beamBody.push_back( PlacedTubs( bos.str(),
      //                                      supportBeamParams.at(bos.str()), // to make sure it exists
      //                                      CLHEP::Hep3Vector(_xCenter, 0., _zCenter+zoff),
      //                                      _beam0_material) );

      // make the first support beam (1) , the other one (2) is a mirror reflection
      // phi0 can be negative, deltaPhi must not

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

          bos.str("TTrackerSupportBeam_");
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

          bos.str("TTrackerSupportServiceEnvelope_");
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

          bos.str("TTrackerSupportServiceSectionEnvelope_");
          bos << std::setw(1) << ibeam << sservice;

          if ( _verbosityLevel > 0 ) {
            cout << __func__ << " making " <<  bos.str() << endl;
          }

          std::string boses =  bos.str();

          bos.str("TTrackerSupportService_");
          bos << std::setw(1) << ibeam << sservice << "_";

          std::string boss =  bos.str();

          if ( _verbosityLevel > 0 ) {
            cout << __func__ << " making " <<  boss << endl;
          }

          for ( int ssservice = 0; ssservice!=_numStations; ++ssservice) {

            bos.str(boss);
            bos << std::setw(2) << ssservice;

            if ( _verbosityLevel > 0 ) {
              cout << __func__ << " making " <<  bos.str() << endl;
            }

            double sHLength = zHalf*(_numStations-ssservice)/_numStations;
            double sOffset  = zoff + zHalf - sHLength;

            if ( _verbosityLevel > 0 ) {
              cout << __func__ << " sHLength, sOffset "
                   << sHLength << ", " << sOffset << endl;
            }

            // the span of one "sub" service of this section
            double deltaPhi  = deltaPhi0/_numStations;
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
      sup._centerPlate = PlacedTubs( "TTrackerSupportCenterPlate", centerPlateTubs, centerSpot, _centerPlateMaterial);
    }

    { // This is the inner ring
      TubsParams innerRingTubs( _innerRingInnerRadius, _innerRingOuterRadius, _innerRingHalfLength,0., _panelPhi);
      sup._innerRing = PlacedTubs( "TTrackerSupportInnerRing", innerRingTubs, CLHEP::Hep3Vector(0.,0.,0), _innerRingMaterial );
    }

    { // This is the channel info for the inner ring
      double outerRadius        = _innerRingInnerRadius + _channelDepth;
      double channelHalfLength  = (_layerHalfSpacing-_strawOuterRadius) + 2.*_strawOuterRadius;
      TubsParams innerChannelTubs( _innerRingInnerRadius, outerRadius, channelHalfLength, 0., _panelPhi );
      sup._innerChannelUpstream   = PlacedTubs( "TTrackerSupportInnerChannelUpstream",   innerChannelTubs, CLHEP::Hep3Vector(0.,0.,-_channelZOffset+_panelZOffset), _envelopeMaterial );
      sup._innerChannelDownstream = PlacedTubs( "TTrackerSupportInnerChannelDownstream", innerChannelTubs, CLHEP::Hep3Vector(0.,0., _channelZOffset-_panelZOffset), _envelopeMaterial );
    }

    { // This is the Outer ring
      double halfLength = (_innerRingHalfLength - _centerPlateHalfLength)/2.;
      double dz         = _centerPlateHalfLength + halfLength;
      if ( _ttVersion > 3 ) {
	halfLength *= 2.0;
	dz = -_centerPlateHalfLength;
      }
      TubsParams outerRingTubs( _outerRingInnerRadius, _outerRingOuterRadius, halfLength,0.0, _panelPhi);
      sup._outerRingUpstream   = PlacedTubs ( "TTrackerSupportOuterRingUpstream",   outerRingTubs, CLHEP::Hep3Vector(0.,0.,-dz), _outerRingMaterial );
      sup._outerRingDownstream = PlacedTubs ( "TTrackerSupportOuterRingDownstream", outerRingTubs, CLHEP::Hep3Vector(0.,0., dz), _outerRingMaterial );
    }

    { // Cover plate
      double dz = _innerRingHalfLength-_coverHalfLength;
      if ( _ttVersion > 3 ) dz = -dz;
      TubsParams coverTubs( _innerRingOuterRadius, _outerRingInnerRadius, _coverHalfLength,0., _panelPhi);
      sup._coverUpstream   = PlacedTubs( "TTrackerSupportCoverUpstream",   coverTubs, CLHEP::Hep3Vector(0.,0.,-dz), _coverMaterial );
      sup._coverDownstream = PlacedTubs( "TTrackerSupportCoverDownstream", coverTubs, CLHEP::Hep3Vector(0.,0., dz), _coverMaterial );
    }

    { // Gas volume
      double halfLength = (_innerRingHalfLength - _centerPlateHalfLength - 2.*_coverHalfLength)/2.;
      double dz         = _centerPlateHalfLength + halfLength;
      if ( _ttVersion > 3 ) {
	halfLength = _innerRingHalfLength - _centerPlateHalfLength - 2.0 * _coverHalfLength;
	dz = 2.0 * _coverHalfLength - _centerPlateHalfLength;
      }

      TubsParams gasTubs( _innerRingOuterRadius, _outerRingInnerRadius, halfLength, 0., _panelPhi );
      sup._gasUpstream   = PlacedTubs ( "TTrackerSupportGasUpstream",  gasTubs, CLHEP::Hep3Vector(0.,0.,-dz), _electronicsSpaceMaterial );
      sup._gasDownstream = PlacedTubs ( "TTrackerSupportGasDownstream", gasTubs, CLHEP::Hep3Vector(0.,0., dz), _electronicsSpaceMaterial );
    }

    // Positions for the next two are in the coordinates of the electronics space.
    { // G10 (now just "Electronics" )
      TubsParams g10Tubs( _innerRingOuterRadius, _outerRingInnerRadius, _electronicsG10HalfLength, 0., _panelPhi );
      sup._g10Upstream   = PlacedTubs ( "TTrackerSupportElecG10Upstream",   g10Tubs, CLHEP::Hep3Vector(0.,0.,-_electronicsG10HalfLength), _electronicsG10Material );
      sup._g10Downstream = PlacedTubs ( "TTrackerSupportElecG10Downstream", g10Tubs, CLHEP::Hep3Vector(0.,0.,-_electronicsG10HalfLength), _electronicsG10Material );
    }

    {
      TubsParams cuTubs( _innerRingOuterRadius, _outerRingInnerRadius, _electronicsCuHhalfLength, 0., _panelPhi);
      sup._cuUpstream   = PlacedTubs ( "TTrackerSupportElecCuUpstream",   cuTubs, CLHEP::Hep3Vector(0.,0.,_electronicsCuHhalfLength), _electronicsCuMaterial);
      sup._cuDownstream = PlacedTubs ( "TTrackerSupportElecCuDownstream", cuTubs, CLHEP::Hep3Vector(0.,0.,_electronicsCuHhalfLength), _electronicsCuMaterial);
    }

  } // end of makeSupportStructure()
  // *******************************


  // This needs to know the z positions of the support rings
  void TTrackerMaker::makeThinSupportRings(){
    SupportStructure& sup  = _tt->_supportStructure;

    TubsParams thinRingTubs ( _endRingInnerRadius, _outerRingOuterRadius, _midRingHalfLength,
                              _midRingPhi0, _midRingdPhi); // by default, half rings, on the bottom part, but user-configurable

    for ( size_t i=0; i< _midRingSlot.size(); ++i){
      std::ostringstream name;
      int station = _midRingSlot.at(i);
      int ipln1 = station*2+1;
      int ipln2 = ipln1+1;
      if ( ipln2 >= _tt->nPlanes() ){
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

  // Create all of the Tubs objects needed for G4 to describe a straw.
  void TTrackerMaker::makeStrawTubs(){
    std::vector<StrawDetail>& details = _tt->_strawDetails;

    for ( std::vector<StrawDetail>::iterator i=details.begin();
          i != details.end(); ++i ){

      StrawDetail& detail(*i);

      detail._outerMetalThickness  = _wallOuterMetalThickness;
      detail._innerMetal1Thickness = _wallInnerMetal1Thickness;
      detail._innerMetal2Thickness = _wallInnerMetal2Thickness;
      detail._wirePlateThickness   = _wirePlateThickness;
      detail._outerMetalMaterial   = _wallOuterMetalMaterial;
      detail._innerMetal1Material  = _wallInnerMetal1Material ;
      detail._innerMetal2Material  = _wallInnerMetal2Material;
      detail._wirePlateMaterial    = _wirePlateMaterial;

    }
  }

  // Envelope that holds one plane ("TTrackerPlaneEnvelope")
  void TTrackerMaker::computePlaneEnvelope(){

    if ( _supportModel == SupportModel::simple ){
      double halfThick = _tt->_supportParams.halfThickness() +
	2.*_tt->_manifoldHalfLengths[2];
      _tt->_planeEnvelopeParams = TubsParams( _tt->_envelopeInnerRadius,
					      _tt->_supportParams.outerRadius(),
					      halfThick);
    } else if ( _supportModel == SupportModel::detailedv0 ){
      // This is new for version 5, its existence doesn't affect earlier
      // versions.
      _tt->_panelEnvelopeParams = TubsParams( _envelopeInnerRadius,
                                               _outerRingOuterRadius,
                                               _innerRingHalfLength
					      + _panelPadding,
					      0., _panelPhi);
      if ( _ttVersion > 3 ) {
	_tt->_planeEnvelopeParams = TubsParams( _envelopeInnerRadius,
						_outerRingOuterRadius,
						2.0 * (_innerRingHalfLength
						       + _panelPadding )
						+ _planePadding );
      } else {
	_tt->_planeEnvelopeParams = TubsParams( _envelopeInnerRadius,
						_outerRingOuterRadius,
						_innerRingHalfLength);
      } // end of if on version in assigning plane envelope params
    }else{
      throw cet::exception("GEOM")
        << "Unknown value of _supportModel in TTrackerMaker::computePlaneEnvelopeParams "
        << _supportModel
        << "\n";
    }
  }

  // Envelope that holds the full TTracker ("TrackerMother")
  void TTrackerMaker::computeTrackerEnvelope(){

    if ( _supportModel == SupportModel::simple ){

      // Envelope of a single plane.
      TubsParams planeEnvelope = _tt->getPlaneEnvelopeParams();

      // Full length from center to center of the first and last planes.
      double fullLength = _tt->_planes.back().origin().z()-_tt->_planes.front().origin().z();

      // Remember the thickness of the planes.
      double halfLength = fullLength/2. + planeEnvelope.zHalfLength();

      _tt->_innerTrackerEnvelopeParams = TubsParams( planeEnvelope.innerRadius(),
                                                     planeEnvelope.outerRadius(),
                                                     halfLength);

    } else if ( _supportModel == SupportModel::detailedv0 ){

      double fullLength = _tt->_planes.back().origin().z()-_tt->_planes.front().origin().z();
      double halfLength = fullLength/2. + _tt->getPlaneEnvelopeParams().zHalfLength();

      TubsParams val( _envelopeInnerRadius, _outerRingOuterRadius, halfLength);
      _tt->_innerTrackerEnvelopeParams = val;

    } else{

      throw cet::exception("GEOM")
        << "Unknown value of _supportModel in TTrackerMaker::computeTrackerEnvelopeParams "
        << _supportModel
        << "\n";
    }

  } //end TTrackerMaker::computeTrackerEvnvelopeParams


  void TTrackerMaker::recomputeHalfLengths(){

    // This code is only valid for the model detailedv0.
    if ( _supportModel != SupportModel::detailedv0 ) {
      return;
    }

    size_t ipln0(0);
    size_t ipnl0(0);

    std::vector<StrawDetail>& allDetails = _tt->_strawDetails;
    //    size_t originalSize                  = allDetails.size();

    if (_verbosityLevel>2) {
      cout << __func__ << " Initial allDetails size "
           <<  allDetails.size()
           << endl;
    }

    // Step 1: Check that the pattern of _detailIndex is as expected.
    int nBad(0);
    for ( const auto& straw : _tt->_allStraws2 ){

      if (_verbosityLevel>2) {
        cout << __func__ << " checking straw "
             << straw._id.getStraw()
             << " id: "
             << straw.id()
             << " has detailIndex of: "
             << straw._detailIndex
             << endl;
      }

      if ( straw._detailIndex != straw._id.getStraw() ){
        ++nBad;
        cout << "Unexpected value of detailIndex. Straw "
             << straw.id()
             << " has detailIndex of: "
             << straw._detailIndex
             << endl;
      }
    }
    if ( nBad > 0 ){
      throw cet::exception("GEOM")
        << "TTRackerMaker::recomputeHalfLengths: patterm of _detailIndex is not as expected."
        << "\n";
    }

    // Inner and outer radii of the channel in inner ring.
    double rmin = _innerRingInnerRadius;
    double rmax = rmin + _channelDepth;

    // Number of straws that are too short.
    int nShort(0);

    // Step 2: For all layers in PanelId(0,0) recompute the straw lengths.
    //         For layers > 0:
    //            - create a new StrawDetail object to hold the new length.
    //            - Reseat the _detail and _detailIndex objects in the Straw object.
    //    vector<Layer>& lays = _tt->_planes.at(ipln0)._panels.at(ipnl0)._layers;

    Panel& panel = _tt->_planes.at(ipln0)._panels.at(ipnl0);

    for ( size_t ilay=0; ilay<_layersPerPanel; ++ilay){
      // Layer& lay(lays.at(ilay));

      // straws in layer are layed out contiguously, only their numbers increase by 2
      // this is going over straws in the layer
      uint16_t uist = -1;
      for (  auto ist = panel.getStrawPointers().cbegin();
             ist < panel.getStrawPointers().cend(); ++ist) {

        const Straw& straw(**ist);
        uint16_t sn = straw.id().getStraw();
        if ( sn%2 != ilay ) continue;
        int idx             = straw.index().asInt(); // straw index

        if (_verbosityLevel>2) {
          cout << __func__ << " recomputing: ist, idx "
               << ++uist << ", "
               << idx
               << " Straw " << straw._id.getStraw()
               << " id: "
               << straw._id
               << " detail index "
               << straw._detailIndex
               << endl;
        }

        StrawDetail& detail = ( ilay == 0 )
          ? allDetails.at(straw._detailIndex)
          : allDetails.at(straw._detailIndex-1);

        double r0 = straw.getMidPoint().perp();
        double r1 = r0 - _strawOuterRadius;
        double r2 = r0 + _strawOuterRadius;

        // Choose half length so that the outer edge of straw just touches the outer
        // limit of the channel.
        double hlen = sqrt(diff_of_squares( rmax, r2));

        // Active half-length of the straw; it may be shorter than the full length but not longer.
        double activeHalfLen = sqrt( diff_of_squares(rmin,r0) )-_passivationMargin;
        activeHalfLen = std::max( activeHalfLen, .0);
        activeHalfLen = std::min( activeHalfLen, hlen);

        // Check that the inner edge of the straw reaches the support
        double r3 = sqrt(sum_of_squares(hlen,r1));

        if ( r3 < rmin ){
          ++nShort;
          cout << "Straw is too short to reach the inner edge of the support.\n"
               << "Straw; " << straw.id()
               << " Radius at inner corner: " << r3 << " mm\n"
               << "Radius at inner edge of the support ring: " << rmin
               << endl;
        }

        if ( ilay == 0 ){
          detail._halfLength       = hlen;
          detail._activeHalfLength = activeHalfLen;
        } else {
          StrawDetail newDetail       = detail;
          newDetail._halfLength       = hlen;
          newDetail._activeHalfLength = activeHalfLen;

          // we need to use/set the detail index

          newDetail._id               = straw._detailIndex;

          allDetails[straw._id.getStraw()] = newDetail;
          straw._detail      = &allDetails.at(straw._id.getStraw());
          // straw._detailIndex = newDetail._id; // not needed in the new model

        }

        if (_verbosityLevel>2) {
          cout << __func__ << " after recomputing: ist, idx "
               << uist << ", "
               << idx
               << " Straw " << straw._id.getStraw()
               << " id: "
               << straw._id
               << endl;
        }

        const StrawDetail& theDetail = straw.getDetail();

        if ( _verbosityLevel > 2 ){

          cout << "Detail for: " << straw.id() << " " << theDetail.Id()            << endl;
          cout << "           outerTubsParams: " << theDetail.getOuterTubsParams()
               << theDetail.gasMaterialName()             << endl;
          cout << "           wallMother:      " << theDetail.wallMother()
               << theDetail.wallMotherMaterialName()      << endl;
          cout << "           wallOuterMetal:  " << theDetail.wallOuterMetal()
               << theDetail.wallOuterMetalMaterialName()  << endl;
          cout << "           wallCore         " << theDetail.wallCore()
               << theDetail.wallCoreMaterialName()        << endl;
          cout << "           wallInnerMetal1: " << theDetail.wallInnerMetal1()
               << theDetail.wallInnerMetal1MaterialName() << endl;
          cout << "           wallInnerMetal2: " << theDetail.wallInnerMetal2()
               << theDetail.wallInnerMetal2MaterialName() << endl;
          cout << "           wireMother:      " << theDetail.wireMother()
               << theDetail.wireMotherMaterialName()      << endl;
          cout << "           wirePlate:       " << theDetail.wirePlate()
               << theDetail.wirePlateMaterialName()       << endl;
          cout << "           wireCore:        " << theDetail.wireCore()
               << theDetail.wireCoreMaterialName()        << endl;
        }

      }
    }
    if ( nShort > 0 ){
      throw cet::exception("GEOM")
        << "TTRackerMaker::recomputeHalfLengths: some straws are too short.\n"
        << "Probably the answer is to deepen the channel."
        << "\n";
    }

    // Step 4: reseat _detail and _detailIndex for all other straws.
    for ( auto& straw : _tt->_allStraws2 ){

      // These are already done:
      if ( (straw.id().getPlane() == 0 ) &&
           straw.id().getPanel() == 0 ) continue;
      if ( straw.id().getLayer()   ==           0  ) continue;

      // Get the new detail object for this straw.
      int idx = straw.id().getStraw();

      if (_verbosityLevel>2) {
        cout << __func__ << " about to reset "
             << " Straw " << straw._id.getStraw()
             << " id: "
             << straw._id << " using detail: "
             << idx
             << endl;
      }

      StrawDetail const& detail = allDetails.at(idx);

      // Update the info about the detail object.
      straw._detail      = &detail;
      straw._detailIndex = idx;
    }

  } //end TTrackerMaker::recomputeHalfLengths
  double
  TTrackerMaker::panelRotation(int ipnl,int ipln) const {
    if ( _rotationPattern == 5 ){
      int jplane = ipln%4;
      int jpln   = ipnl + jplane*_panelsPerPlane;
      double phi = _panelBaseRotations.at(jpln);
      return phi;
    }
    int jplane = ipln%2;
    int ista = (ipln/2)%2;
    int jpln = ipnl + jplane*_panelsPerPlane;
    double phi = _panelBaseRotations.at(jpln);
    if(ista==1)phi += _oddStationRotation;
    return phi;
  }

  double
  TTrackerMaker::panelZSide(int ipanel, int iplane) const {

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
  // TTrackerMaker::makePanelEBKey(int ipanel, int iplane) {
  // }

} // namespace mu2e
