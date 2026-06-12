//
// Construct and return a Target.
//
// Original author Peter Shanahan
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>


// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GeometryService/inc/StoppingTargetMaker.hh"
#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/RotationY.h"
#include "CLHEP/Vector/RotationZ.h"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e {

  StoppingTargetMaker::StoppingTargetMaker(const CLHEP::Hep3Vector& detSysOrigin, SimpleConfig const& c)
    : _detSysOrigin(detSysOrigin)
    , _rIn(0.)
  {
    int verbosity(c.getInt("stoppingTarget.verbosity",0));

    // positions are in detector coordinates (mm).
    _z0InMu2e            = c.getDouble("stoppingTarget.z0InMu2e");
    _deltaZ        = c.getDouble("stoppingTarget.deltaZ");

    // all outer radii must be specified.
    c.getVectorDouble("stoppingTarget.radii", _rOut);
    _rIn = c.getDouble("stoppingTarget.holeRadius", 0);

    // Downstream code counts on this so test it here.
    if ( _rOut.size() < 1 ){
      throw cet::exception("GEOM")
        << "Specified a stopping target with no foils!\n";
    }

    // halfThicknesses can be repeated from last element specified
    c.getVectorDouble("stoppingTarget.halfThicknesses",_halfThicknesses);
    unsigned int size=_halfThicknesses.size();
    for (unsigned int ii=size; ii<_rOut.size(); ++ii)
      _halfThicknesses.push_back(_halfThicknesses[size-1]);

    // x variations can be repeated from last element specified
    c.getVectorDouble("stoppingTarget.xVars",_xVars);
    size=_xVars.size();
    for (unsigned int ii=size; ii<_rOut.size(); ++ii)
      _xVars.push_back(_xVars[size-1]);

    // y variations can be repeated from last element specified
    c.getVectorDouble("stoppingTarget.yVars",_yVars);
    size=_yVars.size();
    for (unsigned int ii=size; ii<_rOut.size(); ++ii)
      _yVars.push_back(_yVars[size-1]);

    // z variations can be repeated from last element specified
    c.getVectorDouble("stoppingTarget.zVars",_zVars);
    size=_zVars.size();
    for (unsigned int ii=size; ii<_rOut.size(); ++ii)
      _zVars.push_back(_zVars[size-1]);

    // x cosines can be repeated from last element specified
    c.getVectorDouble("stoppingTarget.xCos",_xCos);
    size=_xCos.size();
    for (unsigned int ii=size; ii<_rOut.size(); ++ii)
      _xCos.push_back(_xCos[size-1]);

    // y cosines can be repeated from last element specified
    c.getVectorDouble("stoppingTarget.yCos",_yCos);
    size=_yCos.size();
    for (unsigned int ii=size; ii<_rOut.size(); ++ii)
      _yCos.push_back(_yCos[size-1]);


    // Stopping target material determined from specified target material
    // in globalConstants_01.txt file
    if ( c.hasName("stoppingTarget.materials") ){
      throw cet::exception("GEOM")
        << "Specifying stopping target material in geometry file not allowed!\n"
        << "Material is specified in Mu2eG4/test/globalConstants_01.txt\n ";
    }

    _materials.assign( _rOut.size(), "StoppingTarget_"+GlobalConstantsHandle<PhysicsParams>()->getStoppingTargetMaterial() );

    // Search to see if override material is specifed in geom. file
    const string overrideMaterial = c.getString("stoppingTarget.overrideMaterial","");
    if ( !overrideMaterial.empty() ) {
      _materials.assign( _rOut.size(), overrideMaterial );
      cerr << endl
           << " WARNING: Stopping target material ( "
           << GlobalConstantsHandle<PhysicsParams>()->getStoppingTargetMaterial()
           << " ) overridden with: " << overrideMaterial << "\n"
           << endl
           << "     You are specifying a stopping target material\n"
           << "     that may not match what was chosen in your global\n"
           << "     constants file!  Do not do this unless you know\n"
           << "     what you are doing!\n"
           << endl;
    }

    // material of the target enclosing volume
    _fillMaterial=c.getString("stoppingTarget.fillMaterial");

    // stopping foil target supporting structure
    _foilTarget_supportStructure = c.getBool("stoppingTarget.foilTarget_supportStructure");
    _foilTarget_supportStructure_number = c.getInt("stoppingTarget.foilTarget_supportStructure_number"); // number of support wires per foil
    _foilTarget_supportStructure_angleOffset = c.getDouble("stoppingTarget.foilTarget_supportStructure_angleOffset"); // angle of first support wire wrt vertical
    _foilTarget_supportStructure_radius = c.getDouble("stoppingTarget.foilTarget_supportStructure_radius");// radius of the support wires
    _foilTarget_supportStructure_angleOffset = c.getDouble("stoppingTarget.foilTarget_supportStructure_angleOffset");// angle of first support wire wrt vertical
    _foilTarget_supportStructure_fillMaterial = c.getString("stoppingTarget.foilTarget_supportStructure_fillMaterial");

    if(c.getBool("stoppingTarget.foilTarget_supportStructure_endAtOPA", false)) { //if reaching OPA, get OPA parameters
      double opaR1 = c.getDouble("protonabsorber.outerPAInnerRadius0", 436.0);
      double opaR2 = c.getDouble("protonabsorber.outerPAInnerRadius1", 720.0);
      double opaZCenter = c.getDouble("protonabsorber.outerPAZCenter", 6250.0);
      double opaHL = c.getDouble("protonabsorber.outerPAHalfLength", 2250.0);
      unsigned nfoils = _rOut.size();
      double length = (nfoils-1)*_deltaZ; //length of stopping target
      int side = (opaR1 <= opaR2) ? 1 : -1; //to check if cone widens in Z or not, to find point of closest approach
      double zclosest = _z0InMu2e - side*length/2.;
      if(zclosest < opaZCenter - opaHL) zclosest = opaZCenter - opaHL; //if goes past the opa
      if(zclosest > opaZCenter + opaHL) zclosest = opaZCenter + opaHL; //if goes past the opa
      double rclosest = opaR1 + (opaR2-opaR1)*((zclosest - side*_foilTarget_supportStructure_radius)
                                               - (opaZCenter - opaHL))/(2.*opaHL); //linear radial increase from z = opa_z - opa_hl to opa_z + opa_hl (plus radius of the wire)
      _foilTarget_supportStructure_rOut = rclosest - 0.01; //add small buffer
    } else
       _foilTarget_supportStructure_rOut = 250.; //default value used previously of 500 mm diameter OPA

    // debugging print...
    if ( verbosity > 0 ) PrintConfig();

    // Do the real work.
    BuildIt( );
  }


  void StoppingTargetMaker::BuildIt(){

    // Build the Target Geometry.  This means MU2E internal geometry, not
    // Root, G4, or any other scheme.

    _targ = unique_ptr<StoppingTarget>(new StoppingTarget());

    // create the TargetFoils

    // compute the z-detector position of each foil, starting with nominal
    // position and including any specified variations.
    const int n0 = _rOut.size()/2;
    const double offset = ( _rOut.size()%2 == 1) ?
      _z0InMu2e : _z0InMu2e + _deltaZ/2.;

    for ( vector<double>::size_type i=0;
          i<_rOut.size(); ++i){

      // z position of the center of the foil.
      const double z = offset + (int(i)-n0)*_deltaZ + _zVars[i];

      // Only x and y directional cosines are specified in config.  Z
      // is derived locally

      double zCos=1-_xCos[i]*_xCos[i]-_yCos[i]*_yCos[i];
      if (zCos<0.) {
        throw cet::exception("RANGE") <<"Target Foil "<<i
                                      <<" ZCos^2 is negative, ="<< zCos<<"\n";
      } else {
        zCos=sqrt(zCos);
      }

      _targ->_foils.push_back( TargetFoil( i,
                                           CLHEP::Hep3Vector(_xVars[i] + _detSysOrigin.x(),
                                                             _yVars[i] + _detSysOrigin.y(),
                                                             z),
                                           CLHEP::Hep3Vector(_xCos[i],_yCos[i],zCos),
                                           _rOut[i],
                                           _rIn,
                                           _halfThicknesses[i],
                                           _materials[i],
                                           _detSysOrigin
                                           )
                      );
      // create the TargetFoilSupportStructure
      if (_foilTarget_supportStructure) {
              for ( int j=0; j<_foilTarget_supportStructure_number; ++j){


                      double supportStructure_xPosition, supportStructure_yPosition, supportStructure_zPosition; // variables to account deviation of supporting structure center from foil center

                      supportStructure_xPosition = _xVars[i] + _detSysOrigin.x(); // x position is identical with x position of the foil
                      supportStructure_yPosition = _yVars[i] + _detSysOrigin.y(); // y position is identical with y position of the foil
                      supportStructure_zPosition = z; // z position is identical with z position of the foil

                      CLHEP::Hep3Vector supportStructure_foilCenterInMu2e(supportStructure_xPosition, supportStructure_yPosition, supportStructure_zPosition);

                      // calculate support structure length
                      double temp_foilTarget_supportStructure_length=0;

                      if ( ((_foilTarget_supportStructure_rOut - _rOut[i])) > 0 ) {
                              temp_foilTarget_supportStructure_length = ((_foilTarget_supportStructure_rOut - _rOut[i]));
                              temp_foilTarget_supportStructure_length -= 1; // remove 1mm to avoid overlap problems with mother volume
                      } else {
                              temp_foilTarget_supportStructure_length=0;
                      }
                      // cout << "foil  " << i << "  support structure " << j << "  temp_foilTarget_supportStructure_length = " << temp_foilTarget_supportStructure_length << endl;

                      _targ->_supportStructures.push_back( TargetFoilSupportStructure( j, i,
                                              supportStructure_foilCenterInMu2e,
                                              CLHEP::Hep3Vector(_xCos[i],_yCos[i],zCos),
                                              _foilTarget_supportStructure_radius,
                                              temp_foilTarget_supportStructure_length,
                                              _foilTarget_supportStructure_angleOffset,
                                              _rOut[i], // outer radius of the foil the wire connects to
                                              _foilTarget_supportStructure_fillMaterial,
                                              _detSysOrigin
                                              )
                                      );
              }
      }
    }// foil i


    // calculate the parameters of the enclosing cylinder
    //find the radius - maximum of foil radius + offset from axis
    double radius=-1;
    for (unsigned int ifoil=0; ifoil<_targ->_foils.size(); ifoil++)
      {
        double rtest=_targ->_foils[ifoil].rOut() + _targ->_foils[ifoil].centerInDetectorSystem().perp();
        radius=std::max(radius,rtest);
      }
    // beef it up by a mm
    radius+=1;

    // give it to the Target
    _targ->_radius=radius;

    // set the length to accomodate generous tilts to the first and last foils
    double zmin=_targ->_foils[0].centerInMu2e().z()-5;
    double zmax=_targ->_foils[_targ->_foils.size()-1].centerInMu2e().z()+5;

    _targ->_zLen=zmax-zmin;
    _targ->_centerInMu2e = CLHEP::Hep3Vector(_detSysOrigin.x(), _detSysOrigin.y(), (zmax+zmin)/2.);

    _targ->_fillMaterial=_fillMaterial;

  }//::BuildIt

  void StoppingTargetMaker::PrintConfig()
  {
    // printout the StoppingTargetMaker's understanding of what it needs to build.
    //  for debugging...

    std::cout<<"\n StoppingTargetMaker Input Configuration -----------------"<<std::endl;
    std::cout<<"\n Target System:"<<std::endl;
    std::cout
      <<"Stopping target Z0 in Mu2e ="<<_z0InMu2e
      <<", nominal spacing="<<_deltaZ
      <<", enclosing material="<<_fillMaterial
      <<std::endl;

    std::cout<<"\n Foils:"<<std::endl;
    std::cout
      <<"Total Foils="<<_rOut.size()<<std::endl;
    std::cout <<"Outer Radii="<<std::endl;
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_rOut[itf];
    cout <<std::endl;
    std::cout <<"1/2 Thicknesses="<<std::endl;
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_halfThicknesses[itf];
    cout <<std::endl;
    std::cout <<"Center Variations from Nominal"<<std::endl;
    cout <<"x:";
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_yVars[itf];
    cout <<std::endl;
    cout <<"y:";
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_zVars[itf];
    cout <<std::endl;
    cout <<"z:";
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_xVars[itf];
    cout <<std::endl;
    std::cout <<"Normal Directional Cosine"<<std::endl;
    cout <<"x:";
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_xCos[itf];
    cout <<std::endl;
    cout <<"y:";
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_yCos[itf];
    cout <<std::endl;
    std::cout <<"Material="<<std::endl;
    for (unsigned int itf=0; itf<_rOut.size(); itf++)
      cout <<" "<<_materials[itf];
    cout <<std::endl;
  }

} // namespace mu2e
