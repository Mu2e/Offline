//
// Construct and return MBS
//
//
// Original author KLG
//
// Notes
// see mu2e-doc-1351 for naming conventions etc...
//
// Updated by dnbrow01 on 30/10/2015 using mu2e-doc-1351 v7.
// - this produces MBS "Version 2"

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <sstream>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// Mu2e includes
#include "Offline/GeometryService/inc/MBSMaker.hh"
#include "Offline/MBSGeom/inc/MBS.hh"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  MBSMaker::MBSMaker(SimpleConfig const & _config,
                     double solenoidOffset)
  {
    // if( ! _config.getBool("hasMBS",false) ) return;

    // create an empty MBS
    _mbs = unique_ptr<MBS>(new MBS());

    // access its object through a reference

    MBS & mbs = *_mbs.get();

    parseConfig(_config);

    // now create the specific components - Version 1 first, then Version 2
    mbs._version = _MBSVersion;

    if ( _MBSVersion == 1 ) {

      mbs._zMax = _MBSCZ+_BSTSZ+_BSTSHLength;
      mbs._rMax = _SPBSLOuterRadius;
      if (mbs._rMax<_SPBSCOuterRadius) { mbs._rMax = _SPBSCOuterRadius; }
      if (mbs._rMax<_SPBSROuterRadius) { mbs._rMax = _SPBSROuterRadius; }
      if (mbs._rMax<_SPBSSup1OuterRadius) { mbs._rMax = _SPBSSup1OuterRadius; }
      if (mbs._rMax<_SPBSSup2OuterRadius) { mbs._rMax = _SPBSSup2OuterRadius; }
      //mbs._rMin = *std::min_element(_CLV2InnerRadii.begin(), _CLV2InnerRadii.end());
      mbs._rMin = 0.0; //need to go to zero, since we now have an absorber in the endplug hole
      mbs._totLength = 2.0*_BSTSHLength;


      CLHEP::Hep3Vector _MBSMOffsetInMu2e  = CLHEP::Hep3Vector(-solenoidOffset,0.,_MBSCZ);
      mbs._originInMu2e = _MBSMOffsetInMu2e;

      std::vector<double> MBSMCornersZ, MBSMCornersInnRadii, MBSMCornersOutRadii;
      MBSMCornersZ       .push_back(_BSTSZ-_BSTSHLength);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOuterRadius);

      MBSMCornersZ       .push_back(_SPBSLZ-_SPBSLHLength);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOuterRadius);

      MBSMCornersZ       .push_back(_SPBSLZ-_SPBSLHLength);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(mbs._rMax);

      MBSMCornersZ       .push_back(_SPBSRZ+_SPBSRHLength);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(mbs._rMax);

      MBSMCornersZ       .push_back(_SPBSRZ+_SPBSRHLength);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOuterRadius);

      MBSMCornersZ       .push_back(_BSTSZ+_BSTSHLength);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOuterRadius);

      mbs._pMBSMParams = std::unique_ptr<Polycone>
        (new Polycone(MBSMCornersZ,
                      MBSMCornersInnRadii,
                      MBSMCornersOutRadii,
                      _MBSMOffsetInMu2e,
                      "DSVacuum"));

      CLHEP::Hep3Vector _BSTSOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_BSTSZ);

      // Build a polycone for steel pipe, as opposed to a Tube, so that
      // Versions 1 and 2 work better together
      std::vector<double> BSTSCornersZ, BSTSCornersInnRadii, BSTSCornersOutRadii;
      BSTSCornersZ.push_back(-_BSTSHLength);
      BSTSCornersInnRadii.push_back(_BSTSInnerRadius);
      BSTSCornersOutRadii.push_back(_BSTSOuterRadius);
      BSTSCornersZ.push_back(_BSTSHLength);
      BSTSCornersInnRadii.push_back(_BSTSInnerRadius);
      BSTSCornersOutRadii.push_back(_BSTSOuterRadius);

      mbs._pBSTSParams = std::unique_ptr<Polycone>
        (new Polycone(BSTSCornersZ,
                      BSTSCornersInnRadii,
                      BSTSCornersOutRadii,
                      _BSTSOffsetInMu2e,
                      _BSTSMaterialName));


      CLHEP::Hep3Vector _SPBSSup1OffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_SPBSSup1Z);

      mbs._pSPBSSup1Params = std::unique_ptr<Tube>
        (new Tube(_SPBSSup1MaterialName,
                  _SPBSSup1OffsetInMu2e,
                  _SPBSSup1InnerRadius,
                  _SPBSSup1OuterRadius,
                  _SPBSSup1HLength));

      CLHEP::Hep3Vector _SPBSSup2OffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_SPBSSup2Z);

      mbs._pSPBSSup2Params = std::unique_ptr<Tube>
        (new Tube(_SPBSSup2MaterialName,
                  _SPBSSup2OffsetInMu2e,
                  _SPBSSup2InnerRadius,
                  _SPBSSup2OuterRadius,
                  _SPBSSup2HLength));

      CLHEP::Hep3Vector _SPBSLOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_SPBSLZ);

      mbs._pSPBSLParams = std::unique_ptr<Tube>
        (new Tube(_SPBSLMaterialName,
                  _SPBSLOffsetInMu2e,
                  _SPBSLInnerRadius,
                  _SPBSLOuterRadius,
                  _SPBSLHLength));

      CLHEP::Hep3Vector _SPBSCOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_SPBSCZ);

      mbs._pSPBSCParams = std::unique_ptr<Tube>
        (new Tube(_SPBSCMaterialName,
                  _SPBSCOffsetInMu2e,
                  _SPBSCInnerRadius,
                  _SPBSCOuterRadius,
                  _SPBSCHLength,
                  _SPBSCminAngle,
                  _SPBSCmaxAngle));

      CLHEP::Hep3Vector _SPBSROffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_SPBSRZ);

      mbs._pSPBSRParams = std::unique_ptr<Tube>
        (new Tube(_SPBSRMaterialName,
                  _SPBSROffsetInMu2e,
                  _SPBSRInnerRadius,
                  _SPBSROuterRadius,
                  _SPBSRHLength));

      CLHEP::Hep3Vector _BSTCOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_BSTCZ);

      std::size_t nBSTCSurfs = _BSTCLengths.size();
      if (nBSTCSurfs != _BSTCInnerRadii.size() ||  nBSTCSurfs != _BSTCOuterRadii.size()) {
        throw cet::exception("GEOM")
          << "BSTC has different number of radii and lengths, check the geom config file. \n";
      }
      std::vector<double> BSTCCornersZ, BSTCCornersInnRadii, BSTCCornersOutRadii;
      double BSTCtotLength=0;
      for (std::vector<double>::iterator length_it=_BSTCLengths.begin(); length_it!=_BSTCLengths.end(); ++length_it) {
        BSTCtotLength+=*length_it;
      }
      double tmpBSTCPntZ=/*_BSTCZ*/-0.5*BSTCtotLength;
      for (std::size_t isurf=0; isurf<nBSTCSurfs; ++isurf) {
        BSTCCornersZ.push_back(tmpBSTCPntZ);
        BSTCCornersInnRadii.push_back(_BSTCInnerRadii.at(isurf));
        BSTCCornersOutRadii.push_back(_BSTCOuterRadii.at(isurf));
        tmpBSTCPntZ+=_BSTCLengths.at(isurf);
        BSTCCornersZ.push_back(tmpBSTCPntZ);
        BSTCCornersInnRadii.push_back(_BSTCInnerRadii.at(isurf));
        BSTCCornersOutRadii.push_back(_BSTCOuterRadii.at(isurf));
      }

      mbs._pBSTCParams = std::unique_ptr<Polycone>
        (new Polycone(BSTCCornersZ,
                      BSTCCornersInnRadii,
                      BSTCCornersOutRadii,
                      _BSTCOffsetInMu2e,
                      _BSTCMaterialName));

      CLHEP::Hep3Vector _BSBSOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_BSBSZ);

      int nBSBSSurfs = _BSBSLengths.size();
      if (nBSBSSurfs != (int) _BSBSInnerRadii.size() ||  nBSBSSurfs != (int) _BSBSOuterRadii.size()) {
        throw cet::exception("GEOM")
          << "BSBS has different number of radii and lengths, check the geom config file. \n";
      }
      std::vector<double> BSBSCornersZ, BSBSCornersInnRadii, BSBSCornersOutRadii;
      double BSBStotLength=0;
      for (std::vector<double>::iterator length_it=_BSBSLengths.begin(); length_it!=_BSBSLengths.end(); ++length_it) {
        BSBStotLength+=*length_it;
      }
      double tmpBSBSPntZ=/*_BSBSZ*/-0.5*BSBStotLength;
      for (int isurf=0; isurf<nBSBSSurfs; ++isurf) {
        BSBSCornersZ.push_back(tmpBSBSPntZ);
        BSBSCornersInnRadii.push_back(_BSBSInnerRadii.at(isurf));
        BSBSCornersOutRadii.push_back(_BSBSOuterRadii.at(isurf));
        tmpBSBSPntZ+=_BSBSLengths.at(isurf);
        BSBSCornersZ.push_back(tmpBSBSPntZ);
        BSBSCornersInnRadii.push_back(_BSBSInnerRadii.at(isurf));
        BSBSCornersOutRadii.push_back(_BSBSOuterRadii.at(isurf));
      }

      mbs._pBSBSParams = std::unique_ptr<Polycone>
        (new Polycone(BSBSCornersZ,
                      BSBSCornersInnRadii,
                      BSBSCornersOutRadii,
                      _BSBSOffsetInMu2e,
                      _BSBSMaterialName));


      CLHEP::Hep3Vector _CLV2OffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_CLV2Z);

      int nCLV2Surfs = _CLV2Lengths.size();
      if (nCLV2Surfs != (int) _CLV2InnerRadii.size() ||  nCLV2Surfs != (int) _CLV2OuterRadii.size()) {
        throw cet::exception("GEOM")
          << "CLV2 has different number of radii and lengths, check the geom config file. \n";
      }

      std::vector<double> CLV2CornersZ, CLV2CornersInnRadii, CLV2CornersOutRadii;
      double CLV2totLength=0;
      for (std::vector<double>::iterator length_it=_CLV2Lengths.begin(); length_it!=_CLV2Lengths.end(); ++length_it) {
        CLV2totLength+=*length_it;
      }

      double tmpCLV2PntZ=/*_CLV2Z*/-0.5*CLV2totLength;
      for (int isurf=0; isurf<nCLV2Surfs; ++isurf) {
        CLV2CornersZ.push_back(tmpCLV2PntZ);
        CLV2CornersInnRadii.push_back(_CLV2InnerRadii.at(isurf));
        CLV2CornersOutRadii.push_back(_CLV2OuterRadii.at(isurf));
        tmpCLV2PntZ+=_CLV2Lengths.at(isurf);
        CLV2CornersZ.push_back(tmpCLV2PntZ);
        CLV2CornersInnRadii.push_back(_CLV2InnerRadii.at(isurf));
        CLV2CornersOutRadii.push_back(_CLV2OuterRadii.at(isurf));
      }

      mbs._pCLV2Params = std::unique_ptr<Polycone>
        (new Polycone(CLV2CornersZ,
                      CLV2CornersInnRadii,
                      CLV2CornersOutRadii,
                      _CLV2OffsetInMu2e,
                      _CLV2MaterialName));

    } else if ( _MBSVersion >= 2 ) {

      mbs._zMax = _MBSCZ+_BSTSZ+_BSTSHLength;
      mbs._rMax = _SPBSCOuterRadius;
      //mbs._rMin = *std::min_element(_CLV2InnerRadii.begin(), _CLV2InnerRadii.end());
      mbs._rMin = 0.0; //need to go to zero, since we now have an absorber in the endplug hole
      mbs._totLength = 2.0*_BSTSHLength;


      CLHEP::Hep3Vector _MBSMOffsetInMu2e  = CLHEP::Hep3Vector(-solenoidOffset,0.,_MBSCZ);
      mbs._originInMu2e = _MBSMOffsetInMu2e;

      std::vector<double>::iterator lenIter = _BSTSZLengths.begin();
      std::vector<double> MBSMCornersZ, MBSMCornersInnRadii, MBSMCornersOutRadii;
      double theZ = -_BSTSHLength;
      MBSMCornersZ       .push_back(theZ);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOutRadii[0]);

      theZ += *lenIter;
      MBSMCornersZ       .push_back(theZ);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOutRadii[0]);

      MBSMCornersZ       .push_back(theZ);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(mbs._rMax);

      ++lenIter;
      theZ += *lenIter;
      MBSMCornersZ       .push_back(theZ);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(mbs._rMax);

      MBSMCornersZ       .push_back(theZ);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOutRadii[2]);


      MBSMCornersZ       .push_back(_BSTSHLength);
      MBSMCornersInnRadii.push_back(mbs._rMin);
      MBSMCornersOutRadii.push_back(_BSTSOutRadii[2]);

      mbs._pMBSMParams = std::unique_ptr<Polycone>
        (new Polycone(MBSMCornersZ,
                      MBSMCornersInnRadii,
                      MBSMCornersOutRadii,
                      _MBSMOffsetInMu2e,
                      "DSVacuum"));

      CLHEP::Hep3Vector _BSTSOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_BSTSZ);

      // ******************************
      // Stainless Steel pipe
      // ******************************

      // In the geometry file, we have vectors with radii and lengths.
      // But a polycone requires planes at which the radii are defined,
      // with planes separated by lengths.  So first we have to translate the
      // lengths to plane positions and copy the radii to both ends.

      std::size_t nBSTSSurfs = _BSTSZLengths.size();
      if (nBSTSSurfs != _BSTSInnRadii.size() ||  nBSTSSurfs != _BSTSOutRadii.size()) {
        throw cet::exception("GEOM")
          << "BSTS has different number of radii and lengths, check the geom config file. \n";
      }
      std::vector<double> BSTSPlanes, BSTSInRads, BSTSOutRads;

      double zpos = 0.0;
      for ( unsigned int iSeg = 0; iSeg < _BSTSZLengths.size(); iSeg++ ) {

        // Fixme: ensure that zplanes in the BSTS polycone are not coincident.
        // Be sure to stay inside the MBS mother.
        const double epsilon = (iSeg>0 )? 0.1 : 0.0;

        BSTSPlanes.push_back(zpos - _BSTSHLength + _BSTSZ + epsilon);
        BSTSInRads.push_back(_BSTSInnRadii[iSeg]);
        BSTSOutRads.push_back(_BSTSOutRadii[iSeg]);
        zpos += _BSTSZLengths[iSeg];

        BSTSPlanes.push_back(zpos - _BSTSHLength + _BSTSZ);
        BSTSInRads.push_back(_BSTSInnRadii[iSeg]);
        BSTSOutRads.push_back(_BSTSOutRadii[iSeg]);
      }

      mbs._pBSTSParams = std::unique_ptr<Polycone>
        ( new Polycone(BSTSPlanes,
                       BSTSInRads,
                       BSTSOutRads,
                       _BSTSOffsetInMu2e,
                       _BSTSMaterialName ));

      // *****************************
      //     Holes in the Steel and poly
      //    (Added in version 4)
      // *****************************

      if ( _MBSVersion > 3 ) {
        mbs._holeXDimInSteel     = _BSTSHoleXDim;
        mbs._holeYDimInSteel     = _BSTSHoleYDim;
        mbs._holeZDimInSteel     = _BSTSHoleZDim;
        mbs._holeXDimInUpPoly    = _upPolyHoleXDim;
        mbs._holeYDimInUpPoly    = _upPolyHoleYDim;
        mbs._holeZDimInUpPoly    = _upPolyHoleZDim;
        mbs._holeXDimInDownPoly  = _downPolyHoleXDim;
        mbs._holeYDimInDownPoly  = _downPolyHoleYDim;
        mbs._holeZDimInDownPoly  = _downPolyHoleZDim;
        mbs._holeCentersInSteel  = _BSTSHoleCenters;
        mbs._holeCentersInUpstreamPoly = _upPolyHoleCenters;
        mbs._holeCentersInDownstreamPoly = _downPolyHoleCenters;
      }

      // ******************************
      // Outer HDPE
      // ******************************

      CLHEP::Hep3Vector _SPBSCOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_SPBSCZ);

      mbs._pSPBSCParams = std::unique_ptr<Tube>
        (new Tube(_SPBSCMaterialName,
                  _SPBSCOffsetInMu2e,
                  _SPBSCInnerRadius,
                  _SPBSCOuterRadius,
                  _SPBSCHLength,
                  _SPBSCminAngle,
                  _SPBSCmaxAngle));



      // ******************************
      // Inner HDPE upstream
      // ******************************


      CLHEP::Hep3Vector _BSTCOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_BSTCZ);

      std::size_t nBSTCSurfs = _BSTCLengths.size();
      if (nBSTCSurfs != _BSTCInnerRadii.size() ||  nBSTCSurfs != _BSTCOuterRadii.size()) {
        throw cet::exception("GEOM")
          << "BSTC has different number of radii and lengths, check the geom config file. \n";
      }
      std::vector<double> BSTCCornersZ, BSTCCornersInnRadii, BSTCCornersOutRadii;
      double BSTCtotLength=0;
      for (std::vector<double>::iterator length_it=_BSTCLengths.begin(); length_it!=_BSTCLengths.end(); ++length_it) {
        BSTCtotLength+=*length_it;
      }
      double tmpBSTCPntZ=/*_BSTCZ*/-0.5*BSTCtotLength;
      for (std::size_t isurf=0; isurf<nBSTCSurfs; ++isurf) {

        // Fixme: ensure that zplanes in the BSTC polycone are not coincident.
        // No overlap issues here.
        const double epsilon = (isurf>0 )? 0.1 : 0.0;

        BSTCCornersZ.push_back(tmpBSTCPntZ+epsilon);
        BSTCCornersInnRadii.push_back(_BSTCInnerRadii.at(isurf));
        BSTCCornersOutRadii.push_back(_BSTCOuterRadii.at(isurf));
        tmpBSTCPntZ+=_BSTCLengths.at(isurf);

        BSTCCornersZ.push_back(tmpBSTCPntZ);
        BSTCCornersInnRadii.push_back(_BSTCInnerRadii.at(isurf));
        BSTCCornersOutRadii.push_back(_BSTCOuterRadii.at(isurf));
      }

      mbs._pBSTCParams = std::unique_ptr<Polycone>
        (new Polycone(BSTCCornersZ,
                      BSTCCornersInnRadii,
                      BSTCCornersOutRadii,
                      _BSTCOffsetInMu2e,
                      _BSTCMaterialName));

      // ******************************
      // Inner HDPE downstream
      // ******************************


      CLHEP::Hep3Vector _BSBSOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_BSBSZ);

      int nBSBSSurfs = _BSBSLengths.size();
      if (nBSBSSurfs != (int) _BSBSInnerRadii.size() ||  nBSBSSurfs != (int) _BSBSOuterRadii.size()) {
        throw cet::exception("GEOM")
          << "BSBS has different number of radii and lengths, check the geom config file. \n";
      }
      std::vector<double> BSBSCornersZ, BSBSCornersInnRadii, BSBSCornersOutRadii;
      double BSBStotLength=0;
      for (std::vector<double>::iterator length_it=_BSBSLengths.begin(); length_it!=_BSBSLengths.end(); ++length_it) {
        BSBStotLength+=*length_it;
      }
      double tmpBSBSPntZ=/*_BSBSZ*/-0.5*BSBStotLength;
      for (int isurf=0; isurf<nBSBSSurfs; ++isurf) {
        BSBSCornersZ.push_back(tmpBSBSPntZ);
        BSBSCornersInnRadii.push_back(_BSBSInnerRadii.at(isurf));
        BSBSCornersOutRadii.push_back(_BSBSOuterRadii.at(isurf));
        tmpBSBSPntZ+=_BSBSLengths.at(isurf);
        BSBSCornersZ.push_back(tmpBSBSPntZ);
        BSBSCornersInnRadii.push_back(_BSBSInnerRadii.at(isurf));
        BSBSCornersOutRadii.push_back(_BSBSOuterRadii.at(isurf));
      }

      mbs._pBSBSParams = std::unique_ptr<Polycone>
        (new Polycone(BSBSCornersZ,
                      BSBSCornersInnRadii,
                      BSBSCornersOutRadii,
                      _BSBSOffsetInMu2e,
                      _BSBSMaterialName));


      // ******************************
      // HDPE end plug downstream end
      // ******************************

      CLHEP::Hep3Vector _CLV2OffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.,_CLV2Z);

      int nCLV2Surfs = _CLV2Lengths.size();
      if (nCLV2Surfs != (int) _CLV2InnerRadii.size() ||  nCLV2Surfs != (int) _CLV2OuterRadii.size()) {
        throw cet::exception("GEOM")
          << "CLV2 has different number of radii and lengths, check the geom config file. \n";
      }
      std::vector<double> CLV2CornersZ, CLV2CornersInnRadii, CLV2CornersOutRadii;
      double CLV2totLength=0;
      for (std::vector<double>::iterator length_it=_CLV2Lengths.begin(); length_it!=_CLV2Lengths.end(); ++length_it) {
        CLV2totLength+=*length_it;
      }
      double tmpCLV2PntZ=/*_CLV2Z*/-0.5*CLV2totLength;
      for (int isurf=0; isurf<nCLV2Surfs; ++isurf) {
        CLV2CornersZ.push_back(tmpCLV2PntZ);
        CLV2CornersInnRadii.push_back(_CLV2InnerRadii.at(isurf));
        CLV2CornersOutRadii.push_back(_CLV2OuterRadii.at(isurf));
        tmpCLV2PntZ+=_CLV2Lengths.at(isurf);
        CLV2CornersZ.push_back(tmpCLV2PntZ);
        CLV2CornersInnRadii.push_back(_CLV2InnerRadii.at(isurf));
        CLV2CornersOutRadii.push_back(_CLV2OuterRadii.at(isurf));
      }

      mbs._pCLV2Params = std::unique_ptr<Polycone>
        (new Polycone(CLV2CornersZ,
                      CLV2CornersInnRadii,
                      CLV2CornersOutRadii,
                      _CLV2OffsetInMu2e,
                      _CLV2MaterialName));

      if ( _MBSVersion == 3 || _MBSVersion >= 5 ){
        CLHEP::Hep3Vector _CLV2AbsOffsetInMu2e  = _MBSMOffsetInMu2e + CLHEP::Hep3Vector(0.0,0.0,_CLV2Z+0.5*CLV2totLength-_CLV2AbsHLength);
        mbs._pCLV2ABSParams = std::unique_ptr<Tube>
          (new Tube(_CLV2AbsMaterialName,
                    _CLV2AbsOffsetInMu2e,
                    0.0, //inner radius
                    _CLV2InnerRadii.at(0)-0.01, //outer radius of the plug is the same as the inner radius of the hole, leave a 10 micron gap
                    _CLV2AbsHLength));
      }
      //Shield at upstream end of the MBS to protect the Calorimeter
      if (_MBSVersion == 6 ){
        CLHEP::Hep3Vector _CalShieldRingOffsetInMu2e  = _MBSMOffsetInMu2e +
          CLHEP::Hep3Vector(0.0,0.0,_CalShieldRingZ);
        mbs._pCalShieldRingParams = std::unique_ptr<Tube>
          (new Tube(_CalShieldRingMaterialName,
                    _CalShieldRingOffsetInMu2e,
                    _CalShieldRingInnerRadius,
                    _CalShieldRingOuterRadius,
                    _CalShieldRingHLength));
      }

    } // end of Version 2 - specific constructor
  } // end of constructor for MBSMaker

  void MBSMaker::parseConfig( SimpleConfig const & _config ){

    _verbosityLevel       = _config.getInt("mbs.verbosityLevel");

    _MBSVersion           = _config.getInt("mbs.Version",1);
    _MBSCZ                = _config.getDouble("mbs.MBSCZ");

    if ( _MBSVersion >= 2 ) {
      // BSTS is the steel pipe - now a polycone
      _config.getVectorDouble("mbs.BSTSInnerRadii",_BSTSInnRadii);
      _config.getVectorDouble("mbs.BSTSOuterRadii",_BSTSOutRadii);
      _config.getVectorDouble("mbs.BSTSZLengths",_BSTSZLengths);
      _BSTSOuterRadius = 0.0;
      for ( unsigned int iR = 0; iR < _BSTSOutRadii.size(); iR++ ) {
        if ( _BSTSOutRadii[iR] > _BSTSOuterRadius ) _BSTSOuterRadius = _BSTSOutRadii[iR];
      }
    } else if ( _MBSVersion == 1 ) {
      _BSTSInnerRadius      = _config.getDouble("mbs.BSTSInnerRadius");
      _BSTSOuterRadius      = _config.getDouble("mbs.BSTSOuterRadius");
    }
    // Common to version 1 and 2
    _BSTSHLength          = _config.getDouble("mbs.BSTSHLength");
    _BSTSMaterialName     = _config.getString("mbs.BSTSMaterialName");
    _BSTSZ                = _config.getDouble("mbs.BSTSZrelCntr");


    // Version 1 specific items
    if ( _MBSVersion == 1 ) {
      _SPBSSup1InnerRadius  = _config.getDouble("mbs.SPBSSup1InnerRadius");
      _SPBSSup1OuterRadius  = _config.getDouble("mbs.SPBSSup1OuterRadius");
      _SPBSSup1HLength      = _config.getDouble("mbs.SPBSSup1HLength");
      _SPBSSup1MaterialName = _config.getString("mbs.SPBSSup1MaterialName");
      _SPBSSup1Z            = _config.getDouble("mbs.SPBSSup1ZrelCntr");
      _SPBSSup2InnerRadius  = _config.getDouble("mbs.SPBSSup2InnerRadius");
      _SPBSSup2OuterRadius  = _config.getDouble("mbs.SPBSSup2OuterRadius");
      _SPBSSup2HLength      = _config.getDouble("mbs.SPBSSup2HLength");
      _SPBSSup2MaterialName = _config.getString("mbs.SPBSSup2MaterialName");
      _SPBSSup2Z            = _config.getDouble("mbs.SPBSSup2ZrelCntr");

      _SPBSLInnerRadius     = _config.getDouble("mbs.SPBSLInnerRadius");
      _SPBSLOuterRadius     = _config.getDouble("mbs.SPBSLOuterRadius");
      _SPBSLHLength         = _config.getDouble("mbs.SPBSLHLength");
      _SPBSLMaterialName    = _config.getString("mbs.SPBSLMaterialName");
      _SPBSLZ               = _config.getDouble("mbs.SPBLZrelCntr");

      _SPBSRInnerRadius     = _config.getDouble("mbs.SPBSRInnerRadius");
      _SPBSROuterRadius     = _config.getDouble("mbs.SPBSROuterRadius");
      _SPBSRHLength         = _config.getDouble("mbs.SPBSRHLength");
      _SPBSRMaterialName    = _config.getString("mbs.SPBSRMaterialName");
      _SPBSRZ               = _config.getDouble("mbs.SPBRZrelCntr");

    }
    // Following are common to Version 1 and up
    // Outer HDPE
    _SPBSCInnerRadius     = _config.getDouble("mbs.SPBSCInnerRadius");
    _SPBSCOuterRadius     = _config.getDouble("mbs.SPBSCOuterRadius");
    _SPBSCHLength         = _config.getDouble("mbs.SPBSCHLength");
    _SPBSCminAngle        = _config.getDouble("mbs.SPBSCminAngle")*CLHEP::deg;
    _SPBSCmaxAngle        = _config.getDouble("mbs.SPBSCmaxAngle")*CLHEP::deg;
    _SPBSCMaterialName    = _config.getString("mbs.SPBSCMaterialName");
    _SPBSCZ               = _config.getDouble("mbs.SPBCZrelCntr");

    // Inner HDPE upstream
    _config.getVectorDouble("mbs.BSTCInnerRadii",_BSTCInnerRadii);
    _config.getVectorDouble("mbs.BSTCOuterRadii",_BSTCOuterRadii);
    _config.getVectorDouble("mbs.BSTCLengths",_BSTCLengths);
    _BSTCMaterialName     = _config.getString("mbs.BSTCMaterialName");
    _BSTCZ                = _config.getDouble("mbs.BSTCZrelCntr");
    // Inner HDPE downstream
    _config.getVectorDouble("mbs.BSBSInnerRadii",_BSBSInnerRadii);
    _config.getVectorDouble("mbs.BSBSOuterRadii",_BSBSOuterRadii);
    _config.getVectorDouble("mbs.BSBSLengths",_BSBSLengths);
    _BSBSMaterialName     = _config.getString("mbs.BSBSMaterialName");
    _BSBSZ                = _config.getDouble("mbs.BSBSZrelCntr");
    // HDPE end plug downstream end
    _config.getVectorDouble("mbs.CLV2InnerRadii",_CLV2InnerRadii);
    _config.getVectorDouble("mbs.CLV2OuterRadii",_CLV2OuterRadii);
    _config.getVectorDouble("mbs.CLV2Lengths",_CLV2Lengths);
    _CLV2MaterialName     = _config.getString("mbs.CLV2MaterialName");
    _CLV2Z                = _config.getDouble("mbs.CLV2ZrelCntr");

    _CLV2AbsBuild         = _config.getBool("mbs.CLV2.absorber.build",false);
    _CLV2AbsMaterialName  = _config.getString("mbs.CLV2.absorber.MaterialName","Polyethylene096");
    _CLV2AbsHLength       = _config.getDouble("mbs.CLV2.absorber.halflength",0.0);
    if ( _MBSVersion > 3 ) {
      _nHolesSt        = _config.getInt("mbs.nHolesSteel",0);
      _nHolesUP        = _config.getInt("mbs.nHolesUpstreamPoly",0);
      _nHolesDP        = _config.getInt("mbs.nHolesDownstreamPoly",0);
      _BSTSHoleXDim = _config.getDouble("mbs.steelHoleXDim",0.0);
      _BSTSHoleYDim = _config.getDouble("mbs.steelHoleYDim",0.0);
      _BSTSHoleZDim = _config.getDouble("mbs.steelHoleZDim",0.0);
      _upPolyHoleXDim = _config.getDouble("mbs.upPolyHoleXDim",0.0);
      _upPolyHoleYDim = _config.getDouble("mbs.upPolyHoleYDim",0.0);
      _upPolyHoleZDim = _config.getDouble("mbs.upPolyHoleZDim",0.0);
      _downPolyHoleXDim = _config.getDouble("mbs.downPolyHoleXDim",0.0);
      _downPolyHoleYDim = _config.getDouble("mbs.downPolyHoleYDim",0.0);
      _downPolyHoleZDim = _config.getDouble("mbs.downPolyHoleZDim",0.0);
      CLHEP::Hep3Vector moveIt(0.0,0.0,-_BSTSHLength);
      for ( int ihole = 0; ihole < _nHolesSt; ihole++ ) {
        CLHEP::Hep3Vector tempLoc;
        std::ostringstream sHoleName;
        sHoleName << "mbs.steelHoleCenter" << ihole+1;
        tempLoc = _config.getHep3Vector(sHoleName.str());
        tempLoc+=moveIt;
        _BSTSHoleCenters.push_back(tempLoc);
      }
      for ( int ihole = 0; ihole < _nHolesUP; ihole++ ) {
        CLHEP::Hep3Vector tempLoc;
        std::ostringstream pHoleName;
        pHoleName << "mbs.upPolyHoleCenter" << ihole+1;
        tempLoc = _config.getHep3Vector(pHoleName.str());
        tempLoc+=moveIt;
        _upPolyHoleCenters.push_back(tempLoc);
      }
      for ( int ihole = 0; ihole < _nHolesDP; ihole++ ) {
        CLHEP::Hep3Vector tempLoc;
        std::ostringstream pHoleName;
        pHoleName << "mbs.downPolyHoleCenter" << ihole+1;
        tempLoc =_config.getHep3Vector(pHoleName.str());
        tempLoc+=moveIt;
        _downPolyHoleCenters.push_back(tempLoc);
      }

    }
    //Get parameters for the Calorimeter shield
    if(_MBSVersion == 6) {
      _CalShieldRingInnerRadius     = _config.getDouble("mbs.CalShieldRingInnerRadius");
      _CalShieldRingOuterRadius     = _config.getDouble("mbs.CalShieldRingOuterRadius");
      _CalShieldRingHLength         = _config.getDouble("mbs.CalShieldRingHLength");
      _CalShieldRingMaterialName    = _config.getString("mbs.CalShieldRingMaterialName");
      _CalShieldRingZ               = _config.getDouble("mbs.CalShieldRingZRelCnt");
    }

  } // end parConfig function

} // namespace mu2e
