//
// Make a ExtinctionMonitor.
//
// $Id: ExtMonUCIMaker.cc,v 1.2 2011/12/27 21:53:01 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/27 21:53:01 $

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCIMaker.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCI.hh"

namespace mu2e {

  namespace ExtMonUCI {

    ExtMonMaker::ExtMonMaker(const SimpleConfig& config) 
      : _det(0)
    {
      std::vector<double> envelopeParams;
      std::vector<double> envelopeOrigin;
      config.getVectorDouble("extmon_uci.envelopeParams", envelopeParams, 3);
      config.getVectorDouble("extmon_uci.envelopeOrigin", envelopeOrigin, 3);	
      _det.reset(new ExtMon(envelopeParams, envelopeOrigin));
      _det->_rotation = CLHEP::HepRotation::IDENTITY; 

      _det->_nCols = config.getInt("extmon_uci.nCols");
      config.getVectorDouble("extmon_uci.colOuterHalfLengths", _det->_colOuterHalfLengths, _det->_nCols*3);
      config.getVectorDouble("extmon_uci.colInnerHalfLengths", _det->_colInnerHalfLengths, _det->_nCols*3);
      config.getVectorDouble("extmon_uci.colPosition1", _det->_colPosition1, _det->_nCols*3);
      config.getVectorDouble("extmon_uci.colPosition2", _det->_colPosition2, _det->_nCols*3);

      MakeCols();

      _det->_nMags = config.getInt("extmon_uci.nMags");
      config.getVectorDouble("extmon_uci.magOuterHalfLengths", _det->_magOuterHalfLengths, _det->_nMags*3);
      config.getVectorDouble("extmon_uci.magInnerHalfLengths", _det->_magInnerHalfLengths, _det->_nMags*3);
      config.getVectorDouble("extmon_uci.magPosition1", _det->_magPosition1, _det->_nMags*3);
      config.getVectorDouble("extmon_uci.magPosition2", _det->_magPosition2, _det->_nMags*3);

      MakeMags();

      _det->_nTofStations = config.getInt("extmon_uci.nTofStations");
      _det->_nTofSegments = config.getInt("extmon_uci.nTofSegments");
      config.getVectorDouble("extmon_uci.tofHalfLengths", _det->_tofHalfLengths, _det->_nTofStations*3);
      config.getVectorDouble("extmon_uci.tofPosition", _det->_tofPosition, _det->_nTofStations*3);

      MakeTofs();
    }

    void ExtMonMaker::MakeCols()
    {
      for (int iCol = 0; iCol < (int)_det->_nCols; iCol++)
      {
        //cout << "in col " << iCol << endl;
        _det->_cols.push_back( ExtMonCol(iCol) );
	ExtMonCol &col = _det->_cols.back(); 

        col._position1 = CLHEP::Hep3Vector(_det->_colPosition1[3*iCol], _det->_colPosition1[3*iCol+1], _det->_colPosition1[3*iCol+2]);
        col._position2 = CLHEP::Hep3Vector(_det->_colPosition2[3*iCol], _det->_colPosition2[3*iCol+1], _det->_colPosition2[3*iCol+2]);
        //cout << "pos1 " << col._position1 << endl;
        //cout << "pos2 " << col._position2 << endl;

        col._origin = 0.5*(col._position1 + col._position2);
        col._originLocal = _det->mu2eToExtMonPoint(col._origin); 
        //cout << "origin " << col._origin << endl;
        //cout << "originLocal " << col._originLocal << endl;

        const CLHEP::Hep3Vector interZ(0.0, col._position1.y()-col._position2.y(), col._position1.z()-col._position2.z());
        const CLHEP::Hep3Vector newY = interZ.cross( CLHEP::Hep3Vector(1.0, 0.0, 0.0) ).unit();
        const CLHEP::Hep3Vector newZ = CLHEP::Hep3Vector(col._position1.x() - col._position2.x(),
                                                   col._position1.y() - col._position2.y(),
                                                   col._position1.z() - col._position2.z()).unit();
        const CLHEP::Hep3Vector newX = newY.cross(newZ);
        col._rotation = CLHEP::HepRotation::IDENTITY;
        col._rotation.rotateAxes( newX, newY, newZ );
        col._rotation.invert(); 
        //col._rotation.print(cout);

        for (int iDim = 0; iDim < 3; iDim++)
        {
          col._paramsOuter.push_back(_det->_colOuterHalfLengths[3*iCol+iDim]);
          col._paramsInner.push_back(_det->_colInnerHalfLengths[3*iCol+iDim]);
        }
      }
    }

    void ExtMonMaker::MakeMags()
    {
      for (int iMag = 0; iMag < (int)_det->_nMags; iMag++)
      {
        //cout << "in mag " << iMag << endl;
        _det->_mags.push_back( ExtMonMag(iMag) );
        ExtMonMag &mag = _det->_mags.back();

        mag._position1 = CLHEP::Hep3Vector(_det->_magPosition1[3*iMag], _det->_magPosition1[3*iMag+1], _det->_magPosition1[3*iMag+2]);
        mag._position2 = CLHEP::Hep3Vector(_det->_magPosition2[3*iMag], _det->_magPosition2[3*iMag+1], _det->_magPosition2[3*iMag+2]);
        //cout << "pos1 " << mag._position1 << endl;
        //cout << "pos2 " << mag._position2 << endl;

        mag._origin = 0.5*(mag._position1 + mag._position2);
        mag._originLocal = _det->mu2eToExtMonPoint(mag._origin);
        //cout << "origin " << mag._origin << endl;
        //cout << "originLocal " << mag._originLocal << endl;

        const CLHEP::Hep3Vector interZ(0.0, mag._position1.y()-mag._position2.y(), mag._position1.z()-mag._position2.z());
        const CLHEP::Hep3Vector newY = interZ.cross( CLHEP::Hep3Vector(1.0, 0.0, 0.0) ).unit();
        const CLHEP::Hep3Vector newZ = CLHEP::Hep3Vector(mag._position1.x() - mag._position2.x(),
                                                   mag._position1.y() - mag._position2.y(),
                                                   mag._position1.z() - mag._position2.z()).unit();
        const CLHEP::Hep3Vector newX = newY.cross(newZ);
        mag._rotation = CLHEP::HepRotation::IDENTITY;
        mag._rotation.rotateAxes( newX, newY, newZ );
        mag._rotation.invert();
        //mag._rotation.print(cout);

        for (int iDim = 0; iDim < 3; iDim++)
        {
          mag._paramsOuter.push_back(_det->_magOuterHalfLengths[3*iMag+iDim]);
          mag._paramsInner.push_back(_det->_magInnerHalfLengths[3*iMag+iDim]);
        }
      }
    }

    void ExtMonMaker::MakeTofs()
    {
      for (int iTofSta = 0; iTofSta < (int)_det->_nTofStations; iTofSta++)
      {
        //cout << "in tof station " << iTofSta << endl;
        CLHEP::Hep3Vector tofStaCenter(_det->_tofPosition[3*iTofSta], _det->_tofPosition[3*iTofSta+1], _det->_tofPosition[3*iTofSta+2]);

        for (int iTofSeg = 0; iTofSeg < (int)_det->_nTofSegments; iTofSeg++)
        {
          //cout << "  in tof segment " << iTofSeg << endl;
 
          _det->_tofs.push_back( ExtMonTof(iTofSta, iTofSeg) );
          ExtMonTof &tof = _det->_tofs.back();

          tof.setId( iTofSta*_det->_nTofSegments + iTofSeg );

          CLHEP::Hep3Vector tofSegCenter( (iTofSeg-0.5*(_det->_nTofSegments-1))*2.0*_det->_tofHalfLengths[3*iTofSta], 0.0, 0.0);
          
          tof._origin = tofStaCenter + tofSegCenter;
          tof._originLocal = _det->mu2eToExtMonPoint(tof._origin);
          //cout << "  origin " << tof._origin << endl;
          //cout << "  originLocal " << tof._originLocal << endl;

          tof._rotation = CLHEP::HepRotation::IDENTITY;
          //tof._rotation.print(cout);

          for (int iDim = 0; iDim < 3; iDim++)
          {
            tof._params.push_back(_det->_tofHalfLengths[3*iTofSta+iDim]);
          }
        } // end of Tof Segment
      } // end of Tof Station
    }


  }
}
