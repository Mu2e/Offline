//
// Make a ExtinctionMonitor.
//
// $Id: ExtMonUCIMaker.cc,v 1.1 2011/12/10 00:16:15 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/10 00:16:15 $

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

    void ExtMonMaker::MakeTofs()
    {

    }
  }
}
