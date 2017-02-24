//
//  $Id: TrkExtTrajPoint.cc,v 1.2 2013/02/07 02:09:47 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2013/02/07 02:09:47 $
//
//  Original author MyeongJae Lee
//
//

// C++ includes.
#include <iostream>
#include <string>
//#include <stdexcept>
#include <sstream>

// Framework includes.
#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"

using namespace CLHEP;

using namespace std;

namespace mu2e {

  TrkExtTrajPoint::TrkExtTrajPoint() :
    _x(0),
    _p(0),
    _volid(0),
    _trajPtId(-1),
    _fl(0),
    _ft(0)
  {
    _cov =  HepMatrix(6,6,0);
  }

  TrkExtTrajPoint::~TrkExtTrajPoint() 
  { 
  }

  TrkExtTrajPoint::TrkExtTrajPoint(int id, Hep3Vector & x, Hep3Vector & p, HepMatrix & cov, int volid, double fl, double ft) :
    _x(x),
    _p(p),
    _volid(volid),
    _trajPtId(id),
    _fl(fl),
    _ft(ft)
  {
    if (cov.num_row() !=6 ||cov.num_col() !=6) {
      throw cet::exception("ARGUMENT") << "TrkExtTrajPoint Error: invalid dimension for cov matrix" << endl;
    }
    else {
      _cov = cov;
    }
  }

  TrkExtTrajPoint::TrkExtTrajPoint (const TrkExtTrajPoint & dt) {
    _x = dt._x;
    _p = dt._p;
    _cov = dt._cov;
    _volid = dt._volid;
    _trajPtId = dt._trajPtId;
    _fl = dt._fl;
    _ft = dt._ft;
  }

  TrkExtTrajPoint & TrkExtTrajPoint::operator= (const TrkExtTrajPoint & dt) {
    _x = dt._x;
    _p = dt._p;
    _cov = dt._cov;
    _volid = dt._volid;
    _trajPtId = dt._trajPtId;
    _fl = dt._fl;
    _ft = dt._ft;
    return (*this);
  }

  TrkExtTrajPoint::TrkExtTrajPoint(int id, CLHEP::HepVector & r, int volid, double fl, double ft) {
    _x.set(r[0], r[1], r[2]);
    _p.set(r[3], r[4], r[5]);
    _volid = volid;
    _trajPtId = id;
    _fl = fl;
    _ft = ft;
  }

  HepVector TrkExtTrajPoint::vector() {
    HepVector xp(6,0);
    xp[0] = x();
    xp[1] = y();
    xp[2] = z();
    xp[3] = px();
    xp[4] = py();
    xp[5] = pz();
    return xp;
  }

  void TrkExtTrajPoint::setCovariance(CLHEP::HepMatrix & c) {
    _cov = c;
  }

  CLHEP::Hep3Vector TrkExtTrajPoint::positionError () {
    Hep3Vector error(ex(), ey(), ez());
    return error;
  }

  CLHEP::Hep3Vector TrkExtTrajPoint::momentumError () {
    Hep3Vector error(epx(), epy(), epz());
    return error;
  }


} // end namespace mu2e


