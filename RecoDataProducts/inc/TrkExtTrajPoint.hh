//
//  $Id: TrkExtTrajPoint.hh,v 1.2 2013/02/07 02:09:47 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2013/02/07 02:09:47 $
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtTrajPoint_HH
#define TrkExtTrajPoint_HH

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "GeneralUtilities/inc/safeSqrt.hh"

namespace mu2e {



  class TrkExtTrajPoint {

  public:

    TrkExtTrajPoint() ;
    TrkExtTrajPoint(int id, CLHEP::Hep3Vector & x, CLHEP::Hep3Vector & p, CLHEP::HepMatrix & cov, int volid, double fl, double ft) ;
    TrkExtTrajPoint(int id, CLHEP::HepVector & r, int volid, double fl, double ft);
    ~TrkExtTrajPoint() ; 
    TrkExtTrajPoint (const TrkExtTrajPoint & dt);
    TrkExtTrajPoint & operator = (const TrkExtTrajPoint & dt);

    CLHEP::Hep3Vector const & position () const { return _x; }
    CLHEP::Hep3Vector const & momentum () const { return _p; }
    CLHEP::HepMatrix const & covariance () const { return _cov; }
    CLHEP::Hep3Vector positionError () ;
    CLHEP::Hep3Vector momentumError () ;
    double x() const { return _x.x(); }
    double y() const { return _x.y(); }
    double z() const { return _x.z(); }
    double rho() const { return _x.rho(); }
    double px() const { return _p.x(); }
    double py() const { return _p.y(); }
    double pz() const { return _p.z(); }
    double p() const { return _p.mag(); }
    double ex() const { return safeSqrt(_cov[0][0]); }
    double ey() const { return safeSqrt(_cov[1][1]); }
    double ez() const { return safeSqrt(_cov[2][2]); }
    double epx() const { return safeSqrt(_cov[3][3]); }
    double epy() const { return safeSqrt(_cov[4][4]); }
    double epz() const { return safeSqrt(_cov[5][5]); }
    double covxx() const { return _cov[0][0]; }
    double covxy() const { return _cov[0][1]; }
    double covxz() const { return _cov[0][2]; }
    double covyy() const { return _cov[1][1]; }
    double covyz() const { return _cov[1][2]; }
    double covzz() const { return _cov[2][2]; }
    double covpxpx() const { return _cov[3][3]; }
    double covpxpy() const { return _cov[3][4]; }
    double covpxpz() const { return _cov[3][5]; }
    double covpypy() const { return _cov[4][4]; }
    double covpypz() const { return _cov[4][5]; }
    double covpzpz() const { return _cov[5][5]; }
    int volumeId() const { return _volid; }
    int trajPointId() const {return _trajPtId;}
    double flightLength() const {return _fl;}
    double flightTime() const { return _ft; }

    CLHEP::HepVector vector() ;
    void setCovariance(CLHEP::HepMatrix & c) ;
    void setTrajPointId(int i) { _trajPtId = i; }
    void scaleMomentum(double sf) { _p *= sf; }

  private:
    void calculateHelixParameter();

  private:
    CLHEP::Hep3Vector _x;
    CLHEP::Hep3Vector _p;
    CLHEP::HepMatrix  _cov;
    int _volid;
    int _trajPtId;
    double _fl;  // flight length in mm
    double _ft;  // flight time in ns

  };



} // end namespace mu2e


#endif
