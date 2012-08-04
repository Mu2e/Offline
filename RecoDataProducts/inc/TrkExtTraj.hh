//
//  $Id: TrkExtTraj.hh,v 1.1 2012/08/04 00:16:02 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/04 00:16:02 $
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtTraj_HH
#define TrkExtTraj_HH

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"

namespace mu2e {



  class TrkExtTraj {

  public:
    TrkExtTraj(); 
    TrkExtTraj(const TrkExtTrajPoint & trajPoint);
    TrkExtTraj (const TrkExtTraj & dt) ; 
    ~TrkExtTraj() {;} 
    TrkExtTraj & operator = (const TrkExtTraj & dt) ;

    // vector-like modifier
    void push_back (const TrkExtTrajPoint & trajPoint) ;
    void clear() ;

    // vector-like accessor
    unsigned int size () const { return _pt.size(); }
    TrkExtTrajPoint  & operator [] (int i) ;
    TrkExtTrajPoint    operator [] (int i) const;
    TrkExtTrajPoint & front () { return _pt.front(); }
    TrkExtTrajPoint   front () const { return _pt.front(); }
    TrkExtTrajPoint & back  () { return _pt.back(); }
    TrkExtTrajPoint   back  () const { return _pt.back(); }

    // other modifier
    void setExitCode(int ret) {_ret = ret;}
    void addPAHit (unsigned int idx) { _pahitidx.push_back(idx); }
    void addSTHit (unsigned int idx) { _sthitidx.push_back(idx); }

    // other accessor
    int exitCode() const { return _ret; }
    double flightLength() const { return back().flightLength(); }
    double flightTime() const { return back().flightTime(); }
    std::vector<unsigned int> & paHitIndex () { return _pahitidx; }
    std::vector<unsigned int> & stHitIndex () { return _sthitidx; }
    int id() const { return 0; }  //TODO

    /* // TODO
    CLHEP::Hep3Vector position (double fl) ;
    CLHEP::Hep3Vector momentum (double fl) ;
    CLHEP::HepMatrix  covariance (double fl);
    bool positionMomentumCovariance (double fl, CLHEP::Hep3Vector po, CLHEP::HEp3Vector::mo, CLHEP::HepMatrix co);
*/

  private:
    std::vector<TrkExtTrajPoint> _pt; // ext point info
    int _ret; // exit code
    std::vector<unsigned int> _pahitidx;
    std::vector<unsigned int> _sthitidx;
  };



} // end namespace mu2e


#endif
