//
//  $Id: TrkExtTraj.hh,v 1.2 2012/10/23 00:25:07 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/10/23 00:25:07 $
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtTraj_HH
#define TrkExtTraj_HH

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"
#include <utility>

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
    void addPAHit (unsigned int idx1, unsigned int idx2) { _pahitidx.push_back(std::make_pair(idx1,idx2)); }
    void addSTHit (unsigned int idx1, unsigned int idx2) { _sthitidx.push_back(std::make_pair(idx1,idx2)); }
    void makePASTHitTable (void);

    // other accessor
    int exitCode() const { return _ret; }
    double flightLength() const { return back().flightLength(); }
    double flightTime() const { return back().flightTime(); }
//    std::vector<std::pair<unsignedi int,unsigned int> > & paHitIndex () { return _pahitidx; }
//    std::vector<std::pair<unsigned int,unsigned int> > & stHitIndex () { return _sthitidx; }
    int id() const { return 0; }  //TODO. used in EventDisplay/DataInterface.cc
    unsigned int getNPAHits () const { return _pahitidx.size(); }
    unsigned int getNSTHits () const { return _sthitidx.size(); }
    TrkExtTrajPoint & getPAHit (unsigned int i) { return _pt[_ptidx_pa[i].first]; }
    TrkExtTrajPoint & getSTHit (unsigned int i) { return _pt[_ptidx_st[i].first]; }
    double getDeltapPA (unsigned int i) const { return _deltap_pa[i]; }
    double getDeltapST (unsigned int i) const { return _deltap_st[i]; }
    double getDeltapPA (void) ;
    double getDeltapST (void) ;


  private:
    unsigned int findPASTHit (unsigned int idx, unsigned int start);
//    void checkPASTTable(void);

    /* // TODO
    CLHEP::Hep3Vector position (double fl) ;
    CLHEP::Hep3Vector momentum (double fl) ;
    CLHEP::HepMatrix  covariance (double fl);
    bool positionMomentumCovariance (double fl, CLHEP::Hep3Vector po, CLHEP::HEp3Vector::mo, CLHEP::HepMatrix co);
*/

  private:
    std::vector<TrkExtTrajPoint> _pt; // ext point info
    int _ret; // exit code
    std::vector<std::pair<unsigned int, unsigned int> > _pahitidx;  // trajPointIdx pair for entering and exiting PA/ST
    std::vector<std::pair<unsigned int, unsigned int> > _sthitidx;

    //internal variables
    std::vector<std::pair<unsigned int, unsigned int> > _ptidx_pa;  // index in _pt corresponding to _pa/sthitidx
    std::vector<std::pair<unsigned int, unsigned int> > _ptidx_st;
    std::vector<double> _deltap_pa;
    std::vector<double> _deltap_st;

  };



} // end namespace mu2e


#endif
