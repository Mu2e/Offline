//
//
//  Original author MyeongJae Lee
//
// *Note on point information*
// TrkExtTrajPoint have the point information AFTER the extrapolation step. 
// When the extrapolation crosses a volume change, the "after" step of 
// extrapolation locates nearest to the boundary across the bounday. 
// (The "before" step certainly locates before the boundary.) Therefore,
// the TrkExtTrajPoint information for volume crossing is the information
// after the bounday, but as nearest as possible to the boundary 
// (i.e. the changed volume)
//
// The _pa/sthitidx vectors records the trajPointID of 
// these bounday-changing points.
// Therefore, the "first" TrkExtTrajPoint locates inside the PA/ST volume,
// where the "second" TrkExtTrajPoint locates outside of the PA/ST.
//
// By calling addPA/STHit() in the extrapolation module, the vector of
// _pa/sthitidx are filled with the "first" and "second" trajPointIDs of
// PA/ST hits.
// These are converted to the index of _pt and stored in the vector of 
// _ptidx_pa/st, by calling makePASTHitTable() from the extrapolation module.
//
// Note that the index of _pt and trajPointID are not same. The trajPointID 
// increase by one at each extrapolation steps. The index of _pt increases
// when they are recorded. The recording step is not the same with the 
// extrapolation step, and can be adjusted from the configuration. The _pt
// is filled and recorded when the flight length reaches to recording step, 
// or at the start and stop of extrapolation, or at the volume change, 
// or inside the PA/ST. 
// These are properly implemented in the extrapolation module. 


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
    void setHepid (int hepid) { _hepid = hepid; }

    // other accessor
    int hepid() const {return _hepid; }
    int exitCode() const { return _ret; }
    double flightLength() const { return back().flightLength(); }
    double flightTime() const { return back().flightTime(); }
//    std::vector<std::pair<unsignedi int,unsigned int> > & paHitIndex () { return _pahitidx; }
//    std::vector<std::pair<unsigned int,unsigned int> > & stHitIndex () { return _sthitidx; }
    int id() const { return 0; }   //TODO. used in EventDisplay/DataInterface.cc. 
    //Note on id() : Better to replace to correponding trkId. 
    //It's not possilbe since there is no connection to KalRep. 
    unsigned int getNPAHits () const { return _pahitidx.size(); }
    unsigned int getNSTHits () const { return _sthitidx.size(); }
    const TrkExtTrajPoint & getFirstPAHit (unsigned int i) const ;
    const TrkExtTrajPoint & getFirstSTHit (unsigned int i) const ;
    TrkExtTrajPoint getMeanPAHit (unsigned int i) const ;
    TrkExtTrajPoint getMeanSTHit (unsigned int i) const ;
    double getDeltapPA (unsigned int i) const { return _deltap_pa[i]; }
    double getDeltapST (unsigned int i) const { return _deltap_st[i]; }
    double getDeltapPA (void) const ;
    double getDeltapST (void) const ;

    std::vector<TrkExtTrajPoint> getPointsAtZ (double z, unsigned int idx1 = 0, unsigned int idx2 = 0) const;
    //CLHEP::TrkExtTrajPoint getPointAtFl (double fl) const;  //TODO

  private:
    unsigned int findPASTHit (unsigned int idx, unsigned int start);
    unsigned int findNeighborsZ (double z, std::vector<unsigned int> & idx, unsigned int idx1 =0, unsigned int idx2=0) const; 
    //int findNeighborsFl (double z, unsigned int idx) const; 
//    double interpolate2 (double z, double x1, double x2, double y1, double y2) const;
    double interpolate3 (double z, double x1, double x2, double x3, double y1, double y2, double y3) const;
    TrkExtTrajPoint interpolatePoint(unsigned int first, unsigned int second) const;

  private:
    std::vector<TrkExtTrajPoint> _pt; // ext point info
    int _hepid;
    int _ret; // exit code
    std::vector<std::pair<unsigned int, unsigned int> > _pahitidx;  // trajPointIdx pair for entering and exiting PA/ST
    std::vector<std::pair<unsigned int, unsigned int> > _sthitidx;

    std::vector<std::pair<unsigned int, unsigned int> > _ptidx_pa;  // index in _pt corresponding to _pa/sthitidx
    std::vector<std::pair<unsigned int, unsigned int> > _ptidx_st;
    std::vector<double> _deltap_pa;
    std::vector<double> _deltap_st;

  };



} // end namespace mu2e


#endif
