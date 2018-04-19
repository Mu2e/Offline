#ifndef EarlyPatRecPR_SubEvent_hh
#define EarlyPatRecPR_SubEvent_hh

//=============================================================================
//
// An SubEvent is a collection of StawHitSynopsis corresponding to some
// limited time window within the event.  
//
// $Id: SubEvent.hh,v 1.1 2011/06/05 23:11:35 mf Exp $
// $Author: mf $
// $Date: 2011/06/05 23:11:35 $
//
// Original author: Mark Fischler
//
//=============================================================================

// C++ includes
#include <iostream>
#include <vector>
#include <utility>

// Framework includes.

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitSynopsis.hh"

namespace mu2e {

  // SubEvent:  
  //
  // An SubEvent is a collection of StrawHits, all in one event
  // with timings such that their time of the pulse reaching the wire
  // lies, for all included hits, in a given time window.  The idea is
  // that this window is chosen such that all the hits are consistent 
  // with being part of the same track.  

  class SubEvent 
  {
  public:
    // ctors 
    // SubEventsMaker likes to provide just time range, filling in hits later
    SubEvent( double t_start, double t_end ) 
      : t1(t_start)
      , t2(t_end) 
      , hits()
    { }
    
    // Default ctor    
    SubEvent() : t1(0.0), t2(0.0), hits() { }
    
    // filling the data
    void push_back (  StrawHitSynopsis const & s ) 
    {
      hits.push_back(s);
    }
    
    // main consuming uses of this class:  
    
    std:: vector<StrawHitSynopsis> const & getHits() const {return hits;}
    unsigned int nHits() const { return hits.size(); }
    std::pair <double, double> timeRange() const 
    { 
      return std::pair<double, double> (t1,t2); 
    } 
               
  private:
    // private data
    double t1;
    double t2;
    std:: vector<StrawHitSynopsis> hits;
    
  }; // SubEvent class definition

  // TODO -- either decide not to provide operator<< (os,s) or inline one here. 

} // end of mu2e namespace

#endif // EarlyPatRecPR_SubEvent_hh
