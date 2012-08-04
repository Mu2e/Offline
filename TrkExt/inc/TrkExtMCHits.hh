//
//  $Id: TrkExtMCHits.hh,v 1.1 2012/08/04 00:22:10 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/04 00:22:10 $
//
//  Original author MyeongJae Lee
//
//
#ifndef TrkExtMCHits_HH
#define TrkExtMCHits_HH

// C++ includes.
#include <vector>
#include <string>

// Framework includes.
#include "art/Framework/Principal/Event.h"
#include "MCDataProducts/inc/StepPointMCCollection.hh"


namespace mu2e {

  class TrkExtMCHits {

  public:
    TrkExtMCHits() ;
    TrkExtMCHits(std::string g4ModuleLabel, std::string instanceName, double ethi = 50., double dist = 3.) ;
    ~TrkExtMCHits() ;

    unsigned int readMCHits (const art::Event& event);
    std::vector<StepPointMCCollection> & getClusters() {return _hitcol;}
    int getNClusters() {return _hitcol.size();}
    int getNHits() {return _nhits;}
    void setEnergyThreshold (double e) { _eth = ((e>0)?e:0.);}
    void setClusterHitDistance (double d) { _dist = ((d>0)?d:1.);}

  private:
    std::string _g4ModuleLabel;
    std::string _instanceName;
    std::vector<StepPointMCCollection> _hitcol;
    int _nhits;
    double _eth, _dist;


  };



} // end namespace mu2e


#endif
