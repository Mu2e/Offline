//
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
#include "MCDataProducts/inc/StepPointMC.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "cetlib/map_vector.h"

namespace mu2e {

  class TrkExtMCHits {

  public:
    TrkExtMCHits() ;
    TrkExtMCHits(const art::Event& event, std::string g4ModuleLabel, std::string instanceName, cet::map_vector_key simid, double distcut = 1.) ;
    ~TrkExtMCHits() ;

    const std::vector<StepPointMCCollection> & getClusters() const {return _hitcol;}
    unsigned int getNClusters() const {return _hitcol.size();}
    const StepPointMCCollection & getCluster(unsigned int clust) const {return _hitcol[clust]; }
    unsigned int getNHit (unsigned int clust) const { return _hitcol[clust].size(); }
    const StepPointMC & getHit (unsigned int clust, unsigned int hit) const { return (_hitcol[clust])[hit]; }
    CLHEP::Hep3Vector momentum (unsigned int clust) const ;
    CLHEP::Hep3Vector position (unsigned int clust) const ;
    double time (unsigned int clust) const;
    double deltap (unsigned int clust) const;
    double eDep (unsigned int clust) const;
    double ionizingEdep (unsigned int clust) const;
    double nonIonizingEdep (unsigned int clust) const;

  private:
    double interpolate3 (double z, double x1, double x2, double x3, double y1, double y2, double y3) const; 
    double interpolate2 (double z, double x1, double x2, double y1, double y2) const ;

  private:
    std::string _g4ModuleLabel;
    std::string _instanceName;
    std::vector<StepPointMCCollection> _hitcol;
    double _distcut;


  };



} // end namespace mu2e


#endif
