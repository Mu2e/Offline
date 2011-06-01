//=============================================================================
//
// An Incident is a collection of StawHitSynopsis corresponding to some
// limited time window within the event.  
//
// $Id: Incident.cc,v 1.1 2011/06/01 23:29:47 mf Exp $
// $Author: mf $
// $Date: 2011/06/01 23:29:47 $
//
// Original author: Mark Fischler
//
//=============================================================================

// C++ includes
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>

// Framework includes.
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "art/Persistency/Provenance/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"

// Mu2e includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "Mu2eUtilities/inc/LineSegmentPCA.hh"
#include "Mu2eUtilities/inc/SimParticlesWithHits.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"
#include "DataProducts/inc/DPIndexVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "TrackerGeom/inc/StrawId.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "EarlyPatRec/inc/Incident.hh"

namespace mu2e {

// Incident:  
//
// An Incident is a collection of StrawHits, all in one event
// with timings such that their time of the pulse reaching the wire
// lies, for all included hits, in a given time window.  The idea is
// that this window is chosen such that all the hits are consistent 
// with being part of the same track.  

Incident::Incident(Tracker const & tracker, 
                   std::vector<StrawHit const *> & strawHitPtr,
                   double t1, double t2 )
    : hits()
    , _tt(&tracker)
{
  typedef std::vector<StrawHit const*>::const_iterator Sciter;
  for (Sciter s = strawHitPtr.begin(); s != strawHitPtr.end(); ++s) {
    if (strawHitIsInTimeWindow( **s, t1, t2 )) {
      double t = rough_t0_star (**s);
      hits.push_back (StrawHitSynopsis(*_tt, **s, t));
    }
  }
  populateDetectorEventModel();
} // Incident ctor

bool Incident::strawHitIsInTimeWindow( 
     StrawHit const& h, double t1, double t2 ) const
{
  double t = rough_t0_star (h);
  return ( (t >= t1) && (t <= t2) ); 
} // strawHitIsInTimeWindow

double Incident::rough_t0_star ( StrawHit const& h ) const
{
  const double c = 299.8; // mm/ns
  const double wirePropogationSpeed = 200.0; // mm/ns
  const double betaC = c;
  const double averageSecTheta = 1.635; // based on avg tan theta of 1.294
  double strawHalfLength = _tt->getStraw(h.strawIndex()).getHalfLength();
  double strawZ =  _tt->getStraw(h.strawIndex()).getMidPoint().z();
  double t = h.time();
  double delta_t = h.dt();
  return t - 0.5*delta_t + strawHalfLength/wirePropogationSpeed
         + strawZ*averageSecTheta/betaC;
}

void Incident::populateDetectorEventModel() 
{
    // TODO -- we want the hits structure to be a bit richer than 
    // this; at least sorted by detector element etc.
}

} // end of namespace mu2e
