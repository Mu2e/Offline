#ifndef TRACKER_HH
#define TRACKER_HH
//
// Pure virtual base class that will used by both LTracker or TTracker.
//
// $Id: Tracker.hh,v 1.2 2010/04/18 00:31:55 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/18 00:31:55 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <vector>

#include "GeometryService/inc/Detector.hh"
#include "TrackerGeom/inc/Straw.hh"

namespace mu2e {

  class Tracker: public Detector{

  public:
    // No constructor on purpose.
    ~Tracker();
    
    // Compiler generated copy and assignment constructors
    // should be OK.

    // Check for legal identifiers; add these later.
    //virtual isLegal(DeviceId d) const = 0;    
    //virtual bool isLegal(const SectorId& sid) const = 0;
    //vitrual bool isLegal(const LayerId& lid ) const = 0;
    //virtual bool isLegal(const StrawId& sid) const =0;

    virtual const Straw& getStraw ( const StrawId& sid ) const=0;
    virtual const Straw& getStraw ( StrawIndex i ) const=0;
    virtual const std::deque<Straw>& getAllStraws() const=0;
    virtual const std::vector<StrawDetail>& getStrawDetails() const=0;

  };

} //namespace mu2e

#endif
