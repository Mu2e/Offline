#ifndef TrackerGeom_Tracker_hh
#define TrackerGeom_Tracker_hh
//
// Pure virtual base class that will used by both LTracker or TTracker.
//
// $Id: Tracker.hh,v 1.3 2011/05/17 15:41:37 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:37 $
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

#endif /* TrackerGeom_Tracker_hh */
