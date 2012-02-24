#ifndef TrackerGeom_Tracker_hh
#define TrackerGeom_Tracker_hh
//
// Pure virtual base class that will used by both LTracker or TTracker.
//
// $Id: Tracker.hh,v 1.6 2012/02/24 20:55:48 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 20:55:48 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <vector>

#include "Mu2eInterfaces/inc/Detector.hh"
#include "TrackerGeom/inc/Straw.hh"

namespace mu2e {

  class Tracker: virtual public Detector{

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
