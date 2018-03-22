#ifndef TrackerGeom_Tracker_hh
#define TrackerGeom_Tracker_hh
//
// Pure virtual base class that is used by the TTracker and was used
// by the LTracker, which has since been removed.
//
// $Id: Tracker.hh,v 1.8 2013/03/26 23:28:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/26 23:28:23 $
//
// Original author Rob Kutschke
//

#include <deque>
#include <array>
#include <vector>

#include "Mu2eInterfaces/inc/Detector.hh"
#include "TrackerGeom/inc/Straw.hh"

namespace mu2e {

  class Tracker: virtual public Detector{

  public:
    // No constructor on purpose.
    virtual ~Tracker(){}

    // Compiler generated copy and assignment constructors
    // should be OK.

    // Check for legal identifiers; add these later.
    //virtual isLegal(PlaneId d) const = 0;
    //virtual bool isLegal(const PanelId& pnlid) const = 0;
    //vitrual bool isLegal(const LayerId& layid ) const = 0;
    //virtual bool isLegal(const StrawId& strid) const =0;
    constexpr static int _nttstraws = StrawId::_nplanes *
                                      StrawId::_npanels *
                                      StrawId::_nstraws;

    virtual const Straw& getStraw ( const StrawId& strid ) const=0;
    virtual const Straw& getStraw ( StrawIndex i ) const=0;
    virtual const std::array<Straw,_nttstraws>& getAllStraws() const=0;
    virtual const std::vector<StrawDetail>& getStrawDetails() const=0;
    virtual uint16_t nPlanes() const=0;
    virtual uint16_t nStraws() const=0;

  };

} //namespace mu2e

#endif /* TrackerGeom_Tracker_hh */
