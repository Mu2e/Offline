#ifndef ProductionSolenoidGeom_PSEnclosure_hh
#define ProductionSolenoidGeom_PSEnclosure_hh

//
// $Id: PSEnclosure.hh,v 1.4 2012/06/06 19:29:30 gandr Exp $
// $Author: gandr $
// $Date: 2012/06/06 19:29:30 $
//
// Original author Andrei Gaponenko
//

#include <vector>
#include <ostream>

#include "GeomPrimitives/inc/Tube.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSEnclosureMaker;

  class PSEnclosure : virtual public Detector {

  public:

    const Tube& shell() const { return shell_; }
    const Tube& endPlate() const { return endPlate_; }
    const std::vector<Tube>& windows() const { return windows_; }

    unsigned nWindows() const { return windows_.size(); }

  private:

    friend class PSEnclosureMaker;

    // Private ctr: the class should only be constructed via PSEnclosure::PSEnclosureMaker.
    PSEnclosure(const Tube& shell, const Tube& ep)
      : shell_(shell), endPlate_(ep)
    {};

    // Or read back from persistent storage
    PSEnclosure();
    template<class T> friend class art::Wrapper;

    // The real enclosure shape is shown in docdb-2066 It is
    // approximated here by a cylinder closed with a flat end plate.
    Tube shell_;
    Tube endPlate_;
    std::vector<Tube> windows_;
  };

  std::ostream& operator<<(std::ostream& os, const PSEnclosure& pse);

}

#endif/*ProductionSolenoidGeom_PSEnclosure_hh*/
