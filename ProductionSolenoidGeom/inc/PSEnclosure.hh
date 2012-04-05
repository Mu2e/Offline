#ifndef ProductionSolenoidGeom_PSEnclosure_hh
#define ProductionSolenoidGeom_PSEnclosure_hh

//
// $Id: PSEnclosure.hh,v 1.2 2012/04/05 18:43:39 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/05 18:43:39 $
//
// Original author Andrei Gaponenko
//

#include "GeomPrimitives/inc/Tube.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSEnclosureMaker;

  class PSEnclosure : virtual public Detector {

  public:

    const Tube& shell() const { return shell_; }
    const Tube& vacuum() const { return vacuum_; }
    const Tube& endPlate() const { return endPlate_; }

  private:

    friend class PSEnclosureMaker;

    // Private ctr: the class should only be constructed via PSEnclosure::PSEnclosureMaker.
    PSEnclosure(const Tube& shell, const Tube& vac, const Tube& ep)
      : shell_(shell), vacuum_(vac), endPlate_(ep)
    {};

    // Or read back from persistent storage
    PSEnclosure();
    template<class T> friend class art::Wrapper;

    // The real enclosure shape is shown in docdb-2066 It is
    // approximated here by a cylinder closed with a flat end plate.
    Tube shell_;
    Tube vacuum_;
    Tube endPlate_;
  };
}

#endif/*ProductionSolenoidGeom_PSEnclosure_hh*/
