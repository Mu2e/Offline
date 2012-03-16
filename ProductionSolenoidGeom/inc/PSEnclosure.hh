#ifndef ProductionSolenoidGeom_PSEnclosure_hh
#define ProductionSolenoidGeom_PSEnclosure_hh

//
// $Id: PSEnclosure.hh,v 1.1 2012/03/16 05:09:22 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/16 05:09:22 $
//
// Original author Andrei Gaponenko
//

#include "GeomPrimitives/inc/Tube.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

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

    // The real enclosure shape is shown in docdb-2066 It is
    // approximated here by a cylinder closed with a flat end plate.
    Tube shell_;
    Tube vacuum_;
    Tube endPlate_;
  };
}

#endif/*ProductionSolenoidGeom_PSEnclosure_hh*/
