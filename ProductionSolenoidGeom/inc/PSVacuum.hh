#ifndef ProductionSolenoidGeom_PSVacuum_hh
#define ProductionSolenoidGeom_PSVacuum_hh

//
//
// Original author Andrei Gaponenko
//

#include "Offline/GeomPrimitives/inc/Cone.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSVacuumMaker;

  class PSVacuum : virtual public Detector {

  public:

    const Tube& vacuum() const { return vacuum_; }

    // Properties of the vacuum inside the PS
    double vacuumPressure()          const { return _vacuumPressure;     }
    std::string vacuumG4Material() const { return _vacuumG4Material;   }

    // Conical hard connect between the PS bore tube and the TSu ring plate.
    // tsInterconnect() is only meaningful when hasTSInterconnect() is true.
    bool hasTSInterconnect()         const { return hasTSInterconnect_; }
    const Cone& tsInterconnect()     const { return tsInterconnect_; }

  private:

    friend class PSVacuumMaker;

    // Private ctr: the class should only be constructed via PSVacuum::PSVacuumMaker.
    explicit PSVacuum(const Tube& vac)
      : vacuum_(vac)
    {}

    // Or read back from persistent storage
    PSVacuum() {}
    template<class T> friend class art::Wrapper;

    // The real enclosure shape is shown in docdb-2066 It is
    // approximated here by a cylinder closed with a flat end plate.
    Tube vacuum_;

    // Properties of the vacuum inside the PS.
    double _vacuumPressure;
    std::string _vacuumG4Material;

    bool hasTSInterconnect_ = false;
    Cone tsInterconnect_;

  };

}

#endif/*ProductionSolenoidGeom_PSVacuum_hh*/
