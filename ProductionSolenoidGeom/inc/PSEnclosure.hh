#ifndef ProductionSolenoidGeom_PSEnclosure_hh
#define ProductionSolenoidGeom_PSEnclosure_hh

//
//
// Original author Andrei Gaponenko
//
// Dec 2017 - Dave (Louisville) Brown adds version 2, which makes the
// PSEnclosure a conical frustrum instead of a tube.  This is per as-built
// geometry - see Doc 8047.

#include <vector>
#include <ostream>

#include "GeomPrimitives/inc/Cone.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSEnclosureMaker;

  class PSEnclosure : virtual public Detector {

  public:

    const int  version()    const { return version_; }
    const Tube& shell()     const { return shell_; }
    const Cone& shellCone() const { return shellCone_; }

    const Tube& endPlate()  const { return endPlate_; }
    const std::vector<Tube>& windows() const { return windows_; }

    unsigned nWindows()     const { return windows_.size(); }

    void setExtraOffset( double zOff ) { zOffset_ = zOff; }
    const double& getExtraOffset() const { return zOffset_; }

  protected:
    void setVersion( const int& aVers ) { version_ = aVers; }
  private:

    friend class PSEnclosureMaker;

    // Private ctr: the class should only be constructed via PSEnclosure::PSEnclosureMaker.
    PSEnclosure(const Tube& shell, const Tube& ep)
      : shell_(shell), shellCone_(), endPlate_(ep), version_(1)
    {};
    PSEnclosure (const Cone& shellCone, const Tube& ep )
      : shell_(), shellCone_(shellCone), endPlate_(ep), version_(2)
    {};

    // Or read back from persistent storage
    PSEnclosure();
    template<class T> friend class art::Wrapper;

    // The real enclosure shape is shown in docdb-2066 It is
    // approximated here by a cylinder closed with a flat end plate.
    Tube shell_;
    Cone shellCone_;
    Tube endPlate_;
    int version_;
    std::vector<Tube> windows_;

    // The updated version is shown in docdb-4087.  We use a conical frustrum.
    double zOffset_;

    // The version - allow cylindrical or conical PSEnclosure
  };

  std::ostream& operator<<(std::ostream& os, const PSEnclosure& pse);

}

#endif/*ProductionSolenoidGeom_PSEnclosure_hh*/
