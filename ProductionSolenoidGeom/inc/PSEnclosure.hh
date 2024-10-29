#ifndef ProductionSolenoidGeom_PSEnclosure_hh
#define ProductionSolenoidGeom_PSEnclosure_hh

//
//
// Original author Andrei Gaponenko
//
// Dec 2017 - Dave (Louisville) Brown adds version 2, which makes the
// PSEnclosure a conical frustrum instead of a tube.  This is per as-built
// geometry - see Doc 8047.

// Sept  2021 - Michael MacKenzie and Dave (WKU - formerly Lou) Brown
// add version 3.  This adds more detail and gets closer to the
// representation in doc 8047.  Final code represents a merge of their
// independent work.

#include <vector>
#include <ostream>

#include "Offline/GeomPrimitives/inc/Polycone.hh"
#include "Offline/GeomPrimitives/inc/Cone.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PSEnclosureMaker;

  class PSEnclosure : virtual public Detector {

  public:

    const int  version()    const { return version_; }
    const Tube& shell()     const { return shell_; }
    const Cone& shellCone() const { return shellCone_; }

    void setFlange( Tube const& aFlange ) { flange_ = aFlange; }
    const Tube& flange()    const { return flange_; }

    const Tube& endPlate()  const { return endPlate_; }
    const std::vector<Tube>& windows() const { return windows_; }

    unsigned nWindows()     const { return windows_.size(); }

    const std::vector<bool>& hasWindowFrames() const
    { return hasFrames_; }
    const std::vector<bool>& hasWindowFramesOut() const
    { return hasFramesOut_; }
    const std::vector<Tube>& wFramesIn() const { return  wFramesIn_; }
    const std::vector<Tube>& wFramesOut() const { return  wFramesOut_; }

    const Polycone& endPlatePolycone()  const { return endPlatePolycone_; }
    const std::vector<Tube>& windowPipes() const { return windowPipes_; }

    void setExtraOffset( double zOff ) { zOffset_ = zOff; }
    const double& getExtraOffset() const { return zOffset_; }

    const double& getEndPlateInset() const { return ePinset_; }
    const double& getEndPlateRadius() const { return ePradius_; }

  protected:
    void setVersion( const int& aVers ) { version_ = aVers; }
  public:

    friend class PSEnclosureMaker;

    // Private ctr: the class should only be constructed via PSEnclosure::PSEnclosureMaker.
    PSEnclosure(const Tube& shell, const Tube& ep)
      : shell_(shell), shellCone_(), endPlate_(ep), version_(1)
    {};
    PSEnclosure (const Cone& shellCone, const Tube& ep, int vers=2 )
      : shell_(), shellCone_(shellCone), endPlate_(ep), version_(vers)
    {};
    PSEnclosure (const Cone& shellCone, const Polycone& ep, int vers=3 )
      : shell_(), shellCone_(shellCone), endPlatePolycone_(ep), version_(vers)
    {};
  private:
    // Or read back from persistent storage
    PSEnclosure();
    template<class T> friend class art::Wrapper;

    // The real enclosure shape is shown in docdb-2066 It is
    // approximated here by a cylinder closed with a flat end plate.
    Tube shell_;
    Cone shellCone_;
    Tube endPlate_;
    Polycone endPlatePolycone_;
    int version_;
    std::vector<Tube> windows_;
    std::vector<Tube> windowPipes_;
    double ePinset_;
    double ePradius_;

    // For version 3, include flanges (or "frames") for each window

    std::vector<bool> hasFrames_;
    std::vector<bool> hasFramesOut_;
    std::vector<Tube> wFramesIn_;
    std::vector<Tube> wFramesOut_;

    // The updated version is shown in docdb-4087.  We use a conical frustrum.
    double zOffset_;
    // In version 3 we add a flange and...
    Tube flange_;


  };

  std::ostream& operator<<(std::ostream& os, const PSEnclosure& pse);

}

#endif/*ProductionSolenoidGeom_PSEnclosure_hh*/
