#ifndef ExternalShieldingGeom_ExtShieldUpstream_hh
#define ExternalShieldingGeom_ExtShieldUpstream_hh

//
//
// Original author David Norvil Brown
// University of Louisville, Mu2e Collaboration
// October 2014 - Geometry #14 version implemented.

// The Upstream external shielding is described in WBS 5.4.  The
// geometry is described here for the TS Upstream portion of the shielding.
// The PS external shield is described for the Offline in the
// ProductionSolenoidGeom package.  Details can be found in docdb #4999
// Found in this class:
// TSu Shield
// West Wall Shield
// Upstream cap
// poly shield (a very ad hoc version in the initial implementation)
// Protection Collimator

#include <vector>
#include <ostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class ExtShieldUpstreamMaker;

  class ExtShieldUpstream : virtual public Detector {
    // DNB:  Detector is a very lightweight base class, so no harm
    // inheriting from it.  Would it make sense to have a similar base
    // class for shields, though?  Sort of a fake identifier.
  public:

    // The upstream shielding can be built entirely from boxes of concrete
    // This C++ code only assumes the box shape.  Dimensions, locations,
    // and materials are specified in geometry text files.

    const std::vector<std::vector<double> >& getBoxDimensions() const { return _extShieldBoxDims; }
    const std::vector<std::vector<double> >& getBoxTolerances() const { return _extShieldBoxTols; }
    const std::vector<std::string>& getMaterialNames() const { return _materialNames; }
    const std::vector<CLHEP::Hep3Vector>& getCentersOfBoxes() const { return _centerPositions; }
    const std::vector<std::string>& getOrientations() const { return _orientations; }

  private:

    friend class ExtShieldUpstreamMaker;

    // Private ctr: the class should only be constructed via ExtShieldUpstream::ExtShieldUpstreamMaker.
    ExtShieldUpstream(const std::vector<std::vector<double> >& dims, 
		      const std::vector<std::vector<double> >& tols, 
		      const std::vector<std::string>& mats, 
		      const std::vector<CLHEP::Hep3Vector>& sites, 
		      const std::vector<std::string>& orients )
      : _extShieldBoxDims(dims),
	_extShieldBoxTols(tols),
	_materialNames(mats),
	_centerPositions(sites),
	_orientations(orients)
    { }

    // Or read back from persistent storage
    ExtShieldUpstream();
    template<class T> friend class art::Wrapper;


    // Current description based on Geometry 14, adapted by
    // David Norvil Brown, 

    std::vector< std::vector< double > > _extShieldBoxDims;
    std::vector< std::vector< double > > _extShieldBoxTols;
    std::vector< std::string >           _materialNames;
    std::vector< CLHEP::Hep3Vector >     _centerPositions;
    std::vector< std::string >           _orientations;

  };

  std::ostream& operator<<(std::ostream& os, const ExtShieldUpstream& upstr);

}

#endif/*ExternalShieldingGeom_ExtShieldUpstream_hh*/
