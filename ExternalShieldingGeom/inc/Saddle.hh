#ifndef ExternalShieldingGeom_Saddle_hh
#define ExternalShieldingGeom_Saddle_hh

//
//
// Original author David Norvil Brown
// University of Louisville, Mu2e Collaboration
// March 2015 - Geometry #14 version implemented.

// Saddles for the various cryostats hold them in place.
// The saddles can be described using extrusions.  Copy
// the code from ExtShieldDownstream.
// Details can be found in docdb #4678 


#include <vector>
#include <ostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class SaddleMaker;

  class Saddle : virtual public Detector {
  public:

    // 

    int                                      getVersion() const
    { return _version; }

    const std::vector<std::vector<std::vector<double> > >& getOutlines() const 
    { return _saddleOutlines; }
    const std::vector<double> &              getLengths() const 
    { return _saddleLengths; }
    const std::vector<std::vector<double> >& getTolerances() const 
    { return _saddleBoxTols; }
    const std::vector<std::string>&          getMaterialNames() const 
    { return _materialNames; }
    const std::vector<CLHEP::Hep3Vector>&    getCentersOfBoxes() const 
    { return _centerPositions; }
    const std::vector<std::string>&          getOrientations() const 
    { return _orientations; }
    const std::vector<int>&                  getNHoles() const 
    { return _nHoles; }
    const std::vector<int>&                  getNNotches() const
    { return _nNotches; }
    const std::vector<int>&                  getHoleIndices() const
    { return _holeIndices; }
    const std::vector<CLHEP::Hep3Vector>&    getHoleLocations() const
    { return _holeLocations; }
    const std::vector<double>&               getHoleRadii() const
    { return _holeRadii; }
    const std::vector<double>&               getHoleLengths() const
    { return _holeLengths; }
    const std::vector<std::string>&          getHoleOrientations() const
    { return _holeOrientations; }
    const std::vector<int>&                  getNotchIndices() const
    { return _notchIndices; }
    const std::vector<CLHEP::Hep3Vector>&    getNotchLocations() const
    { return _notchLocations; }
    const std::vector<std::vector<double> >& getNotchDimensions() const
    { return _notchDimensions; }

  private:

    friend class SaddleMaker;

    // Private ctr: the class should only be constructed via Saddle::SaddleMaker.
    Saddle(const int                                version,
	   const std::vector<std::vector<std::vector<double> > >& outlines,
	   const std::vector<double>&               lengths,
	   const std::vector<std::vector<double> >& tols, 
	   const std::vector<std::string>&          mats, 
	   const std::vector<CLHEP::Hep3Vector>&    sites, 
	   const std::vector<std::string>&          orients,
	   const std::vector<int>&                  nHoles,
	   const std::vector<int>&                  nNotches,
	   const std::vector<int>&                  iHole,
	   const std::vector<CLHEP::Hep3Vector>&    locHole,
	   const std::vector<double>&               radHole,
	   const std::vector<double>&               lenHole,
	   const std::vector<std::string>&          oHole,
	   const std::vector<int>&                  iNotch,
	   const std::vector<CLHEP::Hep3Vector>&    locNotch,
	   const std::vector<std::vector<double> >& locDims)
      : _version(version),
	_saddleOutlines(outlines),
	_saddleLengths (lengths),
	_saddleBoxTols (tols),
	_materialNames    (mats),
	_centerPositions  (sites),
	_orientations     (orients),
	_nHoles           (nHoles),
	_nNotches         (nNotches),
	_holeIndices      (iHole),
	_holeLocations    (locHole),
        _holeRadii        (radHole),
	_holeLengths      (lenHole),
	_holeOrientations (oHole),
	_notchIndices     (iNotch),
	_notchLocations   (locNotch),
	_notchDimensions  (locDims)
    { }

    // Or read back from persistent storage
    Saddle();
    template<class T> friend class art::Wrapper;


    // Current description based on Geometry 14, adapted by
    // David Norvil Brown, 
    int _version;
    // The following vectors hold one piece of information per "block"
    std::vector<std::vector< std::vector< double > > > _saddleOutlines;
    std::vector< double >                _saddleLengths;
    std::vector< std::vector< double > > _saddleBoxTols;
    std::vector< std::string >           _materialNames;
    std::vector< CLHEP::Hep3Vector >     _centerPositions;
    std::vector< std::string >           _orientations;
    std::vector< int >                   _nHoles;
    std::vector< int >                   _nNotches;
    std::vector< int >                   _holeIndices;
    // The following vectors hold one piece of information per "hole"
    std::vector< CLHEP::Hep3Vector >     _holeLocations;
    std::vector< double >                _holeRadii;
    std::vector< double >                _holeLengths;
    std::vector< std::string >           _holeOrientations;
    // The following vectors hold one piece of information per "notch"
    std::vector< int >                   _notchIndices;
    std::vector< CLHEP::Hep3Vector >     _notchLocations;
    std::vector< std::vector< double > > _notchDimensions;
  };

  std::ostream& operator<<(std::ostream& os, const Saddle& upstr);

}

#endif/*ExternalShieldingGeom_Saddle_hh*/
