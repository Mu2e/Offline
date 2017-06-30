#ifndef ExternalShieldingGeom_Pipe_hh
#define ExternalShieldingGeom_Pipe_hh

//
//
// Original author David Norvil Brown
// University of Louisville, Mu2e Collaboration
// March 2016 

// Pipes for various services.  Each pipe is really a set of pipes 
// (referred to here as "components") contained
// within the first.  If only one component - it is just a plain single pipe.
// The pipes can be described using G4Tubs.  
// Details can be found in docdb #4678 


#include <vector>
#include <ostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class PipeMaker;

  class Pipe : virtual public Detector {
  public:

    // Straight sections of pipe and pipe bends can be (separately) 
    // represented with the same number of parameters, as long as we
    // assume bends will always be of 90 degrees.  If needed later,
    // can add a "bendAngle" parameter which will be meaningless for a 
    // straight section.

    // ***
    // These represet type-level information
    // ***

    const int &                           getVersion() const
    { return _version; }
    const std::vector<int> &              getNComponentsInPipe() const
    { return _nComponentsInPipe; }
    const std::vector<int> &              getNPipes() const
    { return _nPipes; }
    // If the flavor (see below) is "bend", then the length is really
    // the toroidal radius of the containing pipe
    // In version 2, length is not a property of a type, but of an individual
    // pipe.  Use the same vector for code simplicity.
    const std::vector<double> &           getLengths() const
    { return _lengths; }
    const std::vector<std::string> &      getFlavor() const
    { return _flavors; }
    const std::vector<std::string> &      getFillMaterialNames() const
    { return _fillMaterialNames; }

    // ***
    // The following represent pipe-level (top-level) information
    // ***

    // If a straight section, center means the actual center of the pipe.
    // If a bend, then center means the center of curvature of the toroidal
    // section
    const std::vector<std::vector<CLHEP::Hep3Vector> >&    getCentersOfPipes() const 
    { return _centerPositions; }
    const std::vector<std::vector<std::string> >&          getOrientations() const 
    { return _orientations; }

    // The following are component-level information
    const std::vector<std::vector<double> >&  getInnerRads() const
    { return _innerRads; }
    const std::vector<std::vector<double> >&  getOuterRads() const
    { return _outerRads; }
    const std::vector<std::vector<std::string> >& getMaterialNames() const
    { return _materialNames; }
    const std::vector<std::vector<double> >&  getUOffsets() const
    { return _uOffsets; }
    const std::vector<std::vector<double> >&  getVOffsets() const
    { return _vOffsets; }

  private:

    friend class PipeMaker;

    // Private ctr: the class should only be constructed via Pipe::PipeMaker.
    Pipe(const int&                               version,
	 const std::vector<int>&                  nComponentsInPipe,
	 const std::vector<int>&                  nPipes,
	 const std::vector<double>&               lengths,
	 const std::vector<std::string>&          flavs,
	 const std::vector<std::string>&          fillMats,
	 const std::vector<std::vector<CLHEP::Hep3Vector> >&    sites, 
	 const std::vector<std::vector<std::string> >&          orients,
	 const std::vector<std::vector<double> >& innerRads, 
	 const std::vector<std::vector<double> >& outerRads, 
	 const std::vector<std::vector<std::string> >&          mats, 
	 const std::vector<std::vector<double> >& uOffsets, 
	 const std::vector<std::vector<double> >& vOffsets)      
      : _version          (version), 
	_nComponentsInPipe (nComponentsInPipe),
	_nPipes           (nPipes),
	_lengths          (lengths),
	_flavors          (flavs),
	_fillMaterialNames     (fillMats),
	_centerPositions  (sites),
	_orientations     (orients),
	_innerRads        (innerRads),
	_outerRads        (outerRads),
	_materialNames    (mats),
	_uOffsets         (uOffsets),
	_vOffsets         (vOffsets)
    { }

    // Or read back from persistent storage
    Pipe();
    template<class T> friend class art::Wrapper;


    // Current description based on Geometry 14, adapted by
    // David Norvil Brown, 

    int                                  _version;
    // The following vectors hold information about the pipes
    // The first five are type-level information
    std::vector< int >                   _nComponentsInPipe;
    std::vector< int >                   _nPipes;
    std::vector< double >                _lengths;
    std::vector< std::string >           _flavors;
    std::vector< std::string >           _fillMaterialNames;
    // The next two are pipe-level information
    std::vector< std::vector<CLHEP::Hep3Vector > >     _centerPositions;
    std::vector< std::vector<std::string > >           _orientations;
    // The last five are component-level information
    std::vector< std::vector< double > > _innerRads;
    std::vector< std::vector< double > > _outerRads;
    std::vector< std::vector<std::string > >    _materialNames;
    std::vector< std::vector< double > > _uOffsets;
    std::vector< std::vector< double > > _vOffsets;
  };

  std::ostream& operator<<(std::ostream& os, const Pipe& upstr);

}

#endif/*ExternalShieldingGeom_Pipe_hh*/
