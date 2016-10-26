#ifndef ServicesGeom_ElectronicRack_hh
#define ServicesGeom_ElectronicRack_hh

//
//
// Original author David Norvil Brown
// University of Louisville, Mu2e Collaboration
// March 2016 


#include <vector>
#include <ostream>

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "canvas/Persistency/Common/Wrapper.h"

namespace mu2e {

  class ElectronicRackMaker;

  class ElectronicRack : virtual public Detector {
  public:

    // Electronic racks can be built entirely from boxes of concrete
    // This C++ code only assumes the box shape.  Dimensions, locations,
    // and materials are specified in geometry text files.

    const std::vector<std::vector<double> >& getBoxDimensions() const { return _elecRackBoxDims; }
    const std::vector<std::string>& getMaterialNames() const { return _erMaterialNames; }
    const std::vector<CLHEP::Hep3Vector>& getCentersOfBoxes() const { return _erCenterPositions; }
    const std::vector<std::string>& getOrientations() const { return _erOrientations; }

  private:

    friend class ElectronicRackMaker;

    // Private ctr: the class should only be constructed via ElectronicRack::ElectronicRackMaker.
    ElectronicRack(const std::vector<std::vector<double> >& dims, 
		      const std::vector<std::string>& mats, 
		      const std::vector<CLHEP::Hep3Vector>& sites, 
		      const std::vector<std::string>& orients )
      : _elecRackBoxDims(dims),
	_erMaterialNames(mats),
	_erCenterPositions(sites),
	_erOrientations(orients)
    { }

    // Or read back from persistent storage
    ElectronicRack();
    template<class T> friend class art::Wrapper;


    // Current description based on Geometry 14, adapted by
    // David Norvil Brown, 

    std::vector< std::vector< double > > _elecRackBoxDims;
    std::vector< std::string >           _erMaterialNames;
    std::vector< CLHEP::Hep3Vector >     _erCenterPositions;
    std::vector< std::string >           _erOrientations;

  };

  std::ostream& operator<<(std::ostream& os, const ElectronicRack& upstr);

}

#endif/*ServicesGeom_ElectronicRack_hh*/
