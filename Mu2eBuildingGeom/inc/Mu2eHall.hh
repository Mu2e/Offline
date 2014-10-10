// Geometry of the hall, dirt, etc.
//
// Andrei Gaponenko, 2011

#ifndef MU2EHALL_HH
#define MU2EHALL_HH

#include <map>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class ExtrudedSolid {
  public:
    ExtrudedSolid() : yHalfThickness_(0.) {}

    ExtrudedSolid( const std::string& name,
		   const std::string& mat,
		   CLHEP::Hep3Vector& offset,
		   double yht,
		   std::vector<CLHEP::Hep2Vector>&& v ) 
      : name_(name)
      , material_(mat)
      , offsetFromMu2eOrigin_(offset)
      , yHalfThickness_(yht)
      , vertices_(v) 
    {}
    
    const std::string& getName()     const { return name_; }
    const std::string& getMaterial() const { return material_; }

    const CLHEP::Hep3Vector& getOffsetFromMu2eOrigin() const { return offsetFromMu2eOrigin_; }

    double getYhalfThickness() const { return yHalfThickness_; }

    const std::vector<CLHEP::Hep2Vector>& getVertices() const { return vertices_; }

    CLHEP::Hep2Vector& modifyVertex( int i ) { return vertices_.at( (std::size_t)i); }

  private:
    std::string name_;
    std::string material_;
    CLHEP::Hep3Vector offsetFromMu2eOrigin_;
    double yHalfThickness_;
    std::vector<CLHEP::Hep2Vector> vertices_;
  };


  //================================================================================
  class Mu2eHallMaker;

  class Mu2eHall : virtual public Detector {
  public:

    const std::map<std::string,ExtrudedSolid>& getBldgSolids() const { return bldgSolids_; }
    const std::map<std::string,ExtrudedSolid>& getDirtSolids() const { return dirtSolids_; }

    const ExtrudedSolid& 
    getBldgSolid( const std::string& str ) const {
      std::cout << (bldgSolids_.find( str ) == bldgSolids_.end() )<< std::endl;
      return bldgSolids_.find( str )->second;
    }

    const ExtrudedSolid& 
    getDirtSolid( const std::string& str ) const {
      return dirtSolids_.find( str )->second;
    }
 

    //----------------------------------------------------------------
  private:
    friend class Mu2eHallMaker;

    // Private ctr: the class should be only obtained via the maker
    Mu2eHall(){}

    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    std::map<std::string,ExtrudedSolid> bldgSolids_;
    std::map<std::string,ExtrudedSolid> dirtSolids_;

  };

}

#endif/*MU2EHALL_HH*/
