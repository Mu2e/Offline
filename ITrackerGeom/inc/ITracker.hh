#ifndef ITRACKER_HH
#define ITRACKER_HH

#include <deque>
#include <vector>

#include <boost/shared_array.hpp>

#include "ITrackerGeom/inc/SuperLayer.hh"
#include "GeometryService/inc/Detector.hh"

#include "ITrackerGeom/inc/CellGeometryHandle.hh"

namespace mu2e {

  class ITracker: public Detector{

  friend class ITrackerMaker;

  public:
    ITracker() {}
    ~ITracker() {}

    virtual std::string name() const { return "ITracker";}


    enum GeomType { Hexagonal=2, Square };
    
    double r0()   const { return _r0;}
    double z0()   const { return _z0;}
    double rOut() const { return _rOut;}
    int nSWire()  const { return _nSWire;}
    int nSDeltaWire() const { return _nSDeltaWire;}
    int nRing() const { return _nRing;}

    std::string extFile() const { return _extFile; }
    bool isExternal()     const { return _isExternal; }

    int nSuperLayer()  const { return _nSuperLayer; }

    double zHalfLength() const { return _zHalfLength;}

    double maxEndCapDim() const { return _max_EndCap_dim; }

    GeomType geomType() const { return _geomType; }

    CellGeometryHandle* getCellGeometryHandle() const { return _cellhnd.get(); }

    SuperLayer* getSuperLayer(int n) const throw(cms::Exception) {
    	if (n>=0 && n< _nSuperLayer){
    		return &(_sprlr[n]);
    	}
    	else throw cms::Exception("GEOM")<< "Super Layer number: "<< n <<" not present";
    }

    boost::shared_array<SuperLayer> getSuperLayersArray() const {
    	return _sprlr;
    }

protected:

    // Nominal values.  
    // _r0 = Nominal radius of the center of the sector.
    // _z0 = position of the center of the tracker relative to the origin
    //       of the Mu2e coordinate system.
    double _r0;
    double _z0;
    int _nSWire;
    int _nSDeltaWire;
    int _nRing;

    // Outer radius of a logical volume that will just contain the entire tracker.
    double _rOut;

    // Name of external gdml geometry file description.
    std::string _extFile;
    bool _isExternal;

    int _nSuperLayer;

    double _zHalfLength;
    double _max_EndCap_dim;

    //Cell geometry type: 2:Hexagonal, 3:Square
    GeomType _geomType;

    boost::shared_array<SuperLayer> _sprlr;
    
    std::auto_ptr<CellGeometryHandle> _cellhnd;

  };

} //namespace mu2e

#endif /*ITRACKER_HH*/
