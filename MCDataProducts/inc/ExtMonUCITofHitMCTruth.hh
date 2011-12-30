// $Id: ExtMonUCITofHitMCTruth.hh,v 1.1 2011/12/30 20:31:46 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/30 20:31:46 $
//

#ifndef MCDataProducts_ExtMonUCITofHitMCTruth_hh
#define MCDataProducts_ExtMonUCITofHitMCTruth_hh

// C++ includes
#include <iostream>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

// Mu2e includes

namespace mu2e {

  struct ExtMonUCITofHitMCTruth{

  public:

    ExtMonUCITofHitMCTruth():
      _stationId(-1),
      _segmentId(-1),
      _time(0.),
      _energyDep(0.),
      _trackId(-1),
      _pdgId(-1),
      _position(),
      _momentum(),
      _vertex(),
      _vertexMomentum(),
      _vertexTime(0.0),
      _isPrimary(-1),
      _orgTrackId(-1),
      _orgPdgId(-1),
      _orgVertex(),
      _orgVertexMomentum(),
      _orgTime(0.0)
    {
    }

    ExtMonUCITofHitMCTruth( int   stationId,
                            int   segmentId,
                            float time,
                            float energyDep,
                            int   trackId,
                            int   pdgId,
                            CLHEP::Hep3Vector const&     position,
                            CLHEP::Hep3Vector const&     momentum,
                            CLHEP::Hep3Vector const&     vertex,
                            CLHEP::Hep3Vector const&     vertexMomentum,
                            float vertexTime,
                            int   isPrimary,
                            int   orgTrackId,
                            int   orgPdgId,
                            CLHEP::Hep3Vector const&     orgVertex,
                            CLHEP::Hep3Vector const&     orgVertexMomentum,
                            float orgTime ) :
      _stationId(stationId),
      _segmentId(segmentId),
      _time(time),
      _energyDep(energyDep),
      _trackId(trackId),
      _pdgId(pdgId),
      _position(position),
      _momentum(momentum),
      _vertex(vertex),
      _vertexMomentum(vertexMomentum),
      _vertexTime(vertexTime),
      _isPrimary(isPrimary),
      _orgTrackId(orgTrackId),
      _orgPdgId(orgPdgId),
      _orgVertex(orgVertex),
      _orgVertexMomentum(orgVertexMomentum),
      _orgTime(orgTime)
    {
    }

    // Accessors
    int   stationId() const { return _stationId; }
    int   segmentId() const { return _segmentId; }
    float time()      const { return _time;}
    float energyDep() const { return _energyDep; }
    int   trackId()   const { return _trackId; }
    int   pdgId()     const { return _pdgId; }
    CLHEP::Hep3Vector const&     position()         const { return _position;  }
    CLHEP::Hep3Vector const&     momentum()         const { return _momentum;  }
    CLHEP::Hep3Vector const&     vertex()              const { return _vertex;  }
    CLHEP::Hep3Vector const&     vertexMomentum()      const { return _vertexMomentum;  }
    float vertexTime() const { return _vertexTime; }
    int   isPrimary()  const { return _isPrimary; }
    int   orgTrackId() const { return _orgTrackId; }
    int   orgPdgId()   const { return _orgPdgId; }
    CLHEP::Hep3Vector const&     orgVertex()              const { return _orgVertex;  }
    CLHEP::Hep3Vector const&     orgVertexMomentum()      const { return _orgVertexMomentum;  }
    float orgTime()   const { return _orgTime; }    

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

  private:

    int   _stationId;
    int   _segmentId;
    float _time;             // (ns)
    float _energyDep;        // (MeV)
    int   _trackId;
    int   _pdgId;
    CLHEP::Hep3Vector     _position;
    CLHEP::Hep3Vector     _momentum;
    CLHEP::Hep3Vector     _vertex;
    CLHEP::Hep3Vector     _vertexMomentum;
    float _vertexTime;
    int   _isPrimary;
    int   _orgTrackId;
    int   _orgPdgId;
    CLHEP::Hep3Vector     _orgVertex; 
    CLHEP::Hep3Vector     _orgVertexMomentum;
    float _orgTime;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   ExtMonUCITofHitMCTruth const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* MCDataProducts_ExtMonUCITofHitMCTruth_hh */
