#ifndef StrawMCHit_h
#define StrawMCHit_h 1
//
// A persistable G4 Hit - that is the intersection
// of a simulated track with a straw.
//
// $Id: StrawMCHit.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 

  class StrawMCHit{
  public:
  
    StrawMCHit():
      _trackId(-1),
      _strawIndex(-1),
      _edep(0.),
      _time(0.),
      _position(),
      _momentum(){
    }
    
    StrawMCHit( int trackId,
		int strawIndex,
		double edep,
		double time,
		CLHEP::Hep3Vector const& position,
		CLHEP::Hep3Vector const& momentum
		):
      _trackId(trackId),
      _strawIndex(strawIndex),
      _edep(edep),
      _time(time),
      _position(position),
      _momentum(momentum){
    }
    
    ~StrawMCHit(){
    }
    
    StrawMCHit(const StrawMCHit&);
    const StrawMCHit& operator=(const StrawMCHit&);
    
    void print( std::ostream& ost ) const;

    void print() const { print(std::cout); }

  public:
  
    int trackId()     const { return _trackId;    }
    int strawIndex()  const { return _strawIndex; }
    double eDep()     const { return _edep;       }      
    double time()     const { return _time;       }
    CLHEP::Hep3Vector position() const { return _position;   }
    CLHEP::Hep3Vector momentum() const { return _momentum;   }

  private:
  
    int               _trackId;
    int               _strawIndex;
    double            _edep;
    CLHEP::Hep3Vector _position;
    CLHEP::Hep3Vector _momentum;
    double            _time;
    
  };

} // namespace mu2e

#endif
