#ifndef DataProducts_CaloHit_hh
#define DataProducts_CaloHit_hh

#include <iostream>
#include <vector>

namespace mu2e {

  class CaloHit
  {
     public:

       CaloHit():
          _roId(-1), _time(0.), _energyDep(0.) 
       {}

       CaloHit( int roId, float time, float energyDep ):
         _roId(roId), _time(time),_energyDep(energyDep)
       {}

       int        id()         const { return _roId; }
       float      time()       const { return _time;}
       float      energyDep()  const { return _energyDep; }

       void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

     private:
        int       _roId;
        float     _time;     
        float     _energyDep;
  };

  
  inline std::ostream& operator<<(std::ostream& ost, CaloHit const& hit)
  {
      hit.print(ost,false);
      return ost;
  }
  
  
  typedef std::vector<mu2e::CaloHit> CaloHitCollection;
} 

#endif /* DataProducts_CaloHit_hh */
