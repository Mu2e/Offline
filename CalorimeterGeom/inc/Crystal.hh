//
// Hold information about a crystal
// Neighbors and position are given in the "Mu2e" frame, i.e the global frame
// localID and localPosition are given in the disk frame (i.e. local frame)
// Front face is sometimes abreviated FF (yes, I am lazy)
//
// Original author B Echenard
//

#ifndef CalorimeterGeom_Crystal_hh
#define CalorimeterGeom_Crystal_hh

#include "CLHEP/Vector/ThreeVector.h"
#include <vector>

namespace mu2e {

     class Crystal {

         public:
           Crystal(int localID, int diskID, const CLHEP::Hep3Vector& localPosition,
                   const CLHEP::Hep3Vector& idealLocalPosition, const CLHEP::Hep3Vector& size) :
              localID_               (localID),
              diskID_                (diskID),
              size_                  (size),
              localPosition_         (localPosition),
              idealLocalPosition_    (idealLocalPosition),
              position_              (),
              neighborsGlobal_       (),
              nextNeighborsGlobal_   (),
              neighborsGlobalRaw_    (),
              nextNeighborsGlobalRaw_()
           {}


           int                      localID      ()       const {return localID_;}
           int                      diskID       ()       const {return diskID_;}
           const CLHEP::Hep3Vector& size         ()       const {return size_;}
           const CLHEP::Hep3Vector& localPosition()       const {return localPosition_;}
           const CLHEP::Hep3Vector& idealLocalPosition()  const {return localPosition_;}
           const CLHEP::Hep3Vector& position     ()       const {return position_;}
           const std::vector<int>&  neighbors    ()       const {return neighborsGlobal_;}
           const std::vector<int>&  nextNeighbors()       const {return nextNeighborsGlobal_;}

           void setSize         (const CLHEP::Hep3Vector& size) {size_ = size;}
           void setLocalPosition(const CLHEP::Hep3Vector& pos)  {localPosition_ = pos;}
           void setPosition     (const CLHEP::Hep3Vector& pos)  {position_ = pos;}
           void setNeighbors    (const std::vector<int>& list)  {neighborsGlobal_=list;}
           void setNextNeighbors(const std::vector<int>& list)  {nextNeighborsGlobal_=list;}


         private:
           int               localID_;
           int               diskID_;
           CLHEP::Hep3Vector size_;
           CLHEP::Hep3Vector localPosition_;
           CLHEP::Hep3Vector idealLocalPosition_;
           CLHEP::Hep3Vector position_;
           std::vector<int>  neighborsGlobal_;
           std::vector<int>  nextNeighborsGlobal_;
           std::vector<int>  neighborsGlobalRaw_;
           std::vector<int>  nextNeighborsGlobalRaw_;
     };

}


#endif
