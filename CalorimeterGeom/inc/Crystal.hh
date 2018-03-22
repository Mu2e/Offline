//
// Hold information about a crystal
// Neighbors and position are given in the "Mu2e" frame
// localId and localPosition are given in the disk frame (i.e. local frame)
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

             Crystal(int localId, int diskId, const CLHEP::Hep3Vector &localPosition, 
                     const CLHEP::Hep3Vector &localPositionFrontFace) : 
	        localId_(localId), 
	        diskId_(diskId), 
	        localPosition_(localPosition), 
	        localPositionFrontFace_(localPositionFrontFace), 
	        position_(), 
	        neighbors_(), 
	        nextNeighbors_(),
                neighborsRaw_(),
                nextNeighborsRaw_()
	     {}


             int   localId()                                        const {return localId_;}
             int   diskId()                                         const {return diskId_;}
             const CLHEP::Hep3Vector& localPosition()               const {return localPosition_;}
             const CLHEP::Hep3Vector& localPositionFF()             const {return localPositionFrontFace_;}
             const CLHEP::Hep3Vector& position()                    const  {return position_;}
             const std::vector<int>&  neighbors(bool raw=false)     const {return raw ? neighborsRaw_ :neighbors_;}	     
             const std::vector<int>&  nextNeighbors(bool raw=false) const {return raw ? nextNeighborsRaw_: nextNeighbors_;}	     
           
             
             void setPosition(const CLHEP::Hep3Vector& pos)                      {position_ = pos;}
             void setNeighbors(const std::vector<int>& list, bool raw=false)     {raw ? neighborsRaw_=list : neighbors_=list;}
             void setNextNeighbors(const std::vector<int>& list, bool raw=false) {raw ? nextNeighborsRaw_=list : nextNeighbors_=list;}    


	 private:
	 
	     int                 localId_;
	     int                 diskId_;
             CLHEP::Hep3Vector   localPosition_;
             CLHEP::Hep3Vector   localPositionFrontFace_;
             CLHEP::Hep3Vector   position_;
             std::vector<int>    neighbors_;
             std::vector<int>    nextNeighbors_;
	     std::vector<int>    neighborsRaw_;
	     std::vector<int>    nextNeighborsRaw_;


     };

}


#endif 
