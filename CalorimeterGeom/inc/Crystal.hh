//
// Hold information about a crystal
// Neighbors and position are given in the "Mu2e" frame, i.e the global frame
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

             Crystal(int localId, int diskId, const CLHEP::Hep3Vector& localPosition, const CLHEP::Hep3Vector& size) : 
	        localId_(localId), 
	        diskId_(diskId), 
	        localPosition_(localPosition), 
	        position_(),
                size_(size), 
	        neighborsGlobal_(), 
	        nextNeighborsGlobal_(),
                neighborsGlobalRaw_(),
                nextNeighborsGlobalRaw_()
	     {}


             int localId()                                          const {return localId_;}
             int diskId()                                           const {return diskId_;}
             const CLHEP::Hep3Vector& localPosition()               const {return localPosition_;}
             const CLHEP::Hep3Vector& size()                        const {return size_;}
             
             //pre-computed global position and global neignbors to speed up code
             const CLHEP::Hep3Vector& position()                    const        {return position_;}
             const std::vector<int>&  neighbors(bool raw=false)     const        {return raw ? neighborsGlobalRaw_    : neighborsGlobal_;}	     
             const std::vector<int>&  nextNeighbors(bool raw=false) const        {return raw ? nextNeighborsGlobalRaw_: nextNeighborsGlobal_;}	     
             
             void setPosition(const CLHEP::Hep3Vector& pos)                      {position_ = pos;}
             void setNeighbors(const std::vector<int>& list, bool raw=false)     {raw ? neighborsGlobalRaw_=list     : neighborsGlobal_=list;}
             void setNextNeighbors(const std::vector<int>& list, bool raw=false) {raw ? nextNeighborsGlobalRaw_=list : nextNeighborsGlobal_=list;}    
             void adjustSize(const CLHEP::Hep3Vector& size)                      {size_ = size;}


	 private:
	 
	     int               localId_;
	     int               diskId_;
             CLHEP::Hep3Vector localPosition_;
             CLHEP::Hep3Vector position_;
             CLHEP::Hep3Vector size_;
             std::vector<int>  neighborsGlobal_;
             std::vector<int>  nextNeighborsGlobal_;
	     std::vector<int>  neighborsGlobalRaw_;
	     std::vector<int>  nextNeighborsGlobalRaw_;


     };

}


#endif 
