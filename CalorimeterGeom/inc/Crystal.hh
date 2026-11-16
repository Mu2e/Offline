#ifndef CalorimeterGeom_Crystal_hh
#define CalorimeterGeom_Crystal_hh
//
// Hold information about a crystal
// ID, Neighbors, and position are given in the "Mu2e" global reference frame
// localID and localPosition are given in the local disk refernce frame
//
// Original author B Echenard
//
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>

namespace mu2e {

  class Crystal {
    public:
      static constexpr unsigned invalidID_ = CaloConst::_invalid;

      Crystal(unsigned localID, unsigned diskID,
              const CLHEP::Hep3Vector& localPosition,
              const CLHEP::Hep3Vector& size) :
         ID_(invalidID_),
         localID_(localID),
         diskID_(diskID),
         size_(size),
         localPosition_(localPosition),
         position_(),
         neighbors_(),
         nextNeighbors_()
      {}

      unsigned                       ID           () const {return ID_;}
      unsigned                       localID      () const {return localID_;}
      unsigned                       diskID       () const {return diskID_;}
      const CLHEP::Hep3Vector&       size         () const {return size_;}
      const CLHEP::Hep3Vector&       localPosition() const {return localPosition_;}
      const CLHEP::Hep3Vector&       position     () const {return position_;}
      const std::vector<unsigned>&   neighbors    () const {return neighbors_;}
      const std::vector<unsigned>&   nextNeighbors() const {return nextNeighbors_;}

      void setID           (unsigned ID)                       {ID_ = ID;}
      void setLocalID      (unsigned localID)                  {localID_ = localID;}
      void setDiskID       (unsigned diskID)                   {diskID_ = diskID;}
      void setLocalPosition(const CLHEP::Hep3Vector& pos)      {localPosition_ = pos;}
      void setPosition     (const CLHEP::Hep3Vector& pos)      {position_ = pos;}
      void setNeighbors    (const std::vector<unsigned>& list) {neighbors_ = list;}
      void setNextNeighbors(const std::vector<unsigned>& list) {nextNeighbors_ = list;}

    private:
      unsigned              ID_;
      unsigned              localID_;
      unsigned              diskID_;
      CLHEP::Hep3Vector     size_;
      CLHEP::Hep3Vector     localPosition_;
      CLHEP::Hep3Vector     position_;
      std::vector<unsigned> neighbors_;
      std::vector<unsigned> nextNeighbors_;
  };
}
#endif
