#ifndef MCDataProducts_MixingSummary_hh
#define MCDataProducts_MixingSummary_hh
//
// Status information from one call to an event mixing module.
//
// $Id: MixingSummary.hh,v 1.1 2011/06/30 04:39:13 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/30 04:39:13 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "MCDataProducts/inc/StatusG4.hh"

// Includes from art and its tool chain.
#include "art/Framework/IO/ProductMix/MixContainerTypes.h"

// C++ includes
//#include <iosfwd>
#include <iostream>
#include <vector>

namespace mu2e {

  class MixingSummary {

  public:

    // Identify the different StepPointMCCollections.
    enum stepIndex { tracker, virtualDetector, stoppingTarget, crv, calorimeter, caloReadout, stepIndexSize };

    MixingSummary():
      status_(),
      eventStatus_(),
      eventIDs_(),
      genOffsets_(),
      simOffsets_(),
      stepOffsets_(stepIndexSize,std::vector<size_t>()),
      pointTrajectoryOffsets_(){
      std::cerr << "Default constructed mixing summary: "
                << stepOffsets_.size()       << " "
                << stepOffsets_.at(0).size() << " "
                << stepOffsets_.at(5).size() << " "
                << std::endl;
    }

    // Accept compiler written versions of d'tor, copy c'tor and assignment operator.

    // Accessors/Modifier pairs

    StatusG4   status() const { return status_; }
    StatusG4&  status() { return status_; }

    art::EventIDSequence const& eventIDs() const { return eventIDs_;}
    art::EventIDSequence&       eventIDs()       { return eventIDs_;}

    std::vector<size_t> const& genOffsets() const { return genOffsets_; }
    std::vector<size_t>&       genOffsets()       { return genOffsets_; }

    std::vector<size_t> const& simOffsets() const { return simOffsets_; }
    std::vector<size_t>&       simOffsets()       { return simOffsets_; }

    std::vector<size_t> const& stepOffsets(stepIndex i) const { return stepOffsets_.at(i); }
    std::vector<size_t>&       stepOffsets(stepIndex i)       { return stepOffsets_.at(i); }
    std::vector<size_t> const& trackerOffsets()         const { return stepOffsets_[tracker];}
    std::vector<size_t> const& virtualDetectorOffsets() const { return stepOffsets_[virtualDetector];}
    std::vector<size_t> const& stoppingTargetOffsets()  const { return stepOffsets_[stoppingTarget];}
    std::vector<size_t> const& crvOffsets()             const { return stepOffsets_[crv];}
    std::vector<size_t> const& calorimeterOffsets()     const { return stepOffsets_[calorimeter];}
    std::vector<size_t> const& caloReadoutOffsets()     const { return stepOffsets_[caloReadout];}

    std::vector<size_t> const& pointTrajectoryOffsets( stepIndex ) const { return pointTrajectoryOffsets_; }
    std::vector<size_t>&       pointTrajectoryOffsets( stepIndex )       { return pointTrajectoryOffsets_; }

    std::vector<StatusG4> const& eventStatus() const { return eventStatus_;}

    void swap( MixingSummary& rhs ){
      std::swap( status_,                 rhs.status_                 );
      std::swap( eventStatus_,            rhs.eventStatus_            );
      std::swap( eventIDs_,               rhs.eventIDs_               );
      std::swap( genOffsets_,             rhs.genOffsets_             );
      std::swap( simOffsets_,             rhs.simOffsets_             );
      std::swap( stepOffsets_,            rhs.stepOffsets_            );
      std::swap( pointTrajectoryOffsets_, rhs.pointTrajectoryOffsets_ );
    }

    //void print ( std::ostream&, bool newLine=true ) const;

  private:

    // Overall status of the merge operation.
    StatusG4 status_;

    // A Copy of the StatusG4 object from each input event.
    std::vector<StatusG4> eventStatus_;

    // EventIDs from the input file.
    art::EventIDSequence eventIDs_;

    // Position within the output collection at which each input colleciton starts.
    std::vector<size_t> genOffsets_;
    std::vector<size_t> simOffsets_;

    // Key value within the output collection at which each input colleciton starts.
    std::vector<std::vector<size_t> >stepOffsets_;
    std::vector<size_t>              pointTrajectoryOffsets_;

  };

  /*
  inline std::ostream& operator<<( std::ostream& ost,
                                   MixingSummary const& stat){
    stat.print(ost,false);
    return ost;
  }
  */


} // end namespace mu2e

#endif /* MCDataProducts_MixingSummary_hh */
