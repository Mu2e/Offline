#ifndef MCDataProducts_MixingSummary_hh
#define MCDataProducts_MixingSummary_hh
//
// Status information from one call to an event mixing module.
//
// $Id: MixingSummary.hh,v 1.2 2011/10/12 20:04:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/10/12 20:04:08 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepInstanceName.hh"

// Includes from art and its tool chain.
#include "art/Framework/IO/ProductMix/MixTypes.h"

// C++ includes
#include <iosfwd>

namespace mu2e {

  class MixingSummary {

  public:

    typedef StepInstanceName::enum_type enum_type;

    // Identify the different StepPointMCCollections.
    //enum stepIndex { tracker, virtualDetector, timeVD, stoppingTarget, crv, calorimeter, caloReadout, stepIndexSize };

    MixingSummary():
      status_(),
      eventStatus_(),
      eventIDs_(),
      genSizes_(),
      simDeltas_(),
      stepSizes_(StepInstanceName::size(),std::vector<size_t>()),
      pointTrajectoryDeltas_(){
    }

    // Accept compiler written versions of d'tor, copy c'tor and assignment operator.

    // Accessors/Modifier pairs

    StatusG4   status() const { return status_; }
    StatusG4&  status() { return status_; }

    art::EventIDSequence const& eventIDs() const { return eventIDs_;}
    art::EventIDSequence&       eventIDs()       { return eventIDs_;}

    std::vector<size_t> const& genSizes() const { return genSizes_; }
    std::vector<size_t>&       genSizes()       { return genSizes_; }

    std::vector<size_t> const& simDeltas() const { return simDeltas_; }
    std::vector<size_t>&       simDeltas()       { return simDeltas_; }

    std::vector<size_t> const& stepSizes(enum_type i) const { return stepSizes_.at(i); }
    std::vector<size_t>&       stepSizes(enum_type i)       { return stepSizes_.at(i); }

    std::vector<size_t> const& stepSizes(std::string name) const { return stepSizes_.at(StepInstanceName::findByName(name).id()); }
    std::vector<size_t>&       stepSizes(std::string name)       { return stepSizes_.at(StepInstanceName::findByName(name).id()); }

    std::vector<size_t> const& trackerSizes()         const { return stepSizes_[StepInstanceName::tracker];}
    std::vector<size_t> const& virtualDetectorSizes() const { return stepSizes_[StepInstanceName::virtualdetector];}
    std::vector<size_t> const& timeVDSizes()          const { return stepSizes_[StepInstanceName::timeVD];}
    std::vector<size_t> const& stoppingTargetSizes()  const { return stepSizes_[StepInstanceName::stoppingtarget];}
    std::vector<size_t> const& crvSizes()             const { return stepSizes_[StepInstanceName::CRV];}
    std::vector<size_t> const& calorimeterSizes()     const { return stepSizes_[StepInstanceName::calorimeter];}
    std::vector<size_t> const& caloReadoutSizes()     const { return stepSizes_[StepInstanceName::calorimeterRO];}

    std::vector<size_t> const& pointTrajectoryDeltas() const { return pointTrajectoryDeltas_; }
    std::vector<size_t>&       pointTrajectoryDeltas()       { return pointTrajectoryDeltas_; }

    std::vector<StatusG4> const& eventStatus() const { return eventStatus_;}
    std::vector<StatusG4>&       eventStatus()       { return eventStatus_;}

    void swap( MixingSummary& rhs ){
      std::swap( status_,                rhs.status_                );
      std::swap( eventStatus_,           rhs.eventStatus_           );
      std::swap( eventIDs_,              rhs.eventIDs_              );
      std::swap( genSizes_,              rhs.genSizes_              );
      std::swap( simDeltas_,             rhs.simDeltas_             );
      std::swap( stepSizes_,             rhs.stepSizes_             );
      std::swap( pointTrajectoryDeltas_, rhs.pointTrajectoryDeltas_ );
    }

    void print ( std::ostream& ) const;

  private:

    // A rollup of all of the contributing StatusG4 objects.
    StatusG4 status_;

    // A copy of the StatusG4 object from each input event.
    std::vector<StatusG4> eventStatus_;

    // EventIDs from the input file.
    art::EventIDSequence eventIDs_;

    // For collections that are std::vector<T>, record the size of each input collection.
    // For collections that are cet::map_vector<T>, record the delta of each input collection;
    // the delta is one more that the last valid key in the cet::map_vector ( or 0 if the
    // collection is empty).
    std::vector<size_t>               genSizes_;
    std::vector<size_t>               simDeltas_;
    std::vector<std::vector<size_t> > stepSizes_;
    std::vector<size_t>               pointTrajectoryDeltas_;

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   MixingSummary const& summary){
    summary.print(ost);
    return ost;
  }


} // end namespace mu2e

#endif /* MCDataProducts_MixingSummary_hh */
