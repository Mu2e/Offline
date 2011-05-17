//
// A crudely calibrated hit in a straw. See header for full details.
//
// $Id: CrudeStrawHit.cc,v 1.8 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author Rob Kutschke

// C++ includes
#include <ostream>

// Framework includes.
#include "art/Framework/Core/Event.h"
#include "cetlib/exception.h"

// Mu2e includes
#include "ToyDP/inc/CrudeStrawHit.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"

using namespace std;

namespace mu2e {

  CrudeStrawHit::CrudeStrawHit( StrawIndex        strawIndex_,
                                float             driftDistance_,
                                float             driftTime_,
                                float             sigmaD_,
                                float             energy_,
                                DPIndex const&    precursorIndex_,
                                art::Event const* event_
                                ):
    precursorType(unpackedDigi),
    strawIndex(strawIndex_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    precursorIndices(1,precursorIndex_),
    trueDriftDistance(0.),
    stepPointMCPointersValid(false),
    stepPointMCPointers()
  {
    if ( event_){
      resolveTransients( *event_ );
    }
  }

  CrudeStrawHit::CrudeStrawHit( StrawIndex                  strawIndex_,
                                float                       driftDistance_,
                                float                       driftTime_,
                                float                       sigmaD_,
                                float                       energy_,
                                precursor_type              precursorType_,
                                std::vector<DPIndex> const& precursorIndices_,
                                float                       trueDriftDistance_,
                                art::Event const*           event_
                                ):
    precursorType(precursorType_),
    strawIndex(strawIndex_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    precursorIndices(precursorIndices_),
    trueDriftDistance(trueDriftDistance_),
    stepPointMCPointersValid(false),
    stepPointMCPointers()
  {
    if ( event_){
      resolveTransients( *event_ );
    }
  }

  CrudeStrawHit::CrudeStrawHit( StrawIndex        strawIndex_,
                                float             driftDistance_,
                                float             driftTime_,
                                float             sigmaD_,
                                float             energy_,
                                precursor_type    precursorType_,
                                DPIndex  const&   precursorIndex_,
                                float             trueDriftDistance_,
                                art::Event const* event_
                                ):
    precursorType(precursorType_),
    strawIndex(strawIndex_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    precursorIndices(1,precursorIndex_),
    trueDriftDistance(trueDriftDistance_),
    stepPointMCPointersValid(false),
    stepPointMCPointers()
  {
    if ( event_){
      resolveTransients( *event_ );
    }

  }

  // Return the pointers to the precursors of this hit.
  // Throw if they are not available or if we have overridden the safety check.
  std::vector<StepPointMC const *> const& CrudeStrawHit::getStepPointMCs( bool override ) const{
    if ( stepPointMCPointersValid || override ) {
      return stepPointMCPointers;
    }
    throw cet::exception("ProductNotFound")
      << "Cannot compute pointers to StepPointMC without an event being supplied.";
  }  

  // Populate the transient data members.
  // This will get more complicated as different sorts of precursors become available.
  void CrudeStrawHit::resolveTransients( art::Event const& event) const{

    if (stepPointMCPointersValid ) return;

    if ( precursorType != stepPointMC ) {
      throw cet::exception("ProductNotFound")
        << "Cannot compute pointers to StepPointMC from a precursor of type: "
        << precursorType;
    }

    resolveDPIndices<StepPointMCCollection>( event, precursorIndices, stepPointMCPointers);
    stepPointMCPointersValid = true;
  }


  void CrudeStrawHit::print( ostream& ost, bool doEndl ) const {

    ost << "CrudeStraw Hit:"
        << " pretyp: "    << precursorType
        << " id: "        << strawIndex
        << " d: "         << driftDistance
        << " t: "         << driftTime
        << " s: "         << sigmaD
        << " e: "         << energy
        << " pre: ";

    formatPrecursorIndices(ost);

    ost << " tru: "       << trueDriftDistance;

    if ( doEndl ){
      ost << endl;
    }
    
  }
  
  void CrudeStrawHit::formatPrecursorIndices( ostream& ost) const {

    if ( precursorIndices.empty() ) return;
    
    if ( precursorIndices.size() == 1 ){
      ost << precursorIndices[0];
    } else{
      ost << "(";
      for ( vector<int>::size_type i=0;
            i<precursorIndices.size(); ++i ){
        ost << precursorIndices[i];
        if ( i < precursorIndices.size()-1 ){
          ost << ",";
        }
      }
      ost << ")";
    }
    
  }

} // namespace mu2e
