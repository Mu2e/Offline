//
// A crudely calibrated hit in a straw. See header for full details.
//
// $Id: CrudeStrawHit.cc,v 1.2 2009/10/22 16:31:08 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/22 16:31:08 $
//
// Original author Rob Kutschke

// C++ includes
#include <ostream>

// Framework includes.
#include "FWCore/Framework/interface/Event.h"

// Mu2e includes
#include "ToyDP/inc/CrudeStrawHit.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"

using namespace std;

namespace mu2e {


  CrudeStrawHit::CrudeStrawHit():
    precursorType(undefined),
    strawIdx(StrawIndex::fromInt(0)),
    driftDistance(0.),
    driftTime(0.),
    sigmaD(0.),
    energy(0.),
    precursorIndices(),
    trueDriftDistance(0.)
  {
  }


  CrudeStrawHit::CrudeStrawHit( StrawIndex     strawIdx_,
                                float          driftDistance_,
                                float          driftTime_,
                                float          sigmaD_,
                                float          energy_,
				DPIndex const& precursorIndex_
                                ):
    precursorType(unpackedDigi),
    strawIdx(strawIdx_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    precursorIndices(1,precursorIndex_),
    trueDriftDistance(0.)
  {
  }

  CrudeStrawHit::CrudeStrawHit( StrawIndex                  strawIdx_,
                                float                       driftDistance_,
                                float                       driftTime_,
                                float                       sigmaD_,
                                float                       energy_,
				precursor_type              precursorType_,
				std::vector<DPIndex> const& precursorIndices_,
                                float                       trueDriftDistance_
                                ):
    precursorType(precursorType_),
    strawIdx(strawIdx_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    precursorIndices(precursorIndices_),
    trueDriftDistance(trueDriftDistance_)
  {
  }

  CrudeStrawHit::CrudeStrawHit( StrawIndex       strawIdx_,
                                float            driftDistance_,
                                float            driftTime_,
                                float            sigmaD_,
                                float            energy_,
				precursor_type   precursorType_,
				DPIndex  const&  precursorIndex_,
                                float            trueDriftDistance_
                                ):
    precursorType(precursorType_),
    strawIdx(strawIdx_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    precursorIndices(1,precursorIndex_),
    trueDriftDistance(trueDriftDistance_)
  {
  }

  void CrudeStrawHit::print( ostream& ost, bool doEndl ) const {

    ost << "CrudeStraw Hit:"
	<< " pretyp: "    << precursorType
        << " id: "        << strawIdx
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
