//
// A crudely calibrated hit in a straw. See header for full details.
//
// $Id: CrudeStrawHit.cc,v 1.1 2009/10/09 13:31:32 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/09 13:31:32 $
//
// Original author Rob Kutschke

// C++ includes
#include <iostream>
#include <sstream>

// Mu2e incldues
#include "ToyDP/inc/CrudeStrawHit.hh"

using namespace std;

namespace mu2e {


  CrudeStrawHit::CrudeStrawHit():
    precursor(-1),
    strawIdx(-1),
    driftDistance(0.),
    driftTime(0.),
    sigmaD(0.),
    energy(0.),
    altPrecursor(),
    trueDriftDistance(0.)
  {
  }

  CrudeStrawHit::CrudeStrawHit( int         precursor_,
                                StrawIndex  strawIdx_,
                                float       driftDistance_,
                                float       driftTime_,
                                float       sigmaD_,
                                float       energy_
                                ):
    precursor(precursor_),
    strawIdx(strawIdx_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    altPrecursor(),
    trueDriftDistance(0.)
  {
  }

  CrudeStrawHit::CrudeStrawHit( StrawIndex               strawIdx_,
                                float                    driftDistance_,
                                float                    driftTime_,
                                float                    sigmaD_,
                                float                    energy_,
                                std::vector<int> const & altPrecursor_,
                                float                    trueDriftDistance_
                                ):
    precursor(-1),
    strawIdx(strawIdx_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    altPrecursor(altPrecursor_),
    trueDriftDistance(trueDriftDistance_)
  {
  }

  CrudeStrawHit::CrudeStrawHit( StrawIndex strawIdx_,
                                float      driftDistance_,
                                float      driftTime_,
                                float      sigmaD_,
                                float      energy_,
                                int        altPrecursor_,
                                float      trueDriftDistance_
                                ):
    precursor(-1),
    strawIdx(strawIdx_),
    driftDistance(driftDistance_),
    driftTime(driftTime_),
    sigmaD(sigmaD_),
    energy(energy_),
    altPrecursor(),
    trueDriftDistance(trueDriftDistance_)
  {
    altPrecursor.push_back(altPrecursor_);
  }

  void CrudeStrawHit::print( ostream& ost ) const {
    ost << toString() << endl;
  }

  string CrudeStrawHit::toString() const {

    ostringstream oss;

    oss << "CrudeStraw Hit:"
	<< " pre: "    << precursor
        << " id: "     << strawIdx
        << " d: "      << driftDistance
        << " t: "      << driftTime
        << " s: "      << sigmaD
        << " e: "      << energy
        << " alt: "    << formatAltPrecursors()
        << " tru: "    << trueDriftDistance;
    
    return oss.str();
    
  }
  
  string CrudeStrawHit::formatAltPrecursors() const {
    if ( altPrecursor.size() == 0 ) return string();
    
    ostringstream oss;
    if ( altPrecursor.size() == 1 ){
      oss << altPrecursor[0];
    } else{
      oss << "(";
      for ( vector<int>::size_type i=0;
            i<altPrecursor.size(); ++i ){
        oss << altPrecursor[i];
        if ( i < altPrecursor.size()-1 ){
          oss << ",";
        }
      }
      oss << ")";
    }

    return oss.str();
  }

} // namespace mu2e
