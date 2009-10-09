#ifndef CrudeStrawHit_h
#define CrudeStrawHit_h 1
//
// A persistable class representing a crude hit on a straw.
// Crude means that represents a hit straight off of the detector
// with calibration applied.  It is not yet assciated with any
// cluster, or with any track, so it cannot have full calibration
// applied.
//
// Notes:
// 1) About precurors:
//    a) For data from the experiment, this hit will arise from the
//       processing of exactly 1 unpacked digitized waveform.  An index to that
//       waveform is stored in precursor.
//
//    b) For MC events that have been through the full MC chain, including
//       creation of digis, the precursor will also index to exactly 1 unpacked
//       digitized waveform.
//
//    c) For MC events for which there was no simulation of digis, the precursor
//       will be -1 and the alternate precursor, altPrecursor, will point back to 
//       one or more StepPointMC objects.
//
//    d) One can also imagine middle fidelty MC for which the precursor will be -1
//       and altPrecursor will point back to 1 or more objects that are 
//       intermediate between StepPointMC and unpacked digitized waveforms.
//
// 2) For the LTracker this is a dense index 0...(N-1).  For the others it is
//    to be defined.
// 
// 
// $Id: CrudeStrawHit.hh,v 1.1 2009/10/09 13:31:32 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/10/09 13:31:32 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iostream>
#include <vector>
#include <string>

// Mu2e includes
#include "LTrackerGeom/inc/StrawIndex.hh"

namespace mu2e { 

  struct CrudeStrawHit{

  public:

    // Root requires us to have a default c'tor.
    CrudeStrawHit();

    // Constructor for a hit either from data or from the
    // full MC chain.
    CrudeStrawHit( int        precursor_,
		   StrawIndex strawIdx_,
		   float      driftDistance_,
		   float      driftTime_,
		   float      sigmaD_,
		   float      energy_
		   );


    // Constructor from an MC without digis.
    CrudeStrawHit( StrawIndex               strawIdx_,
		   float                    driftDistance_,
		   float                    driftTime_,
		   float                    sigmaD_,
		   float                    energy_,
		   std::vector<int> const & altPrecursor_,
		   float                    trueDriftDistance_
		   );

    // A special case of the previous c'tor when there is only one 
    // altPrecursor.
    CrudeStrawHit( StrawIndex  strawIdx_,
		   float       driftDistance_,
		   float       driftTime_,
		   float       sigmaD_,
		   float       energy_,
                   int         altPrecursor_,
		   float       trueDriftDistance_
		   );
    
    // Accept compiler generated versions of:
    //   d'tor
    //   copy c'tor 
    //   assignment operator

    // Formatted information as a string as a printed to a stream.
    std::string toString() const;

    void print( std::ostream& ost ) const;

    inline void print() const { 
      print(std::cout); 
    }

    // Data members:
    int        precursor;     // See note 1.
    StrawIndex strawIdx;      // See note 2.
    float      driftDistance; // (mm)
    float      driftTime;     // (ns)
    float      sigmaD;        // (mm)
    float      energy;        // ( TBD: keV, MIP ?)

    // Data members for objects created in MC processing.
    std::vector<int>  altPrecursor;      // See note 1.
    float             trueDriftDistance; // (mm)

  private:
    std::string formatAltPrecursors() const;
    
  };

  inline std::ostream& operator<<( std::ostream& ost,
				   CrudeStrawHit const& hit){
    ost << hit.toString();
    return ost;
  }
  

} // namespace mu2e

#endif
