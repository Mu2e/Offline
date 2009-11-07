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
//       processing of exactly 1 unpacked digitized waveform.
//
//    b) For MC events that have been through the full MC chain, including
//       creation of digis, this hit will arise from the processing of 
//       exactly 1 unpacked digitized waveform.
//
//    c) For MC events for which there was no simulation of digis, this
//       event can arise from several different sources.  One example
//       is that it could come directly from StepPointMC objects. In
//       this case there hit could come from more than 1 precursor objects.
//
//    To represent the above information there are two data members.
//       precursor_type precursorType;
//       std::vector<DPIndex> precursorIndices;
//
//    The first tells us something about how these hits were made;
//    see the enum precursor_type for possible values.
//    The second is vector of DPIndex objects.  A DPIndex specifies
//    a data product by its product ID and an index into that data product.
//
//
// 2) For the LTracker this is a dense index 0...(N-1).  For the others it is
//    to be defined.
// 
// 
// $Id: CrudeStrawHit.hh,v 1.4 2009/11/07 01:09:25 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/07 01:09:25 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iosfwd>
#include <vector>
#include <string>

// Mu2e includes
#include "LTrackerGeom/inc/StrawIndex.hh"
#include "ToyDP/inc/DPIndex.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"

// Forward declarations
namespace edm{
  class Event;
}

namespace mu2e { 

  struct CrudeStrawHit{

  public:

    enum precursor_type { undefined, unpackedDigi, stepPointMC};

    // Data members:
    precursor_type precursorType; // See note 1.
    StrawIndex     strawIndex;    // See note 2.
    float          driftDistance; // (mm)
    float          driftTime;     // (ns)
    float          sigmaD;        // (mm)
    float          energy;        // ( TBD: keV, MIP ?)

    // Data members for objects created in MC processing.
    std::vector<DPIndex> precursorIndices;  // See note 1.
    float                trueDriftDistance; // (mm)

    // Root requires us to have a default c'tor.
    // Use the BOOST_STATIC_ASSERT to make this visible only to GenReflex.
    CrudeStrawHit():
      precursorType(undefined),
      strawIndex(StrawIndex(-1)),
      driftDistance(0.),
      driftTime(0.),
      sigmaD(0.),
      energy(0.),
      precursorIndices(),
      trueDriftDistance(0.){
    }


    // Constructor for a hit that came from an unpacked digi, either 
    // from data or from the full MC chain.
    CrudeStrawHit( StrawIndex     strawIndex_,
		   float          driftDistance_,
		   float          driftTime_,
		   float          sigmaD_,
		   float          energy_,
		   DPIndex const& precursorIndex_
		   );


    // Constructor from MC.
    CrudeStrawHit( StrawIndex                  strawIndex_,
		   float                       driftDistance_,
		   float                       driftTime_,
		   float                       sigmaD_,
		   float                       energy_,
		   precursor_type              precursorType_,
		   std::vector<DPIndex> const& precursorIndices_,
		   float                       trueDriftDistance_
		   );

    // A special case of the previous c'tor when there is only one  altPrecursor.
    CrudeStrawHit( StrawIndex       strawIndex_,
		   float            driftDistance_,
		   float            driftTime_,
		   float            sigmaD_,
		   float            energy_,
		   precursor_type   precursorType_,
		   DPIndex const&   precursorIndex_,
		   float            trueDriftDistance_
		   );
    
    // Accept compiler generated versions of:
    //   d'tor
    //   copy c'tor 
    //   assignment operator

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    // Return the pointers to the precursors of this hit.
    std::vector<StepPointMC const *> const& getStepPointMC( edm::Event const& event) const{
      resolveTransients(event);
      return stepPointMCPointers;
    }

    // Fill a std::vector with (pointers to const) of the precursors of this hit.
    void getStepPointMC( edm::Event const&    event, 
			 std::vector<StepPointMC const*>& v ) const{
      resolveTransients(event);
      v.insert(v.end(), stepPointMCPointers.begin(), stepPointMCPointers.end());
    }

    // Return the pointers to the precursors of this hit.
    // This will throw if the pointers are not valid.
    std::vector<StepPointMC const *> const& getStepPointMC() const;

    // Resolve all of the transient information in this object.
    void resolveTransients( edm::Event const& event) const;

    // Call this after readback to ensure that the object is in a safe state.
    // In a future release of root, this will not be necessary.
    void resetTransients() const{
      stepPointMCPointers.clear();
      stepPointMCPointersValid = false;
    }

    bool checkValid() const { return stepPointMCPointersValid;}
    std::vector<StepPointMC const*> const& checkPointers() const {
      return stepPointMCPointers;
    }


  private:

    // Private data is not persistable and needs to be recomputed from the
    // persistent data after read in.  Use the bool to do lazy evaluation when
    // possible.

    mutable bool stepPointMCPointersValid;
    mutable std::vector<StepPointMC const *> stepPointMCPointers;

    // Helper function to aid in printing.
    void formatPrecursorIndices( std::ostream& ost ) const;
    
  };

  inline std::ostream& operator<<( std::ostream& ost,
				   CrudeStrawHit const& hit){
    hit.print(ost,false);
    return ost;
  }
  

} // namespace mu2e

#endif
