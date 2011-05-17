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
// $Id: CrudeStrawHit.hh,v 1.8 2011/05/17 15:36:01 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:36:01 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <iosfwd>
#include <vector>
#include <string>

// Mu2e includes
#include "TrackerGeom/inc/StrawIndex.hh"
#include "ToyDP/inc/DPIndex.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/resolveDPIndices.hh"

// Forward declarations
namespace art{
  class Event;
}

namespace mu2e { 

  struct CrudeStrawHit{

  public:

    enum precursor_type { undefined, unpackedDigi, stepPointMC, saltAndPepper};

    // Data members for both data and MC events.
    precursor_type precursorType; // See note 1.
    StrawIndex     strawIndex;    // See note 2.
    float          driftDistance; // (mm)
    float          driftTime;     // (ns)
    float          sigmaD;        // (mm)
    float          energy;        // ( TBD: keV, MIP ?)

    std::vector<DPIndex> precursorIndices;  // See note 1.

    // Data members for objects created in MC processing.
    float  trueDriftDistance; // (mm)

    // Root requires us to have a default c'tor.
    // 2 suggestions on how to make it not available to others:
    //   - Use the BOOST_STATIC_ASSERT to make this visible only to GenReflex.
    //   - CrudeStrawHit( TRootIOCtor * ) .... see email from Philippe Nov 6/09.
    CrudeStrawHit():
      precursorType(undefined),
      strawIndex(StrawIndex(-1)),
      driftDistance(0.),
      driftTime(0.),
      sigmaD(0.),
      energy(0.),
      precursorIndices(),
      trueDriftDistance(0.),
      stepPointMCPointersValid(false),
      stepPointMCPointers(){
    }

    // Constructor for a hit that came from an unpacked digi, either 
    // from data or from the full MC chain.
    CrudeStrawHit( StrawIndex         strawIndex_,
                   float              driftDistance_,
                   float              driftTime_,
                   float              sigmaD_,
                   float              energy_,
                   DPIndex const&     precursorIndex_,
                   art::Event const * event_ = 0
                   );


    // Constructor from MC.
    CrudeStrawHit( StrawIndex                  strawIndex_,
                   float                       driftDistance_,
                   float                       driftTime_,
                   float                       sigmaD_,
                   float                       energy_,
                   precursor_type              precursorType_,
                   std::vector<DPIndex> const& precursorIndices_,
                   float                       trueDriftDistance_,
                   art::Event const *          event_ = 0
                   );

    // A special case of the previous c'tor when there is only one precursor.
    CrudeStrawHit( StrawIndex        strawIndex_,
                   float             driftDistance_,
                   float             driftTime_,
                   float             sigmaD_,
                   float             energy_,
                   precursor_type    precursorType_,
                   DPIndex const&    precursorIndex_,
                   float             trueDriftDistance_,
                   art::Event const* event_ = 0
                   );
    
    // Accept compiler generated versions of:
    //   d'tor
    //   copy c'tor 
    //   assignment operator

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    // Return the pointers to the precursors of this hit.
    std::vector<StepPointMC const *> const& getStepPointMC( art::Event const& event) const{
      resolveTransients(event);
      return stepPointMCPointers;
    }

    // Fill a std::vector with (pointers to const) of the precursors of this hit.
    void getStepPointMC( art::Event const&    event, 
                         std::vector<StepPointMC const*>& v ) const{
      resolveTransients(event);
      v.insert(v.end(), stepPointMCPointers.begin(), stepPointMCPointers.end());
    }

    // Return the pointers to the precursors of this hit.
    // This will throw if the pointers are not valid.
    // An empty collection can be a valid return.
    // The override is for debugging.
    std::vector<StepPointMC const *> const& getStepPointMCs( bool override = false) const;

    // Use this to check validity if there is a chance getStepPointMC might throw.
    bool stepPointMCsValid() const { return stepPointMCPointersValid;}

    // Resolve all of the transient information in this object.
    void resolveTransients( art::Event const& event) const;

  private:

    // Private data is not persistable and needs to be recomputed from the
    // persistent data after read in.
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
