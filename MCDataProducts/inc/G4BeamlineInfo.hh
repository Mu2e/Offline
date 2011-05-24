#ifndef MCDataProducts_G4BeamlineInfo_hh
#define MCDataProducts_G4BeamlineInfo_hh
//
// Simple holder for extra data available from G4Beamline generator
//
#include <iostream>

namespace mu2e {

  class G4BeamlineInfo {

  public:

    G4BeamlineInfo(): _event_id(0), _track_id(0), _weight(0.0), _time(0.0) { }

    // Constructor for a hit that came from an unpacked digi, either
    // from data or from the full MC chain.
    G4BeamlineInfo( int event_id, int track_id, float weight, float time ) :
      _event_id(event_id), _track_id(track_id), _weight(weight), _time(time) { }

    // Accessors
    int   eventId() const { return _event_id; }
    int   trackId() const { return _track_id; }
    float weight()  const { return _weight; }
    float time()    const { return _time; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

    // Need this to create wrapper
    void swap(G4BeamlineInfo &);

  private:

    int   _event_id; // Event id from G4BL data file
    int   _track_id; // Track id from G4BL data file
    float _weight;   // Weight   from G4BL data file
    float _time;     // Time     from G4BL data file

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   G4BeamlineInfo const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* MCDataProducts_G4BeamlineInfo_hh */
