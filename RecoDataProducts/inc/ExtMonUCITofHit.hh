#ifndef RecoDataProducts_ExtMonUCITofHit_hh
#define RecoDataProducts_ExtMonUCITofHit_hh

//

// C++ includes
#include <iostream>

namespace art {
  class ProductID;
}

namespace mu2e {

  class ExtMonUCITofHit{

  public:

    ExtMonUCITofHit():
      _stationId(-1),
      _segmentId(-1),
      _time(0.),
      _energyDep(0.)
    {}

    ExtMonUCITofHit(int stationId, int segmentId, double time, double energy):
      _stationId(stationId),
      _segmentId(segmentId),
      _time(time),
      _energyDep(energy)
    {}

    // Accessors

    int              stationId() const { return _stationId; }
    int              segmentId() const { return _segmentId; }
    float            time()      const { return _time; }
    float            energyDep() const { return _energyDep; }

    // Accept compiler generated versions of d'tor, copy c'tor, assignment operator.

    void setEnergyDep(double energy);

    // Print contents of the object.
    void print( std::ostream& ost = std::cout, bool doEndl = true ) const;

   private:

    int              _stationId;
    int              _segmentId;
    float            _time;             // (ns)
    float            _energyDep;        // (MeV)

  };

  inline std::ostream& operator<<( std::ostream& ost,
                                   ExtMonUCITofHit const& hit){
    hit.print(ost,false);
    return ost;
  }

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonUCITofHit_hh */
