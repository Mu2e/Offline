#ifndef RecoDataProducts_CrvCoincidenceCluster_hh
#define RecoDataProducts_CrvCoincidenceCluster_hh
//
//
// Contact person Ralf Ehrlich
//

#include "Offline/DataProducts/inc/CRVId.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include "CLHEP/Vector/ThreeVector.h"

#include <array>
#include <vector>

namespace mu2e
{
  class CrvCoincidenceCluster
  {
    public:

    CrvCoincidenceCluster() {}

    CrvCoincidenceCluster(int crvSectorType, double startTime, double endTime, float PEs,
              const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses, float slope, const std::vector<int> &layers,
              const std::array<size_t,CRVId::nSidesPerBar> &sideHits, const std::array<float, CRVId::nSidesPerBar> &sidePEs,
              const std::array<double,CRVId::nSidesPerBar> &sideTimes, bool twoReadoutSides, double avgHitTime, const CLHEP::Hep3Vector &avgHitPos) :
              _crvSectorType(crvSectorType), _startTime(startTime), _endTime(endTime), _PEs(PEs),
              _crvRecoPulses(crvRecoPulses), _slope(slope), _layers(layers),
              _sideHits(sideHits), _sidePEs(sidePEs), _sideTimes(sideTimes),
              _twoReadoutSides(twoReadoutSides), _avgHitTime(avgHitTime), _avgHitPos(avgHitPos) {}

    int                                         GetCrvSectorType() const {return _crvSectorType;}
    double                                      GetStartTime() const     {return _startTime;}
    double                                      GetEndTime() const       {return _endTime;}
    float                                       GetPEs() const           {return _PEs;}
    const std::vector<art::Ptr<CrvRecoPulse> > &GetCrvRecoPulses() const {return _crvRecoPulses;}
    float                                       GetSlope() const         {return _slope;}
    const std::vector<int>                     &GetLayers() const        {return _layers;}

    const std::array<size_t, CRVId::nSidesPerBar>   &GetSideHits() const        {return _sideHits;}
    const std::array<float, CRVId::nSidesPerBar>    &GetSidePEs() const         {return _sidePEs;}
    const std::array<double, CRVId::nSidesPerBar>   &GetSideTimes() const       {return _sideTimes;}
    bool                                             HasTwoReadoutSides() const {return _twoReadoutSides;}
    double                                           GetAvgHitTime() const      {return _avgHitTime;}
    const CLHEP::Hep3Vector                         &GetAvgHitPos() const       {return _avgHitPos;}

    // allow setting the Ptrs for reco compression
    void SetCrvRecoPulses(std::vector<art::Ptr<CrvRecoPulse> > const& pulses) { _crvRecoPulses = pulses; }

    private:

    int                                  _crvSectorType{0};
    double                               _startTime{0};
    double                               _endTime{0};
    float                                _PEs{0};
    std::vector<art::Ptr<CrvRecoPulse> > _crvRecoPulses;
    float                                _slope{0};
    std::vector<int>                     _layers;

    std::array<size_t, CRVId::nSidesPerBar>  _sideHits{0};  //number of hits on both readout sides of the modules of this cluster
    std::array<float, CRVId::nSidesPerBar>   _sidePEs{0};   //number of PEs on both readout sides
    std::array<double, CRVId::nSidesPerBar>  _sideTimes{0}; //average pulse times on both readout sides //entries are only valid, if the corresponding entries in sidePEs > 0;
    bool                                     _twoReadoutSides{false};  //indicates, if the average hit time and hit position are based on both readout sides
                                                                       //if not, the center of counters is used for the longitudinal position
    double                                   _avgHitTime{0};       //average hit time (based on the times of both readout sides, if available)
    CLHEP::Hep3Vector                        _avgHitPos;           //average hit position (based on the times of both readout sides, if available)
  };
  typedef std::vector<mu2e::CrvCoincidenceCluster> CrvCoincidenceClusterCollection;
}

#endif /* RecoDataProducts_CrvCoincidenceCluster_hh */
