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
              const std::array<double,CRVId::nSidesPerBar> &sideTimes, double avgHitTime, const CLHEP::Hep3Vector &avgHitPos,
               bool hitPosAndTimeCalculated) :
              _crvSectorType(crvSectorType), _startTime(startTime), _endTime(endTime), _PEs(PEs),
              _crvRecoPulses(crvRecoPulses), _slope(slope), _layers(layers),
              _sideHits(sideHits), _sidePEs(sidePEs), _sideTimes(sideTimes),
              _avgHitTime(avgHitTime), _avgHitPos(avgHitPos),
              _hitPosAndTimeCalculated(hitPosAndTimeCalculated) {}

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
    double                                           GetAvgHitTime() const      {return _avgHitTime;}
    const CLHEP::Hep3Vector                         &GetAvgHitPos() const       {return _avgHitPos;}

    bool HitPosAndTimeCalculated() const  {return _hitPosAndTimeCalculated;}

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
    std::array<double, CRVId::nSidesPerBar>  _sideTimes{0}; //average pulse times on both readout sides
                                                            //entries are only valid, if the corresponding entries in sidePEs > 0;
    double                                   _avgHitTime{0};//-hitPosAndTimeCalculated true (correlated hits at both readout sides):
                                                            //likely time when the particle hit the counter at the avgHitPos
                                                            //based on time difference between both sides and the fiber signal speed.
                                                            //-otherwise:
                                                            //likely time when the particle hit the counter under the assumption that the
                                                            //particle hit the counter at the avgHitPos (PE-weighted average center of all
                                                            //involved counters) and based on the signal travel time from the SiPM to
                                                            //the avgHitPos. if hits on both sides: use the average
                                                            //-an additional average time due to electronics response and
                                                            //processes such as scintillation/WLS decay times is subtracted.
    CLHEP::Hep3Vector                        _avgHitPos;    //-hitPosAndTimeCalculated true (correlated hits at both readout sides):
                                                            //likely position where the particle hit the counter
                                                            //based on time difference between both sides and the fiber signal speed.
                                                            //-otherwise:
                                                            //PE-weighted average center of all involved counters.
    bool _hitPosAndTimeCalculated{false};         //indicates that the hits at both readout sides seem to be correlated, i.e. they could have been
                                                  //caused by the same particle (based on the time difference, signal speed, and distance)
  };
  typedef std::vector<mu2e::CrvCoincidenceCluster> CrvCoincidenceClusterCollection;
}

#endif /* RecoDataProducts_CrvCoincidenceCluster_hh */
