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

    CrvCoincidenceCluster(int crvSectorType, const CLHEP::Hep3Vector &avgCounterPos, double startTime, double endTime, float PEs,
              const std::vector<art::Ptr<CrvRecoPulse> > &crvRecoPulses, float slope, const std::vector<int> &layers,
              const std::array<size_t,CRVId::nSidesPerBar> &sideHits, const std::array<float, CRVId::nSidesPerBar> &sidePEs,
              const std::array<double,CRVId::nSidesPerBar> &sideTimes) :
              _crvSectorType(crvSectorType), _avgCounterPos(avgCounterPos), _startTime(startTime), _endTime(endTime), _PEs(PEs),
              _crvRecoPulses(crvRecoPulses), _slope(slope), _layers(layers),
              _sideHits(sideHits), _sidePEs(sidePEs), _sideTimes(sideTimes) {}

    int                                         GetCrvSectorType() const {return _crvSectorType;}
    const CLHEP::Hep3Vector                    &GetAvgCounterPos() const {return _avgCounterPos;}
    double                                      GetStartTime() const     {return _startTime;}
    double                                      GetEndTime() const       {return _endTime;}
    float                                       GetPEs() const           {return _PEs;}
    const std::vector<art::Ptr<CrvRecoPulse> > &GetCrvRecoPulses() const {return _crvRecoPulses;}
    float                                       GetSlope() const         {return _slope;}
    const std::vector<int>                     &GetLayers() const        {return _layers;}

    const std::array<size_t, CRVId::nSidesPerBar>   &GetSideHits() const    {return _sideHits;}
    const std::array<float, CRVId::nSidesPerBar>    &GetSidePEs() const     {return _sidePEs;}
    const std::array<double, CRVId::nSidesPerBar>   &GetSideTimes() const   {return _sideTimes;}

    // allow setting the Ptrs for reco compression
    void SetCrvRecoPulses(std::vector<art::Ptr<CrvRecoPulse> > const& pulses) { _crvRecoPulses = pulses; }

    private:

    int                                  _crvSectorType;
    CLHEP::Hep3Vector                    _avgCounterPos;
    double                               _startTime;
    double                               _endTime;
    float                                _PEs;
    std::vector<art::Ptr<CrvRecoPulse> > _crvRecoPulses;
    float                                _slope;
    std::vector<int>                     _layers;

    std::array<size_t, CRVId::nSidesPerBar>  _sideHits;  //number of hits on both sides of the modules of this cluster
    std::array<float, CRVId::nSidesPerBar>   _sidePEs;   //number of PEs on both sides
    std::array<double, CRVId::nSidesPerBar>  _sideTimes; //average pulse times on both sides //entries are only valid, if the corresponding entries in sidePEs > 0;
  };
  typedef std::vector<mu2e::CrvCoincidenceCluster> CrvCoincidenceClusterCollection;
}

#endif /* RecoDataProducts_CrvCoincidenceCluster_hh */
