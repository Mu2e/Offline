#include <iostream>
#include <map>
#include <string>

#include "CRVAnalysis/inc/CRVAnalysisInfo.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "art/Framework/Principal/Event.h"

namespace mu2e
{
  class CRVAnalysis
  {
    public:
    CRVAnalysis(const std::string &g4ModuleLabel, const std::string &crvCoincidenceModuleLabel);
    void FillCRVInfoMCStructure(const art::Event& event, CRVAnalysisInfo &info);
   
    private:
    struct CoincidenceHit
    {
      double _time;
      int    _PEs;
      mu2e::CRSScintillatorBarIndex _counter;
      CoincidenceHit(double time, int PEs, mu2e::CRSScintillatorBarIndex counter) : _time(time), _PEs(PEs), _counter(counter) {}
      bool operator==(const CoincidenceHit &hit) const
      {
        return(_time==hit._time && _PEs==hit._PEs && _counter==hit._counter);
      }
    };
    typedef std::vector<CoincidenceHit> CoincidenceCluster;

    struct CrvPlane
    {
      int _coordinate;
      double _position;
      CrvPlane() : _coordinate(-1), _position(NAN) {} //needed by map [] access operator
      CrvPlane(int coordinate, double position) : _coordinate(coordinate), _position(position) {}
    };
    std::map<int,CrvPlane> _crvPlanes;

    CLHEP::Hep3Vector _detSysOrigin;

    bool FindCrvPlaneCrossings(const art::Event& event, const cet::map_vector_key& particleKey, const int &sectorType,
                               CLHEP::Hep3Vector &point, CLHEP::Hep3Vector &direction);

    std::string _g4ModuleLabel;
    std::string _crvCoincidenceModuleLabel;
  };

}


