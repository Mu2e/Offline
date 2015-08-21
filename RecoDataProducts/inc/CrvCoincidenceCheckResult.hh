#ifndef RecoDataProducts_CrvCoincidenceCheckResult_hh
#define RecoDataProducts_CrvCoincidenceCheckResult_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

#include <vector>

namespace mu2e 
{
  class CrvCoincidenceCheckResult
  {
    public:

    struct CoincidenceCombination
    {
      double                        _time[3];
      int                           _PEs[3];
      mu2e::CRSScintillatorBarIndex _counters[3];
      int                           _SiPMs[3];
      CoincidenceCombination() {}
    };

    struct CoincidenceHit
    {
      double                        _time;
      int                           _PEs;
      mu2e::CRSScintillatorBarIndex _counter;
      int                           _SiPM;
      CoincidenceHit(double time, int PEs, mu2e::CRSScintillatorBarIndex counter, int SiPM): _time(time), _PEs(PEs), _counter(counter), _SiPM(SiPM) {}
      bool operator<(const CoincidenceHit& rhs) const  //time ordered, 
      {                                                //identical hits will be removed, 
        if(_time!=rhs._time) return _time < rhs._time; //but allow for hits with the same time but on different counters/SiPMs
        if(_counter!=rhs._counter) return _counter < rhs._counter;
        return _SiPM < rhs._SiPM;
      }
    };

    struct DeadTimeWindow
    {
      double _startTime, _endTime;
      std::vector<CoincidenceHit> _hits;
      DeadTimeWindow(double start, double end, const std::vector<CoincidenceHit> &hits) : _startTime(start), _endTime(end), _hits(hits) {}
    };

    CrvCoincidenceCheckResult() {}

    const bool CoincidenceFound() const 
    {
      return (!_coincidenceCombinations.empty());
    }

    const std::vector<CoincidenceCombination> &GetCoincidenceCombinations() const
    {
      return _coincidenceCombinations;
    }

    std::vector<CoincidenceCombination> &GetCoincidenceCombinations()
    {
      return _coincidenceCombinations;
    }

    std::vector<DeadTimeWindow> GetDeadTimeWindows(double leadingTime, double trailingTime) const;

    private:

    std::vector<CoincidenceCombination> _coincidenceCombinations;
  };
}

#endif /* RecoDataProducts_CrvCoincidenceCheckResult_hh */
