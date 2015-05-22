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

    struct DeadTimeWindow
    {
      double _startTime, _endTime;
      DeadTimeWindow(double start, double end) : _startTime(start), _endTime(end) {}
    };

    CrvCoincidenceCheckResult() {}

    void SetCoincidence(bool coincidence)
    {
      _coincidence = coincidence;
    }

    const bool GetCoincidence() const 
    {
      return _coincidence;
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

    bool                                _coincidence;
    std::vector<CoincidenceCombination> _coincidenceCombinations;
  };
}

#endif /* RecoDataProducts_CrvCoincidenceCheckResult_hh */
