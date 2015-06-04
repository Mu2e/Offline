#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"

#include <set>

namespace mu2e 
{

  std::vector<CrvCoincidenceCheckResult::DeadTimeWindow> CrvCoincidenceCheckResult::GetDeadTimeWindows(double leadingTime, double trailingTime) const
  {
    std::vector<CrvCoincidenceCheckResult::DeadTimeWindow> deadTimeWindows;
    std::set<double> times;
    std::vector<CoincidenceCombination>::const_iterator iter;
    for(iter=_coincidenceCombinations.begin(); iter!=_coincidenceCombinations.end(); iter++)
    {
      times.insert(iter->_time, iter->_time+3);
    }

    std::set<double>::const_iterator t=times.begin();
    while(t!=times.end())
    {
      double startTime = *t - leadingTime;
      double endTime = *t + trailingTime;
      while(++t != times.end())
      {
        if(endTime > *t - leadingTime) endTime = *t + trailingTime;
        else break;
      }
      deadTimeWindows.emplace_back(startTime,endTime);
    }

    return(deadTimeWindows);
  }

}
