#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"

#include <set>

namespace mu2e 
{

  std::vector<CrvCoincidenceCheckResult::DeadTimeWindow> CrvCoincidenceCheckResult::GetDeadTimeWindows(double leadingTime, double trailingTime) const
  {
    std::vector<DeadTimeWindow> deadTimeWindows;

    //create time ordered set of single coincidence hits and remove redundant hits
    std::set<CoincidenceHit> hits;

    std::vector<CoincidenceCombination>::const_iterator iter;
    for(iter=_coincidenceCombinations.begin(); iter!=_coincidenceCombinations.end(); iter++)
    {
      for(int i=0; i<3; i++) hits.emplace(iter->_time[i], iter->_PEs[i], iter->_counters[i], iter->_SiPMs[i]);
    }

    std::set<CoincidenceHit>::const_iterator h=hits.begin();
    while(h!=hits.end())
    {
      double startTime = h->_time - leadingTime;
      double endTime = h->_time + trailingTime;
      std::vector<CoincidenceHit> hitsToStore;
      hitsToStore.push_back(*h);
      while(++h != hits.end())
      {
        if(endTime >= h->_time - leadingTime) 
        {
          endTime = h->_time + trailingTime;
          hitsToStore.push_back(*h);
        }
        else break;
      }
      deadTimeWindows.emplace_back(startTime,endTime,hitsToStore);
    }

    return(deadTimeWindows);
  }

}
