//forms clusters of adjacent straws in the L-Tracker for pattern recognition
//
// $Id: HitCluster.cc,v 1.5 2011/05/18 16:11:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 16:11:17 $
//
//original author R. Bernstein
//
#include "HitCluster/inc/HitCluster.hh"
using namespace std;
namespace mu2e{
  namespace hitcluster{

    HitCluster::~HitCluster(){};

    HitCluster::hitNeighbours HitCluster::findHitNeighbours()
    {
      vector<StrawIndex> const& nearest = _straw->nearestNeighboursByIndex();
      hitNeighbours nearbyStraws;

      //put candidate under study at front of vector
      int firstStraw = _hit->volumeId();
      Candidate _cand;
      _cand.id = firstStraw;
      _cand.hitPointer = _hit;
      nearbyStraws.push_back(_cand);

      //and now find everyone else
      for ( vector<int>::size_type ihit =0;
            ihit<nearest.size(); ++ihit ){

        StrawIndex idx = nearest[ihit];
        //               for  (vector<StepPointMCCollection>::size_type ithhit=0; ithhit < (*_hits)->size(); ++ithhit)
               for  (vector<StepPointMCCollection>::size_type ithhit=0; ithhit < _hits->size(); ++ithhit)
          {
            //            const StepPointMC&  nextHit = (**_hits)[ithhit];//dereference handle to get pointer,then dereference pointer
            const StepPointMC&  nextHit = (*_hits)[ithhit];//dereference pointer
            if ( nextHit.volumeId() == idx.asUint() ){
              Candidate _nextCand;
              _nextCand.id = idx.asInt();
              _nextCand.hitPointer = _hit;
              nearbyStraws.push_back(_nextCand);
            }
          }
      }
    return nearbyStraws;
  }

    void HitCluster::matchAndMerge(bool& match,std::vector<HitCluster>& finalClusters)
    {
      //this is owned by a trialCluster, so I can just get its list of straws and start comparing
      match = false;
      for (std::vector<HitCluster>::size_type ithCluster = 0; (ithCluster <= finalClusters.size()-1) && !match;
           ++ithCluster)
        {
          vector<Candidate> firstSet = finalClusters.at(ithCluster).getStraws();
          for (hitNeighbours::size_type ithStraw = 0; (ithStraw <= firstSet.size()-1) && !match; ++ithStraw)
            {
              int lookingForThisStraw = firstSet.at(ithStraw).id;
              //loop over straws in our test cluster, and note the &&!match -- once there's a match
              //we immediately merge the clusters and there's no reason to continue the loop.
              for (hitNeighbours::size_type jthStraw = 0; (jthStraw <= listOfStraws.size()-1) && !match; ++jthStraw)
                {
                  int trialStraw = listOfStraws.at(jthStraw).id;
                  if (trialStraw == lookingForThisStraw)
                    {
                      //omg, a match; this cluster should get merged
                      match = true;
                      //gobble up all straws in this cluster and stick them in the finalCluster.  It would be
                      //nice to not have to go through the whole list but there's no guarantee the first straw
                      //will be the one that matches.
                      //doing it this way will have duplicates; I can check now, but then I will but repeating
                      //searches every time I look at a new trial cluster.  So I will save this to the end
                      //and crunch down just once.
                      for (hitNeighbours::size_type eatThisStraw = 0; eatThisStraw <= listOfStraws.size()-1;
                           ++eatThisStraw)
                        {
                          finalClusters.at(ithCluster).addStraw(listOfStraws.at(eatThisStraw));
                        }
                    }
                }
            }
        }
      return;
    }


    void HitCluster::cleanUpDuplicates()
    {
      //cleanup the list of straws by removing duplicates; this is easiest for containers
      //since erase screws up iterators

      sort(listOfStraws.begin(),listOfStraws.end());
      hitNeighbours::iterator where = unique(listOfStraws.begin(),listOfStraws.end());
      listOfStraws.erase(where,listOfStraws.end());
      return;
    }






    HitCluster::hitNeighbours HitCluster::getStraws()
    {
      return listOfStraws;
    }

    void HitCluster::addStraw(Candidate addThisStraw)
    {
      listOfStraws.push_back(addThisStraw);
      return;
    }

    int HitCluster::getStraw(int ithStraw)
    {
      return listOfStraws.at(ithStraw).id;
    }



  }     // namespace hitcluster
}       // namespace mu2e
