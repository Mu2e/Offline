#include "Offline/TrkHitReco/inc/Chi2Clusterer.hh"
#include <algorithm>
#include <vector>
#include <queue>

namespace mu2e
{
  Chi2Clusterer::Chi2Clusterer(const Config& config) :
    dhit_             (config.hitDistance()),
    dseed_            (config.seedDistance()),
    minClusterHits_   (config.minClusterHits()),
    maxHitdt_         (config.maxHitTimeDiff()),
    maxNiter_         (config.maxCluIterations()),
    bkgmask_          (config.bkgmsk()),
    sigmask_          (config.sigmsk()),
    testflag_         (config.testflag()),
    diag_             (config.diag())
  {

    tbin_     = 1.0;
    chi2Cut_  = std::pow(((dseed_+dhit_)/2),2);
    distMethodFlag_ = BkgCluster::chi2;

  }

  //---------------------------------------------------------------------------------------
  void Chi2Clusterer::init() {}


  //----------------------------------------------------------------------------------------------------------
  void Chi2Clusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol, int iev)
  {

    std::vector<Chi2BkgHit> BkgHits;
    BkgHits.reserve(chcol.size());
    if (chcol.empty()) return;

    initClustering(chcol, BkgHits);
    doClustering(chcol, clusters, BkgHits);

    long unsigned int minchits = minClusterHits_;
    auto removePred = [minchits](auto& cluster) {return cluster.hits().size() < minchits;};
    clusters.erase(std::remove_if(clusters.begin(),clusters.end(),removePred),clusters.end());

    auto transPred = [&BkgHits] (const int i) {return BkgHits[i].chidx_;};
    for (auto& cluster: clusters) std::transform(cluster.hits().begin(),cluster.hits().end(),cluster.hits().begin(),transPred);

    if(diag_>1) dump(clusters,BkgHits);

  }


  //----------------------------------------------------------------------------------------------------------------------
  void Chi2Clusterer::initClustering(const ComboHitCollection& chcol, std::vector<Chi2BkgHit>& BkgHits)
  {
     float maxTime(0);
     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       BkgHits.emplace_back(Chi2BkgHit(ich));
       maxTime = std::max(maxTime,chcol[ich].correctedTime());
     }
     tbin_ = (maxTime+1.0)/float(numBuckets_);

     auto resPred = [&chcol](const Chi2BkgHit& x, const Chi2BkgHit& y) {return chcol[x.chidx_].wireRes() < chcol[y.chidx_].wireRes();};
     std::sort(BkgHits.begin(),BkgHits.end(),resPred);
  }


  //----------------------------------------------------------------------------------------------------------------------
  void Chi2Clusterer::doClustering(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<Chi2BkgHit>& BkgHits)
  {
    unsigned niter(0);
    float startChi2(0);
    float deltaChi2(0);
    //Check average Chi2 change. This is set to middle point of hit and seed radii(dhit_+dseed_/2)^2
    while ( deltaChi2 < chi2Cut_ && niter < maxNiter_ ) {
      ++niter;
      formClusters(chcol, clusters, BkgHits);

      float aveChi2(0);
      float sumChi2(0);
      unsigned nofCl(0);
      for (const auto& cluster: clusters) {
        if(cluster.hits().size() > 1) {
          sumChi2 += cluster.points().chisquared();
          ++nofCl;
        }
      }
      aveChi2 = sumChi2/nofCl;
      if(niter == 1) startChi2 = aveChi2;
      deltaChi2 = aveChi2 - startChi2;
    }
    if(diag_>0) std::cout<<"Nof_iter\t"<<niter<<"\tTotal_nof_clusters\t"<<clusters.size()<<"\tAverage_delta_chi2\t"<<deltaChi2<<"\n";

  }

  //-------------------------------------------------------------------------------------------------------------------
  // ChiClusterer strategy;
  // Create a new cluster if there is none. For new hits, check how much cluster chi2 is going to change when a new hit is added for all existing clusters.
  // For the minimum change, if the change is within dhit_ parameter, add it to the cluster. Create a new cluster if change is greater than dseed_.
  // Remove clusters with only single combo hit and try again to place them into existing clusters.
  unsigned Chi2Clusterer::formClusters(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<Chi2BkgHit>& BkgHits)
  {
    std::array<std::vector<int>, numBuckets_> clusterIndices;
    for (size_t ic=0;ic<clusters.size();++ic) {
      int itime  = int(clusters[ic].time()/tbin_);
      clusterIndices[itime].emplace_back(ic);
    }

    int ditime(int(maxHitdt_/tbin_));
    for (size_t ihit=0;ihit<BkgHits.size();++ihit) {

      auto& hit  = BkgHits[ihit];
      const auto& chit = chcol[hit.chidx_];

      if(hit.clusterIdx_ == -1){
        // -- find closest cluster. restrict search to clusters close in time using the clusterIndices structure
        int minc(-1);
        float mindist(dseed_ + 1.0f);
        int itime = int(chit.correctedTime()/tbin_);
        int imin = std::max(0,itime-ditime);
        int imax = std::min(numBuckets_,itime+ditime+1);

        for (int i=imin;i<imax;++i) {
          for (const auto& ic : clusterIndices[i]) {
            float dist = distance(clusters[ic],chcol[hit.chidx_]);
            if (dist < mindist) {mindist = dist;minc = ic;}
          }
        }

        // -- either add hit to existing cluster, form new cluster, or do nothing if hit is "in between"
        if (mindist < dhit_)
          clusters[minc].addHit(ihit,TwoDPoint(chit.pos(),chit.uDir(),chit.uVar(),chit.vVar()),chit.correctedTime(),chit.nStrawHits());
        else if (mindist > dseed_) {
          minc = clusters.size();
          clusters.emplace_back(BkgCluster(TwoDPoint(chit.pos(),chit.uDir(),chit.uVar(),chit.vVar()),chit.correctedTime(),chit.nStrawHits()));
          clusters[minc].addHit(ihit);
          int itime = int(chit.correctedTime()/tbin_);
          clusterIndices[itime].emplace_back(minc);
          clusters[minc].setDistanceMethod(distMethodFlag_);
        }
        else minc = -1;
      }
    }

    for(auto cluster = clusters.begin(); cluster != clusters.end();){
      if((*cluster).hits().size() <= 1){
        BkgHits[(*cluster).hits().at(0)].clusterIdx_ = -1;
        cluster = clusters.erase(cluster);
      }
      else{
        for (auto& hit : (*cluster).hits()){
          BkgHits[hit].clusterIdx_ = std::distance(clusters.begin(),cluster);
        }
        ++cluster;
      }
    }
    return clusters.size();
  }

  //---------------------------------------------------------------------------------------
  float Chi2Clusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const
  {
    return std::sqrt(cluster.points().dChi2(TwoDPoint(hit.pos(),hit.uDir(),hit.uVar(),hit.vVar()))) - std::sqrt(cluster.points().chisquared());
  }

  //-------------------------------------------------------------------------------------------
  void Chi2Clusterer::dump(const std::vector<BkgCluster>& clusters, const std::vector<Chi2BkgHit>& BkgHits)
  {
    int iclu(0);
    for (auto& cluster: clusters) {
      std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  "<<cluster.hits().size()<<"\n"
      <<"Cluster chi2 "<<cluster.points().chisquared()<<" cluster consistency "<<cluster.points().consistency()<<"  - ";
      for (auto& hit : cluster.hits()) std::cout<<BkgHits[hit].chidx_<<" ";
      std::cout<<std::endl;
      ++iclu;
    }
  }

}
