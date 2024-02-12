#include "Offline/TrkHitReco/inc/Chi2Clusterer.hh"
#include <algorithm>
#include <vector>
#include <queue>

namespace mu2e
{
  Chi2Clusterer::Chi2Clusterer(const std::optional<Config> config) :
    dhit_             (config.value().hitDistance()),
    dseed_            (config.value().seedDistance()),
    minClusterHits_   (config.value().minClusterHits()),
    maxHitdt_         (config.value().maxHitTimeDiff()),
    maxNiter_         (config.value().maxCluIterations()),
    bkgmask_          (config.value().bkgmsk()),
    sigmask_          (config.value().sigmsk()),
    testflag_         (config.value().testflag()),
    diag_             (config.value().diag())
  {

    tbin_     = 0.;
    tmin_     = 1800.;
    tmax_     = 0.;
    chi2Cut_  = 1.;
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

    auto transPred = [&BkgHits] (const int i) {return BkgHits[i].chidx_;};
    for (auto& cluster: clusters) std::transform(cluster.hits().begin(),cluster.hits().end(),cluster.hits().begin(),transPred);

    if(diag_>1) dump(clusters,BkgHits);

  }


  //----------------------------------------------------------------------------------------------------------------------
  void Chi2Clusterer::initClustering(const ComboHitCollection& chcol, std::vector<Chi2BkgHit>& BkgHits)
  {
     float sumDriftTime(0);
     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       BkgHits.emplace_back(Chi2BkgHit(ich));
       tmin_ = std::min(tmin_,chcol[ich].correctedTime());
       tmax_ = std::max(tmax_,chcol[ich].correctedTime());
       sumDriftTime += chcol[ich].driftTime();
     }
     tbin_ = (sumDriftTime/chcol.size())/2.;
     numBuckets_ = (tmax_ - tmin_ + tbin_)/tbin_;

     auto resPred = [&chcol](const Chi2BkgHit& x, const Chi2BkgHit& y) {return chcol[x.chidx_].wireRes() < chcol[y.chidx_].wireRes();};
     std::sort(BkgHits.begin(),BkgHits.end(),resPred);
  }


  //----------------------------------------------------------------------------------------------------------------------
  void Chi2Clusterer::doClustering(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<Chi2BkgHit>& BkgHits)
  {
    unsigned niter(0);
    float prevChi2(0);
    float deltaChi2(1000);
    //Check average Chi2 change.
    while ( deltaChi2 > chi2Cut_ && niter < maxNiter_ ) {
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
      deltaChi2 = aveChi2 - prevChi2;
      prevChi2 = aveChi2;
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
    std::vector<std::vector<size_t>> clusterIndices;
    clusterIndices.resize(numBuckets_);
    for (size_t ic=0;ic<clusters.size();++ic) {
      size_t itime  = size_t((clusters[ic].time()-tmin_)/tbin_);
      clusterIndices[itime].emplace_back(ic);
      clusterIndices[itime].shrink_to_fit();
    }

    size_t ditime(size_t(maxHitdt_/tbin_));
    for (size_t ihit=0;ihit<BkgHits.size();++ihit) {

      auto& hit  = BkgHits[ihit];
      const auto& chit = chcol[hit.chidx_];

      if(hit.clusterIdx_ == -1){
        // -- find closest cluster. restrict search to clusters close in time using the clusterIndices structure
        int minc(-1);
        float mindist(dseed_ + 1.0f);
        size_t itime = size_t((chit.correctedTime()-tmin_)/tbin_);
        size_t imin = std::max((size_t)0,itime-ditime);
        size_t imax = std::min(numBuckets_,itime+ditime+1);

        for (size_t i=imin;i<imax;++i) {
          for (const auto& ic : clusterIndices[i]) {
            float dist = distance(clusters[ic],chcol[hit.chidx_]);
            if (dist < mindist) {mindist = dist;minc = ic;}
          }
        }

        // -- either add hit to existing cluster, form new cluster, or do nothing if hit is "in between"
        if (mindist < dhit_)
          clusters[minc].addHit(ihit,TwoDPoint(chit.pos(),chit.uDir(),chit.uVar(),chit.vVar()));
        else if (mindist > dseed_) {
          minc = clusters.size();
          clusters.emplace_back(BkgCluster(TwoDPoint(chit.pos(),chit.uDir(),chit.uVar(),chit.vVar()),distMethodFlag_));
          clusters[minc].addHit(ihit);
          size_t itime = size_t((chit.correctedTime()-tmin_)/tbin_);
          clusterIndices[itime].emplace_back(minc);
        }
        else minc = -1;
      }
    }

    for(auto cluster = clusters.begin(); cluster != clusters.end();){
      if((*cluster).hits().size() <= minClusterHits_){
        BkgHits[(*cluster).hits().at(0)].clusterIdx_ = -1;
        cluster = clusters.erase(cluster);
      }
      else{
        //update cluster idx and cluster time
        float weight(0.),wtime(0.);
        for (auto& hit : (*cluster).hits()){
          BkgHits[hit].clusterIdx_ = std::distance(clusters.begin(),cluster);
          weight += chcol[BkgHits[hit].chidx_].nStrawHits();
          wtime += chcol[BkgHits[hit].chidx_].correctedTime()*chcol[BkgHits[hit].chidx_].nStrawHits();
        }
        (*cluster)._time = wtime/weight;
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
    size_t iclu(0);
    for (auto& cluster: clusters) {
      std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  "<<cluster.hits().size()<<"\n"
      <<"Cluster chi2 "<<cluster.points().chisquared()<<" cluster consistency "<<cluster.points().consistency()<<"  - ";
      for (auto& hit : cluster.hits()) std::cout<<BkgHits[hit].chidx_<<" ";
      std::cout<<std::endl;
      ++iclu;
    }
  }

}
