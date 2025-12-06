#include "Offline/TrkHitReco/inc/DBSClusterer.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include <algorithm>
#include <vector>
#include <queue>


namespace mu2e
{
  DBSClusterer::DBSClusterer(const std::optional<Config> config) :
    DBSminExpand_     (config.value().DBSminN()),
    deltaTime_        (config.value().hitDeltaTime()),
    deltaZ_           (config.value().hitDeltaZ()),
    deltaXY2_         (config.value().hitDeltaXY()*config.value().hitDeltaXY()),
    minClusterHits_   (config.value().minClusterHits()),
    bkgmask_          (config.value().bkgmsk()),
    sigmask_          (config.value().sigmsk()),
    testflag_         (config.value().testflag()),
    kerasW_           (config.value().kerasWeights()),
    diag_             (config.value().diag())
  {}


  //---------------------------------------------------------------------------------------
  void DBSClusterer::init() {
    //Add Classifier init here (see TNTClsuterer for example)
  }


  //----------------------------------------------------------------------------------------------------------
  void DBSClusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol)
  {
     if (chcol.empty()) return;

     std::vector<unsigned> idx;
     idx.reserve(chcol.size());

     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       idx.emplace_back(ich);
     }
     std::sort(idx.begin(),idx.end(),[&chcol](auto i, auto j){return chcol[i].correctedTime() < chcol[j].correctedTime();});


     //--- DBScan algorithm
     const unsigned        noiseID(chcol.size()+1u);
     const unsigned        unprocessedID(chcol.size()+2u);
     unsigned              currentClusterID(0);
     std::queue<unsigned>  inspect;
     std::vector<unsigned> hitToCluster(idx.size(),unprocessedID);
     std::vector<unsigned> neighbors;
     neighbors.reserve(32);

     for (size_t i=0;i<idx.size();++i) {

       // if a point has already been assigned to a cluster, continue
       if ( hitToCluster[i] != unprocessedID) continue;

       // If the neighborhood is too sparse, assign it to noise
       unsigned nNeighbors = findNeighbors(i, idx, chcol, neighbors);
       if (nNeighbors < DBSminExpand_) {
         hitToCluster[i] = noiseID;
         continue;
       }

       hitToCluster[i] = currentClusterID;
       BkgCluster thisCluster;
       thisCluster.addHit(idx[i]);
       // Extend the cluster by adding / expanding around neighbors
       for (const auto& j : neighbors) inspect.push(j);

       while (!inspect.empty()){
         auto j = inspect.front();
         inspect.pop();

         if (hitToCluster[j] == noiseID) {
           hitToCluster[j] = currentClusterID;
           thisCluster.addHit(idx[j]);
         }
         if (hitToCluster[j] != unprocessedID) continue;

         hitToCluster[j] = currentClusterID;
         thisCluster.addHit(idx[j]);

         auto nNeighbors = findNeighbors(j,idx,chcol,neighbors);
         if (nNeighbors < DBSminExpand_) continue;

         for (const auto& k : neighbors) {
           if (hitToCluster[k]==unprocessedID || hitToCluster[k]==noiseID) inspect.push(k);
         }
       }

       if (thisCluster.hits().size() < minClusterHits_) continue;
       clusters.push_back(std::move(thisCluster));
       ++currentClusterID;
     }

     /*
     //Make noise hits into their own cluster if we need to
     for (size_t i=0;i<idx.size();++i) {
       if (hitToCluster[i] != noiseID) continue;
       clusters.emplace_back(BkgCluster());
       clusters.back().addHit(idx[i]);
     }
     */

     //Calculate the cluster properties
     for (auto& cluster : clusters) calculateCluster(cluster, chcol);

     if (diag_>1) dump(clusters);
  }



  //---------------------------------------------------------------------------------------
  // Find the neighbors of given a point - can use any suitable distance function
  unsigned DBSClusterer::findNeighbors(unsigned ihit, const std::vector<unsigned>& idx, const ComboHitCollection& chcol, std::vector<unsigned>& neighbors)
  {
    // Define number of neighbors as number of straw hits in the neighboring point.
    // Commneted version is number of neighbor is 1 regardless of number of straw hits

    //unsigned nNeighbors(0);
    unsigned nNeighbors(chcol[idx[ihit]].nStrawHits()-1);
    neighbors.clear();
    float time0 = chcol[idx[ihit]].correctedTime();
    float x0    = chcol[idx[ihit]].pos().x();
    float y0    = chcol[idx[ihit]].pos().y();
    float z0    = chcol[idx[ihit]].pos().z();

    unsigned istart(ihit);
    while (istart>0 && time0 - chcol[idx[istart]].correctedTime() < deltaTime_) --istart;
    if    (time0 - chcol[idx[istart]].correctedTime()             > deltaTime_) ++istart;

    for (size_t j=istart; j<idx.size(); ++j){
      if (j==ihit) continue;
      if (chcol[idx[j]].correctedTime() - time0 > deltaTime_) break;
      if (std::abs(chcol[idx[j]].pos().z()- z0) > deltaZ_)    continue;

      float dist = (chcol[idx[j]].pos().x()-x0)*(chcol[idx[j]].pos().x()-x0) +
                   (chcol[idx[j]].pos().y()-y0)*(chcol[idx[j]].pos().y()-y0);
      if (dist > deltaXY2_) continue;

      neighbors.emplace_back(j);
      nNeighbors += chcol[idx[j]].nStrawHits();
      //++nNeighbors;
    }

    return nNeighbors;
  }


  //---------------------------------------------------------------------------------------
  // this is only used for diagnosis at this point
  float DBSClusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const
  {
    float psep_x = hit.pos().x()-cluster.pos().x();
    float psep_y = hit.pos().y()-cluster.pos().y();
    return sqrt(psep_x*psep_x+psep_y*psep_y);
    // alterntively, could use the chi2 distance
    //return std::sqrt(cluster.points().dChi2(TwoDPoint(hit.pos(),hit.uDir(),hit.uVar(),hit.vVar())))
    //         - std::sqrt(cluster.points().chisquared());
  }



  //---------------------------------------------------------------------------------------
  void DBSClusterer::calculateCluster(BkgCluster& cluster, const ComboHitCollection& chcol)
  {
    if (cluster.hits().empty()) {cluster.time(0.0f);cluster.pos(XYZVectorF(0.0f,0.0f,0.0f));return;}

    if (cluster.hits().size()==1) {
      int idx = cluster.hits().at(0);
      cluster.time(chcol[idx].correctedTime());
      cluster.pos(XYZVectorF(chcol[idx].pos().x(),chcol[idx].pos().y(),chcol[idx].pos().z()));
      XYZVectorF hitpos(chcol[idx].pos().x(), chcol[idx].pos().y(), chcol[idx].pos().z());
      cluster.addHitPosition(hitpos);
      return;
    }

    float sumWeight(0),crho(0),cphi(0),ctime(0), cz(0);
    for (auto& idx : cluster.hits()) {
      float weight = chcol[idx].nStrawHits();
      float dt     = chcol[idx].correctedTime();
      float dr     = sqrtf(chcol[idx].pos().perp2());
      float dz     = chcol[idx].pos().z();
      XYZVectorF hitpos(chcol[idx].pos().x(), chcol[idx].pos().y(), chcol[idx].pos().z());
      cluster.addHitPosition(hitpos);
      float dp     = chcol[idx].phi();
      if (dp > M_PI)  dp -= 2*M_PI;
      if (dp < -M_PI) dp += 2*M_PI;

      ctime    += dt*weight;
      crho     += dr*weight;
      cphi     += dp*weight;
      cz       += dz*weight;
      sumWeight += weight;
    }

    crho  /= sumWeight;
    cphi  /= sumWeight;
    ctime /= sumWeight;
    cz    /= sumWeight;
    cluster.time(ctime);
    cluster.pos(XYZVectorF(crho*cos(cphi),crho*sin(cphi),cz));
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::classifyCluster(BkgCluster& cluster, const ComboHitCollection& chcol){

     //code logic to classify cluster with MVA

     cluster.setKerasQ(-1.0);
  }





  //-------------------------------------------------------------------------------------------
  void DBSClusterer::dump(const std::vector<BkgCluster>& clusters)
  {
    int iclu(0);
    for (auto& cluster: clusters) {
      std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  "<<cluster.hits().size()<<"  - ";
      for (auto& hit : cluster.hits()) std::cout<<hit<<" ";
      std::cout<<std::endl;
      ++iclu;
    }
  }

}
