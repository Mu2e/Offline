#include "Offline/TrkHitReco/inc/DBSClusterer.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/TrkHitReco/inc/TrainBkgDiag.hxx"

#include <algorithm>
#include <vector>
#include <limits>


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
  {
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::init() {
    // Add Classifier init here (see TNTClusterer for example)
     ConfigFileLookupPolicy configFile;
     auto kerasWgtsFile = configFile(kerasW_);
     sofiePtr_          = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>(kerasWgtsFile);
  }


  //----------------------------------------------------------------------------------------------------------
  void DBSClusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol)
  {
     if (chcol.empty()) return;

     std::vector<unsigned> idx; // list of combo hit IDs
     idx.reserve(chcol.size());
     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       idx.emplace_back(ich);
     }
     // Sort combo hits which are not flagged as background according to their corrected Time
     std::sort(idx.begin(),idx.end(),[&chcol](auto i, auto j){
       return chcol[i].correctedTime() < chcol[j].correctedTime();
     });

     const unsigned        noiseID(chcol.size()+1u);
     const unsigned        unprocessedID(chcol.size()+2u);
     unsigned              currentClusterID(0);
     // Using DFS(stack) instead of BFS(queue)
     // For DBSCAN, cluster membership of core points is invariant,
     // but traversal order may affect assignment of border points
     // in rare ambiguous cases. This change was validated to not
     // impact physics performance while improving cache locality.
     std::vector<unsigned> inspect;
     inspect.reserve(idx.size());
     std::vector<unsigned> hitToCluster(idx.size(),unprocessedID); // number of combohits used for clustering, cluster ID
     std::vector<unsigned> neighbors;
     neighbors.reserve(256);
     clusters.reserve(std::max(16UL, idx.size()/10));
     unsigned nNeighbors = 0;
     for (size_t i=0;i<idx.size();++i) {
       // If a point has already been assigned to a cluster, continue
       if ( hitToCluster[i] != unprocessedID) continue;

       // If the neighborhood is too sparse, assign it to noise
       nNeighbors = findNeighbors(i, idx, chcol, neighbors);
       if (nNeighbors < DBSminExpand_) {
         hitToCluster[i] = noiseID;
         continue;
       }

       hitToCluster[i] = currentClusterID;
       BkgCluster thisCluster;
       thisCluster.addHit(idx[i]);
       // Extend the cluster by adding/expanding around neighbors
       inspect.clear();
       for (const auto& j : neighbors) inspect.push_back(j);

       while (!inspect.empty()){
         auto j = inspect.back();
         inspect.pop_back();

         if (hitToCluster[j] == noiseID) {
           hitToCluster[j] = currentClusterID;
           thisCluster.addHit(idx[j]);
         }
         if (hitToCluster[j] != unprocessedID) continue;

         hitToCluster[j] = currentClusterID;
         thisCluster.addHit(idx[j]);

         nNeighbors = findNeighbors(j,idx,chcol,neighbors);
         if (nNeighbors >= DBSminExpand_){
           for (const auto& k : neighbors) {
             if (hitToCluster[k] == unprocessedID || hitToCluster[k] == noiseID)
               inspect.push_back(k);
           }
         }
       }
       if (thisCluster.hits().size() >= minClusterHits_){
         clusters.push_back(std::move(thisCluster));
         ++currentClusterID;
       }
     }
     // Calculate the cluster properties
     for (auto& cluster : clusters) calculateCluster(cluster, chcol);
  }

  //---------------------------------------------------------------------------------------
  // Find the neighbors of given a point - can use any suitable distance function
  unsigned DBSClusterer::findNeighbors(unsigned ihit, const std::vector<unsigned>& idx, const ComboHitCollection& chcol, std::vector<unsigned>& neighbors)
  {
    neighbors.clear();
    const auto& hit0 = chcol[idx[ihit]];
    float time0 = hit0.correctedTime();
    float x0    = hit0.pos().x();
    float y0    = hit0.pos().y();
    float z0    = hit0.pos().z();
    unsigned   nNeighbors = 0;
    if (hit0.nStrawHits() > 0) nNeighbors = hit0.nStrawHits() - 1;
    float minTime = time0 - deltaTime_;
    // Use binary search for O(log N) entry into the time-sorted index vector.
    // idx[i] contains the hit index; chcol[idx[i]] is sorted by correctedTime.
    auto it_start = std::partition_point(idx.begin(), idx.end(), [&chcol, minTime](unsigned i){
      return chcol[i].correctedTime() < minTime;
    });
    size_t istart = std::distance(idx.begin(),it_start);
    for (size_t j = istart; j < idx.size(); ++j){
      if (j == ihit) continue;
      const auto& hitj = chcol[idx[j]];
      float dt = hitj.correctedTime() - time0;
      if (dt > deltaTime_) break;
      // Time is constrained by partition_point (backward) and the break (forward)
      // Now check Spatial constraints
      if (std::abs(hitj.pos().z() - z0) > deltaZ_) continue;
      float dx = hitj.pos().x() - x0;
      float dy = hitj.pos().y() - y0;
      float distsq = (dx*dx) + (dy*dy);
      if (distsq <= deltaXY2_){
        neighbors.emplace_back(j);
        nNeighbors += hitj.nStrawHits();
      }
    }
    return nNeighbors;
  }


  //---------------------------------------------------------------------------------------
  // This is only used for diagnosis at this point
  float DBSClusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const
  {
    float psep_x = hit.pos().x()-cluster.pos().x();
    float psep_y = hit.pos().y()-cluster.pos().y();
    return sqrt(psep_x*psep_x+psep_y*psep_y);
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::calculateCluster(BkgCluster& cluster, const ComboHitCollection& chcol)
  {
    if (cluster.hits().empty()) {cluster.time(0.0f);cluster.pos(XYZVectorF(0.0f,0.0f,0.0f));return;}
    if (cluster.hits().size()==1) {
      int idx = cluster.hits().at(0);
      cluster.time(chcol[idx].correctedTime());
      cluster.edep(chcol[idx].energyDep());
      cluster.pos(XYZVectorF(chcol[idx].pos().x(),chcol[idx].pos().y(),chcol[idx].pos().z()));
      XYZVectorF hitpos(chcol[idx].pos().x(), chcol[idx].pos().y(), chcol[idx].pos().z());
      cluster.addHitPosition(hitpos);
      return;
    }
    float sumWeight(0),crho(0),ctime(0), cz(0), cedep(0), cphi(0);
    float phi_ref = chcol[cluster.hits().at(0)].phi();
    for (auto& hitIdx : cluster.hits()) {
      float weight = chcol[hitIdx].nStrawHits();
      float dt     = chcol[hitIdx].correctedTime();
      float dr     = sqrtf(chcol[hitIdx].pos().perp2());
      float dz     = chcol[hitIdx].pos().z();
      float edep   = chcol[hitIdx].energyDep();

      XYZVectorF hitpos = chcol[hitIdx].pos();
      cluster.addHitPosition(hitpos);

      float dp   = hitpos.phi();
      float dphi = dp - phi_ref;
      if (dphi > M_PI)  dphi -= 2*M_PI;
      if (dphi < -M_PI) dphi += 2*M_PI;
      float correctedPhi = phi_ref + dphi;

      ctime    += dt*weight;
      crho     += dr*weight;
      cphi     += correctedPhi*weight;
      cz       += dz*weight;
      cedep    += edep*weight;
      sumWeight += weight;
    }
    cphi  /= sumWeight;
    crho  /= sumWeight;
    ctime /= sumWeight;
    cz    /= sumWeight;
    cedep /= sumWeight; //Weighted average energy deposition of a cluster

    if (cphi > M_PI)  cphi -= 2*M_PI;
    if (cphi < -M_PI) cphi += 2*M_PI;

    cluster.time(ctime);
    cluster.pos(XYZVectorF(crho*cos(cphi),crho*sin(cphi),cz));
    cluster.edep(cedep);
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::classifyCluster(BkgCluster& cluster, const ComboHitCollection& chcol){

    // Code logic to classify cluster with MVA
    // Clusters with less than 3 combo hits have a default keras quality of 0.0
    // and they are not flagged as background clusters
    if(cluster.hits().size() < 3) {
      cluster.setKerasQ(0.0);
      return;
    }
    // find averages
    double sqrSumDeltaTime(0.),sqrSumDeltaX(0.), sqrSumDeltaY(0.), sqrSumDeltaPhi(0.);
    unsigned nhits(0);
    float zmin = std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();
    float phimin = std::numeric_limits<float>::max();
    float phimax = -std::numeric_limits<float>::max();
    float phiclust = cluster.pos().phi();
    if(phiclust > M_PI) phiclust -=2*M_PI;
    if(phiclust < -M_PI) phiclust +=2*M_PI;
    // Safe: Clusters with < 3 hits are returned earlier
    unsigned nchits = cluster.hits().size();
    for (const auto& chit : cluster.hits()) {
      const auto& hit = chcol[chit];
      nhits += hit.nStrawHits();
      float hZ = hit.pos().Z();
      if (hZ < zmin) zmin = hZ;
      if (hZ > zmax) zmax = hZ;
      float dx = hit.pos().x() - cluster.pos().x();
      float dy = hit.pos().y() - cluster.pos().y();
      float dt = hit.correctedTime() - cluster.time();
      sqrSumDeltaX    += dx*dx;
      sqrSumDeltaY    += dy*dy;
      sqrSumDeltaTime += dt*dt;
      float phihit = hit.phi();
      float dphi_rel = phihit- phiclust;
      if(dphi_rel > M_PI)  dphi_rel -= 2*M_PI;
      if(dphi_rel < -M_PI) dphi_rel += 2*M_PI;
      if(dphi_rel < phimin) phimin = dphi_rel;
      if(dphi_rel > phimax) phimax = dphi_rel;
      sqrSumDeltaPhi  += dphi_rel*dphi_rel;
    }
    // Fill mva input variables
    std::array<float,7> kerasvars;
    kerasvars[0] = cluster.pos().Rho(); // cluster rho, cyl coor
    kerasvars[1] = zmax - zmin; // zdiff
    kerasvars[2] = phimax - phimin; // phidiff;
    kerasvars[3] = nhits;
    kerasvars[4] = std::sqrt((sqrSumDeltaX+sqrSumDeltaY)/nchits); // RMS of cluster rho
    kerasvars[5] = std::sqrt(sqrSumDeltaTime/nchits); // RMS of cluster time
    kerasvars[6] = std::sqrt(sqrSumDeltaPhi/nchits); // RMS of cluster phi
    std::vector<float> kerasout = sofiePtr_->infer(kerasvars.data());
    cluster.setKerasQ(kerasout[0]);
  }

}
