#include "Offline/TrkHitReco/inc/DBSClusterer.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TMath.h"
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
  {
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::init() {
    //Add Classifier init here (see TNTClsuterer for example)
     ConfigFileLookupPolicy configFile;
     auto kerasWgtsFile = configFile(kerasW_);
     sofiePtr_          = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>(kerasWgtsFile);
  }


  //----------------------------------------------------------------------------------------------------------
  void DBSClusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol)
  {
     if (chcol.empty()) return;

     std::vector<unsigned> idx; //list of combo hit IDs
     idx.reserve(chcol.size());
     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       idx.emplace_back(ich);
     }
     //Sort combo hits which are not flagged as background according to their corrected Time
     std::sort(idx.begin(),idx.end(),[&chcol](auto i, auto j){
       return chcol[i].correctedTime() < chcol[j].correctedTime();
     });

     const unsigned        noiseID(chcol.size()+1u);
     const unsigned        unprocessedID(chcol.size()+2u);
     unsigned              currentClusterID(0);
     //std::queue<unsigned>  inspect;
     std::vector<unsigned> inspect;
     inspect.reserve(idx.size());
     std::vector<unsigned> hitToCluster(idx.size(),unprocessedID); //number of combohits used for clustering, cluster ID
     std::vector<unsigned> neighbors;
     neighbors.reserve(256);
     clusters.reserve(std::max(16UL, idx.size()/10));
     for (size_t i=0;i<idx.size();++i) {
       // if a point has already been assigned to a cluster, continue
       //int nNeighbors = 0;
       if ( hitToCluster[i] != unprocessedID) continue;

       // If the neighborhood is too sparse, assign it to noise
       int nNeighbors = findNeighbors(i, idx, chcol, neighbors);
       if (nNeighbors < DBSminExpand_) {
         hitToCluster[i] = noiseID;
         continue;
       }

       hitToCluster[i] = currentClusterID;
       BkgCluster thisCluster;
       thisCluster.addHit(idx[i]);
       // Extend the cluster by adding / expanding around neighbors
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
     //Calculate the cluster properties
     for (auto& cluster : clusters) calculateCluster(cluster, chcol);
     if (diag_>1) dump(clusters);
  }



  //---------------------------------------------------------------------------------------
  // Find the neighbors of given a point - can use any suitable distance function
  int DBSClusterer::findNeighbors(unsigned ihit, const std::vector<unsigned>& idx, const ComboHitCollection& chcol, std::vector<unsigned>& neighbors)
  {
    neighbors.clear();
    const auto& hit0 = chcol[idx[ihit]];
    float time0 = hit0.correctedTime();
    float x0    = hit0.pos().x();
    float y0    =  hit0.pos().y();
    float z0    = hit0.pos().z();
    unsigned nNeighbors = hit0.nStrawHits()-1;
    float minTime = time0 - deltaTime_;
    auto it_start = std::lower_bound(idx.begin(), idx.end(), minTime, [&chcol](unsigned i, float val){
      return chcol[i].correctedTime() < val;
    });
    size_t istart = std::distance(idx.begin(),it_start);
    for (size_t j = istart; j < idx.size(); ++j){
      if (j == ihit) continue;
      const auto& hitj = chcol[idx[j]];
      float dt = hitj.correctedTime() - time0;
      if (dt > deltaTime_) break;
      if (std::abs(hitj.pos().z() - z0) > deltaZ_) continue;
      float dx = hitj.pos().x() - x0;
      float dy = hitj.pos().y() - y0;
      float distsq = (dx*dx) + (dy*dy);
      if (distsq <= deltaXY2_){
        neighbors.emplace_back(j);
        nNeighbors += hitj.nStrawHits();
      }
    }
    /*unsigned istart(ihit);
    while (istart>0 && time0 - chcol[idx[istart]].correctedTime() < deltaTime_) --istart;
    if    (time0 - chcol[idx[istart]].correctedTime()             > deltaTime_) ++istart;
    for (size_t j=istart; j<idx.size(); ++j){
      if (j==ihit) continue;
      if (chcol[idx[j]].correctedTime() - time0 > deltaTime_) break;     // deltaTime = 15
      if (std::abs(chcol[idx[j]].pos().z()- z0) > deltaZ_)    continue; //deltaZ = 800
      float dist = (chcol[idx[j]].pos().x()-x0)*(chcol[idx[j]].pos().x()-x0) +
                   (chcol[idx[j]].pos().y()-y0)*(chcol[idx[j]].pos().y()-y0);
      if (dist > deltaXY2_) continue;
      neighbors.emplace_back(j);
      nNeighbors += chcol[idx[j]].nStrawHits();
      }*/
    return nNeighbors;
  }


  //---------------------------------------------------------------------------------------
  // this is only used for diagnosis at this point
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
    for (auto& idx : cluster.hits()) {
      float weight = chcol[idx].nStrawHits();
      float dt     = chcol[idx].correctedTime();
      float dr     = sqrtf(chcol[idx].pos().perp2());
      float dz     = chcol[idx].pos().z();
      float edep   = chcol[idx].energyDep();
      XYZVectorF hitpos(chcol[idx].pos().x(), chcol[idx].pos().y(), chcol[idx].pos().z());
      cluster.addHitPosition(hitpos);
      float dp     = chcol[idx].phi();
      ctime    += dt*weight;
      crho     += dr*weight;
      cphi     += dp*weight;
      cz       += dz*weight;
      cedep    += edep*weight;
      sumWeight += weight;
    }
    cphi  /= sumWeight;
    crho  /= sumWeight;
    ctime /= sumWeight;
    cz    /= sumWeight;
    cedep /= sumWeight;
    cluster.time(ctime);
    cluster.pos(XYZVectorF(crho*cos(cphi),crho*sin(cphi),cz));
    cluster.edep(cedep);
  }

  //---------------------------------------------------------------------------------------
  void DBSClusterer::classifyCluster(BkgCluster& cluster, const ComboHitCollection& chcol){

    //code logic to classify cluster with MVA
    // count hits and planes
    std::array<int,StrawId::_nplanes> hitplanes{0};
    for (const auto& chit : cluster.hits()) {
    const ComboHit& ch = chcol[chit];
    hitplanes[ch.strawId().plane()] += ch.nStrawHits();
    }

    unsigned npexp(0),np(0),nhits(0);
    //std::vector<float> hz;
    //std::vector<float> hphi;
    //float phidiff(0.0);
    int ipmin(0),ipmax(StrawId::_nplanes-1);
    while (hitplanes[ipmin]==0 && ipmin<StrawId::_nplanes) ++ipmin;
    while (hitplanes[ipmax]==0 and ipmax>0)                --ipmax;
    int fp(ipmin),lp(ipmin-1),pgap(0);
    for (int ip = ipmin; ip <= ipmax; ++ip) {
    npexp++; // should use TTracker to see if plane is physically present FIXME!
    if (hitplanes[ip]> 0){
      ++np;
      if(lp > 0 && ip - lp -1 > pgap)pgap = ip - lp -1;
      if(ip > lp)lp = ip;
      if(ip < fp)fp = ip;
      lp = ip;
    }
    nhits += hitplanes[ip];
    }
    if(nhits < 1 || np < 2) return;
    // find averages
    double sqrSumDeltaTime(0.),sqrSumDeltaX(0.), sqrSumDeltaY(0.), sqrSumDeltaPhi(0.);
    float zmin = std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();
    float phimin = std::numeric_limits<float>::max();
    float phimax = -std::numeric_limits<float>::max();
    float phiclust = cluster.pos().phi();
    if(phiclust > M_PI) phiclust -=2*M_PI;
    if(phiclust < -M_PI) phiclust +=2*M_PI;
    unsigned nchits = cluster.hits().size();
    for (const auto& chit : cluster.hits()) {
      const auto& hit = chcol[chit];
      float hZ = hit.pos().Z();
      if (hZ < zmin) zmin = hZ;
      if (hZ > zmax) zmax = hZ;
      float dx = hit.pos().x() - cluster.pos().x();
      float dy = hit.pos().y() - cluster.pos().y();
      float dt = hit.time() - cluster.time();
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
    // fill mva input variables
    std::array<float,7> kerasvars;
    kerasvars[0]  = cluster.pos().Rho(); // cluster rho, cyl coor
    kerasvars[1]  = zmax - zmin; //zdiff; //fp;// first plane hit
    kerasvars[2]  = phimax - phimin; //phidiff;
    kerasvars[3]  = nhits;
    kerasvars[4]  = std::sqrt((sqrSumDeltaX+sqrSumDeltaY)/nchits); // RMS of cluster rho
    kerasvars[5] = std::sqrt(sqrSumDeltaTime/nchits);// RMS of cluster time
    kerasvars[6] = std::sqrt(sqrSumDeltaPhi/nchits);// RMS of cluster phi
    std::vector<float> kerasout = sofiePtr_->infer(kerasvars.data());
    cluster.setKerasQ(kerasout[0]);
    if (diag_>0)std::cout << "kerasout = " << kerasout[0] << std::endl;
    //cluster.setKerasQ(-1.0);
  }





  //-------------------------------------------------------------------------------------------
  void DBSClusterer::dump(const std::vector<BkgCluster>& clusters)
  {
    int iclu(0);
    for (auto& cluster: clusters) {
      for (auto& hit : cluster.hits()) std::cout<<hit<<" ";
      std::cout<<std::endl;
      ++iclu;
    }
  }

}
