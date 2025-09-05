#include "Offline/TrkHitReco/inc/TNTClusterer.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TMath.h"

#include <algorithm>
#include <vector>
#include <queue>


namespace mu2e
{
  TNTClusterer::TNTClusterer(const std::optional<Config> config) :
    dhit_          (config.value().hitDistance()),
    dseed_         (config.value().seedDistance()),
    dd_            (config.value().clusterDiameter()),
    dt_            (config.value().clusterTime()),
    minClusterHits_(config.value().minClusterHits()),
    maxHitdt_      (config.value().maxHitTimeDiff()),
    maxDistSum_    (config.value().maxSumDistance()),
    maxNiter_      (config.value().maxCluIterations()),
    useMedian_     (config.value().medianCentroid()),
    comboInit_     (config.value().comboInit()),
    bkgmask_       (config.value().bkgmsk()),
    sigmask_       (config.value().sigmsk()),
    testflag_      (config.value().testflag()),
    minnhits_      (config.value().minActiveHits() ),
    minnp_         (config.value().minNPlanes()),
    kerasW_        (config.value().kerasWeights()),
    diag_          (config.value().diag())
  {
    float minerr (config.value().minHitError());
    float maxdist(config.value().maxDistance());
    float trms   (config.value().timeRMS());

    tbin_     =1.0;
    dd2_      = dd_*dd_;
    maxwt_    = 1.0f/minerr;
    md2_      = maxdist*maxdist;
    trms2inv_ = 1.0f/trms/trms;
    distMethodFlag_ = BkgCluster::spatial;

  }


  //---------------------------------------------------------------------------------------
  void TNTClusterer::init()
  {
     ConfigFileLookupPolicy configFile;
     auto kerasWgtsFile = configFile(kerasW_);
     sofiePtr_          = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>(kerasWgtsFile);
  }


  //----------------------------------------------------------------------------------------------------------
  void TNTClusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol)
  {
    std::vector<BkgHit> BkgHits;
    BkgHits.reserve(chcol.size());
    if (chcol.empty()) return;

    initClustering(chcol, BkgHits);
    doClustering(chcol, clusters, BkgHits);

    long unsigned int minchits = minClusterHits_;
    auto removePred = [minchits](auto& cluster) {return cluster.hits().size() < minchits;};
    clusters.erase(std::remove_if(clusters.begin(),clusters.end(),removePred),clusters.end());

    auto transPred = [&BkgHits] (const int i) {return BkgHits[i].chidx_;};
    for (auto& cluster: clusters) std::transform(cluster.hits().begin(),cluster.hits().end(),cluster.hits().begin(),transPred);
  }


  //----------------------------------------------------------------------------------------------------------------------
  void TNTClusterer::initClustering(const ComboHitCollection& chcol, std::vector<BkgHit>& BkgHits)
  {
     float maxTime(0);
     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       BkgHits.emplace_back(BkgHit(ich));
       maxTime = std::max(maxTime,chcol[ich].correctedTime());
     }
     tbin_ = (maxTime+1.0)/float(numBuckets_);

     auto resPred = [&chcol](const BkgHit& x, const BkgHit& y) {return chcol[x.chidx_].wireRes() < chcol[y.chidx_].wireRes();};
     if (comboInit_) std::sort(BkgHits.begin(),BkgHits.end(),resPred);
  }


  //----------------------------------------------------------------------------------------------------------------------
  void TNTClusterer::doClustering(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& BkgHits)
  {
    unsigned niter(0);
    float odist(2.0f*maxDistSum_);
    float tdist(0.0f);
    while ( std::abs(odist - tdist) > maxDistSum_ && niter < maxNiter_ ) {
      ++niter;
      formClusters(chcol, clusters, BkgHits);

      odist = tdist;
      tdist = 0.0f;
      for (const auto& cluster: clusters) {
        for (const auto& cidx : cluster.hits()) tdist += BkgHits[cidx].distance_;
      }
    }
  }


  //-------------------------------------------------------------------------------------------------------------------
  // loop over hits, re-affect them to their original cluster if they are still within the radius, otherwise look at
  // candidate clusters to check if they could be added. If not, make a new cluster.
  // speed up: don't update clusters who haven't changed + cache cluster id in a given time window in clusterIndex (array if vectors)
  //
  unsigned TNTClusterer::formClusters(const ComboHitCollection& chcol, std::vector<BkgCluster>& clusters, std::vector<BkgHit>& BkgHits)
  {
    std::array<std::vector<int>, numBuckets_> clusterIndices;
    for (size_t ic=0;ic<clusters.size();++ic) {
      int itime  = int(clusters[ic].time()/tbin_);
      clusterIndices[itime].emplace_back(ic);
      clusters[ic].clearHits();
    }

    int ditime(int(maxHitdt_/tbin_));
    unsigned nchanged(0);
    for (size_t ihit=0;ihit<BkgHits.size();++ihit) {

      auto& hit  = BkgHits[ihit];
      const auto& chit = chcol[hit.chidx_];

      // -- if hit is ok, reassign it right away
      if (hit.distance_ < dhit_) {
        clusters[hit.clusterIdx_].addHit(ihit);
        continue;
      }

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
      if (mindist < dhit_) {
        clusters[minc].addHit(ihit);
      }
      else if (mindist > dseed_) {
        minc = clusters.size();
        clusters.emplace_back(chit.pos(),chit.correctedTime(),distMethodFlag_);
        clusters[minc].addHit(ihit);
        int itime = int(chit.correctedTime()/tbin_);
        clusterIndices[itime].emplace_back(minc);
      }
      else{
        BkgHits[ihit].distance_ = 10000.0f;
        minc = -1;
      }

      // -- update cluster flag and hit->cluster pointer if association has changed
      if (minc != -1) {
        if (hit.clusterIdx_ != minc) {
          ++nchanged;
          if (hit.clusterIdx_ != -1) clusters[hit.clusterIdx_]._flag = BkgClusterFlag::update;
          clusters[minc]._flag = BkgClusterFlag::update;
        }
        hit.clusterIdx_ = minc;
      }
      else {
        if (hit.clusterIdx_ != -1) {
          ++nchanged;
          clusters[hit.clusterIdx_]._flag = BkgClusterFlag::update;
        }
        hit.clusterIdx_ = -1;
      }
    }

    // -- update cluster and hit distance if needed
    for (auto& cluster : clusters) {
      if (cluster._flag == BkgClusterFlag::update) {
        cluster._flag = BkgClusterFlag::unchanged;
        updateCluster(cluster, chcol, BkgHits);
        if (cluster.hits().size()==1) {BkgHits[cluster.hits().at(0)].distance_ = 0.0f;}
        else {
          for (auto& hit : cluster.hits()) BkgHits[hit].distance_ = distance(cluster,chcol[BkgHits[hit].chidx_]);
        }
      }
    }

    return nchanged;
  }


  //-----------------------------------------------------------------------------------------------
  void TNTClusterer::mergeClusters(std::vector<BkgCluster>& clusters, const ComboHitCollection& chcol,
                                   std::vector<BkgHit>& BkgHits, float dt, float dd2)
  {
    unsigned niter(0);
    while (niter < maxNiter_) {
      int nchanged(0);
      for (auto it1 = clusters.begin();it1 != std::prev(clusters.end());++it1) {
        if (it1->hits().empty()) continue;
        for (auto it2 = std::next(it1);it2!=clusters.end();++it2) {
          if (it2->hits().empty()) continue;
          if (std::abs(it1->time() - it2->time()) > dt) continue;
          if ((it1->pos() - it2->pos()).perp2() > dd2)  continue;
          ++nchanged;
          mergeTwoClusters(*it1,*it2);
        }
      }

      ++niter;
      if (diag_>0) std::cout<<"Merge "<<niter<<" "<<nchanged<<"  "<<clusters.size()<<std::endl;
      if (nchanged==0) break;

      std::remove_if(clusters.begin(),clusters.end(),[](auto& cluster) {return cluster.hits().empty();});
      for (auto& cluster : clusters ) updateCluster(cluster, chcol, BkgHits);
    }
    return;
  }

  void TNTClusterer::mergeTwoClusters(BkgCluster& clu1, BkgCluster& clu2 )
  {
    if (clu1.hits().empty() || clu2.hits().empty()) return;
    clu1.hits().insert(clu1.hits().end(),clu2.hits().begin(), clu2.hits().end());
    clu2.hits().clear();
  }


  //---------------------------------------------------------------------------------------
  // only count differences if they are above the natural hit size (drift time, straw size)
  float TNTClusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const
  {
    float retval(0.0f);
    float psep_x = hit.pos().x()-cluster.pos().x();
    float psep_y = hit.pos().y()-cluster.pos().y();
    float d2     = psep_x*psep_x+psep_y*psep_y;

    if (d2 > md2_) {return dseed_+1.0f;}

    float dt = std::abs(hit.correctedTime()-cluster.time());
    if (dt > maxHitdt_) {return dseed_+1.0f;}


    if (dt > dt_) {float tdist = dt -dt_;retval = tdist*tdist*trms2inv_;}
    if (d2 > dd2_) {
      //This is equivalent to but faster than the commented lines
      //XYZVectorF that(-hit.uDir2D().y(),hit.uDir2D().x(),0.0);
      //float dw = std::max(0.0f,hit.uDir2D().Dot(psep)-dd_)/hit.posRes(ComboHit::wire);
      //float dp = std::max(0.0f,that.Dot(psep)-dd_)*maxwt_;//maxwt = 1/minerr
      float hwx = hit.uDir2D().x();
      float hwy = hit.uDir2D().y();
      float dw  = std::max(0.0f,(psep_x*hwx+psep_y*hwy-dd_)/hit.posRes(ComboHit::wire));
      float dp  = std::max(0.0f,(hwx*psep_y-hwy*psep_x-dd_)*maxwt_);
      retval += dw*dw + dp*dp;
    }
    return retval;
  }

  void TNTClusterer::updateCluster(BkgCluster& cluster, const ComboHitCollection& chcol, std::vector<BkgHit>& BkgHits)
  {

    if (cluster.hits().empty()) {cluster.time(0.0f);cluster.pos(XYZVectorF(0.0f,0.0f,0.0f));return;}

    if (cluster.hits().size()==1) {
      int idx = BkgHits[cluster.hits().at(0)].chidx_;
      cluster.time(chcol[idx].correctedTime());
      cluster.pos(XYZVectorF(chcol[idx].pos().x(),chcol[idx].pos().y(),0.0f));
      return;
    }

    float crho  = sqrtf(cluster.pos().perp2());
    float cphi  = cluster.pos().phi();
    float ctime = cluster.time();

    if (useMedian_) {
      std::vector<float> racc,pacc,tacc;
      for (auto& hit : cluster.hits()) {
        int idx  = BkgHits[hit].chidx_;
        float dt = chcol[idx].correctedTime() - ctime;
        float dr = sqrtf(chcol[idx].pos().perp2()) - crho;
        float dp = chcol[idx].phi() - cphi;
        if (dp > M_PI)  dp -= 2*M_PI;
        if (dp < -M_PI) dp += 2*M_PI;

        // weight according to the # of hits
        for (int i=0;i<chcol[idx].nStrawHits();++i) {
          racc.emplace_back(dr);
          pacc.emplace_back(dp);
          tacc.emplace_back(dt);
        }
      }

      size_t vecSize = racc.size()/2;
      std::nth_element(racc.begin(),racc.begin()+vecSize,racc.end());
      std::nth_element(pacc.begin(),pacc.begin()+vecSize,pacc.end());
      std::nth_element(tacc.begin(),tacc.begin()+vecSize,tacc.end());

      crho  += racc[vecSize];
      cphi  += pacc[vecSize];
      ctime += tacc[vecSize];
    }
    else
    {
      float sumWeight(0), deltaT(0), deltaP(0), deltaR(0);
      for (auto& hit : cluster.hits()) {
        int   idx    = BkgHits[hit].chidx_;
        float weight = chcol[idx].nStrawHits();
        float dt     = chcol[idx].correctedTime()-ctime;
        float dr     = sqrtf(chcol[idx].pos().perp2()) - crho;

        float dp     = chcol[idx].phi()-cphi;
        if (dp > M_PI)  dp -= 2*M_PI;
        if (dp < -M_PI) dp += 2*M_PI;

        deltaT    += dt*weight;
        deltaR    += dr*weight;
        deltaP    += dp*weight;
        sumWeight += weight;
      }
      crho  += deltaR/sumWeight;
      cphi  += deltaP/sumWeight;
      ctime += deltaT/sumWeight;
    }

    cluster.time(ctime);
    cluster.pos(XYZVectorF(crho*cos(cphi),crho*sin(cphi),0.0f));

  }


  //---------------------------------------------------------------------------------------
  void TNTClusterer::classifyCluster(BkgCluster& cluster, const ComboHitCollection& chcol)
  {
    // count hits and planes
    std::array<int,StrawId::_nplanes> hitplanes{0};
    for (const auto& chit : cluster.hits()) {
    const ComboHit& ch = chcol[chit];
    hitplanes[ch.strawId().plane()] += ch.nStrawHits();
    }

    unsigned npexp(0),np(0),nhits(0);
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


    if(nhits < minnhits_ || np < minnp_) return;

    // find averages
    double sumEdep(0.);
    double sqrSumDeltaTime(0.);
    double sqrSumDeltaX(0.);
    double sqrSumDeltaY(0.);
    double sqrSumQual(0.);
    double sumPitch(0.);
    double sumYaw(0.);
    double sumwPitch(0.);
    double sumwYaw(0.);
    double sumEcc(0.);
    double sumwEcc(0.);
    unsigned nsthits(0.);
    unsigned nchits = cluster.hits().size();
    for (const auto& chit : cluster.hits()) {
      sumEdep         +=  chcol[chit].energyDep()/chcol[chit].nStrawHits();
      sqrSumDeltaX    += std::pow(chcol[chit].pos().x() - cluster.pos().x(),2);
      sqrSumDeltaY    += std::pow(chcol[chit].pos().y() - cluster.pos().y(),2);
      sqrSumDeltaTime += std::pow(chcol[chit].time() - cluster.time(),2);
      auto hdir        = chcol[chit].hDir();
      auto wecc        = chcol[chit].nStrawHits();
      sumEcc          += std::sqrt(1-(chcol[chit].vVar()/chcol[chit].uVar()))*wecc;
      sumwEcc         += wecc;
      if (chcol[chit].flag().hasAllProperties(StrawHitFlag::sline)){

        //quality of SLine fit
        sqrSumQual += std::pow(chcol[chit].qual(),2);

        //angle with Mu2e-Y
        double varPitch = std::pow(TMath::ACos(std::sqrt(chcol[chit].hcostVar())),2);
        double wPitch = 1/varPitch;
        double signPitch = hdir.Y()/std::abs(hdir.Y());
        sumPitch += signPitch*wPitch*hdir.theta();
        sumwPitch += wPitch;

        ROOT::Math::XYZVectorF z = {0,0,1};
        ROOT::Math::XYZVectorF dxdz = {hdir.X(),0,hdir.Z()};
        float magdxdz = std::sqrt(dxdz.Mag2());

        //angle with Mu2e-Z
        double varYaw = std::sqrt(chcol[chit].hphiVar() + varPitch);
        double wYaw = 1/varYaw;
        double signYaw = hdir.X()/std::abs(hdir.X());
        sumYaw += signYaw*wYaw*TMath::ACos(dxdz.Dot(z)/magdxdz);
        sumwYaw += wYaw;

        // # of stereo hits with SLine
        nsthits++;
      }
    }

    // fill mva input variables
    std::array<float,12> kerasvars;
    kerasvars[0]  = cluster.pos().Rho(); // cluster rho, cyl coor
    kerasvars[1]  = fp;// first plane hit
    kerasvars[2]  = lp;// last plane hit
    kerasvars[3]  = pgap;// largest plane gap without hits between planes with hits
    kerasvars[4]  = np;// # of planes hit
    kerasvars[5]  =  static_cast<float>(np)/static_cast<float>(lp - fp);// fraction of planes hit between first and last plane
    kerasvars[6]  = nhits;// sum of straw hits
    kerasvars[7]  = std::sqrt((sqrSumDeltaX+sqrSumDeltaY)/nchits);  // RMS of cluster rho
    kerasvars[8]  = std::sqrt(sqrSumDeltaTime/nchits);// RMS of cluster time
    kerasvars[9]  = nsthits > 0 ? sumPitch/sumwPitch : 0.;
    kerasvars[10] = nsthits > 0 ? sumYaw/sumwYaw : 0.;
    kerasvars[11] = sumEcc/sumwEcc;

    std::vector<float> kerasout = sofiePtr_->infer(kerasvars.data());
    cluster.setKerasQ(kerasout[0]);

    if (diag_>0)std::cout << "kerasout = " << kerasout[0] << std::endl;
  }



  //-------------------------------------------------------------------------------------------
  void TNTClusterer::dump(const std::vector<BkgCluster>& clusters, const std::vector<BkgHit>& BkgHits)
  {
    int iclu(0);
    for (auto& cluster: clusters) {
      std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  "<<cluster.hits().size()<<"  - ";
      for (auto& hit : cluster.hits()) std::cout<<BkgHits[hit].chidx_<<" ";
      std::cout<<std::endl;
      ++iclu;
    }
  }

}
