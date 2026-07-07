#include "Offline/TrkHitReco/inc/Chi2Clusterer.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TMath.h"

#include <algorithm>
#include <vector>
#include <queue>


namespace mu2e
{
  Chi2Clusterer::Chi2Clusterer(const std::optional<Config> config) :
    dhit_          (config.value().hitDistance()),
    dseed_         (config.value().seedDistance()),
    minClusterHits_(config.value().minClusterHits()),
    maxHitdt_      (config.value().maxHitTimeDiff()),
    maxNiter_      (config.value().maxCluIterations()),
    bkgmask_       (config.value().bkgmsk()),
    sigmask_       (config.value().sigmsk()),
    testflag_      (config.value().testflag()),
    minnhits_      (config.value().minActiveHits() ),
    minnp_         (config.value().minNPlanes()),
    kerasW_        (config.value().kerasWeights()),
    useSLine_      (config.value().useSLine()),
    diag_          (config.value().diag())
  {

    tbin_     = 0.;
    tmin_     = 1800.;
    tmax_     = 0.;
    chi2Cut_  = 1.;
    distMethodFlag_ = BkgCluster::chi2;

  }

  //---------------------------------------------------------------------------------------
  void Chi2Clusterer::init()
  {
     ConfigFileLookupPolicy configFile;

     auto kerasWgtsFile = configFile(kerasW_);
     switch ( useSLine_ ){
       case 0 :  sofiePtr1_ = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>(kerasWgtsFile);break;
       case 1 :  sofiePtr2_ = std::make_shared<TMVA_SOFIE_TrainBkgDiagStationChi2SLine::Session>(kerasWgtsFile);break;
     }
  }


  //----------------------------------------------------------------------------------------------------------
  void Chi2Clusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol)
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



  //---------------------------------------------------------------------------------------
  void Chi2Clusterer::classifyCluster(BkgCluster& cluster, const ComboHitCollection& chcol)
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

    std::vector<float> kerasout;
    switch ( useSLine_ ){
     case 0 : kerasout = sofiePtr1_->infer(kerasvars.data());break;
     case 1 : kerasout = sofiePtr2_->infer(kerasvars.data());break;
    }

    cluster.setKerasQ(kerasout[0]);

    if (diag_>0)std::cout << "kerasout = " << kerasout[0] << std::endl;
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
