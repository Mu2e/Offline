//
// Tracker time / phi cluster finder
//
// Original author B. Echenard
//

#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrkPatRec/inc/TimeAndPhiClusterFinder_types.hh"

namespace
{
   struct TimePhiCandidate
   {
       TimePhiCandidate() {strawIdx_.reserve(32);}

       TimePhiCandidate(unsigned nsh, const std::vector<StrawHitIndex>& strawIdx, const art::Ptr<mu2e::CaloCluster>& caloCluster) :
         nsh_(nsh), strawIdx_(strawIdx), caloCluster_(caloCluster)
       {}

       unsigned                    nsh_      = 0;
       std::vector<StrawHitIndex>  strawIdx_ = {};
       art::Ptr<mu2e::CaloCluster> caloCluster_;
   };

   typedef std::vector<TimePhiCandidate> TimePhiCandidateCollection;
}






namespace mu2e {

  class TimeAndPhiClusterFinder : public art::EDProducer
  {
    public:

      using Config_types = TimeAndPhiClusterFinderTypes::Config ;
      using Data_types   = TimeAndPhiClusterFinderTypes::Data_t ;
      struct Config
      {
          using Name    = fhicl::Name;
          using Comment = fhicl::Comment;
          fhicl::Atom<art::InputTag>              comboHitCollection     {Name("ComboHitCollection"),     Comment("ComboHit collection {Name") };
          fhicl::Atom<art::InputTag>              caloClusterCollection  {Name("CaloClusterCollection"),  Comment("Calo cluster collection {Name") };
          fhicl::Table<MVATools::Config>          MVATime                {Name("MVATime"),                Comment("MVA for time cluster cleaning") };
          fhicl::Sequence<std::string>            hsel                   {Name("HitSelectionBits"),       Comment("HitSelectionBits") };
          fhicl::Sequence<std::string>            hbkg                   {Name("HitBackgroundBits"),      Comment("HitBackgroundBits") };
          fhicl::Atom<bool>                       testflag               {Name("TestFlag"),               Comment("Test hit flags") };
          fhicl::Atom<bool>                       usecc                  {Name("UseCaloCluster"),         Comment("Use calorimeter cluster") };
          fhicl::Atom<float>                      ccmine                 {Name("CaloClusterMinE"),        Comment("Minimum energy for calorimeter cluster") };
          fhicl::Atom<float>                      ccweight               {Name("CaloClusterWeight"),      Comment("Weight of cluster in tracker hits") };
          fhicl::Atom<unsigned>                   algoAssignHits         {Name("AlgoAssignHits"),         Comment("Hit assignment algorithm, must be 1 or 2") };
          fhicl::Atom<unsigned>                   minNSHits              {Name("MinNSHits"),              Comment("Minimum number of hits for cluster") };
          fhicl::Atom<float>                      tbin                   {Name("Tbin"),                   Comment("Time histogram bin width") };
          fhicl::Atom<unsigned>                   minTimeYbin            {Name("MinTimeYbin"),            Comment("Minimum number of bins to start recording max for scanning algo") };
          fhicl::Atom<float>                      maxTimeDT              {Name("MaxTimeDT"),              Comment("Max time difference for hits in cluster") };
          fhicl::Atom<bool>                       filterMVA              {Name("FilterMVA"),              Comment("Refine time cluster with NN info") };
          fhicl::Atom<unsigned>                   minHitSelect           {Name("MinHitSelect"),           Comment("Minimun hits to copy full time cluster without any filtering") };
          fhicl::Atom<float>                      minCutMVA              {Name("MinCutMVA"),              Comment("Minimun value of NN output to keep hit") };
          fhicl::Atom<bool>                       splitPhi               {Name("SplitPhi"),               Comment("Split time cluster with phi info") };
          fhicl::Atom<bool>                       reqphi                 {Name("RequirePhi"),             Comment("Don't include time clusters before split into phi clusters") };
          fhicl::Atom<float>                      maxDeltaPhi            {Name("MaxDeltaPhi"),            Comment("Max delta phi between consecutive hits in cluster") };
          fhicl::Atom<int>                        maxNdiff               {Name("MaxNdiff"),               Comment("Max difference of number of hits between duplicated clusters") };
          fhicl::Atom<int>                        diag                   {Name("Diag"),                   Comment("Diag level"), 0 };
          fhicl::Table<Config_types>              diagPlugin             {Name("DiagPlugin"),             Comment("Diag Plugin config")};
      };


      explicit TimeAndPhiClusterFinder(const art::EDProducer::Table<Config>& config);
      ~TimeAndPhiClusterFinder() = default;

      void beginJob() override;
      void produce(art::Event& e) override;


    private:
      int                                             iev_;
      const art::ProductToken<ComboHitCollection>     chToken_;
      const art::ProductToken<CaloClusterCollection>  ccToken_;
      MVATools                                        MVATime_;
      StrawHitFlag                                    hsel_, hbkg_;
      bool                                            testflag_;
      bool                                            usecc_;
      float                                           ccmine_,ccweight_;
      unsigned                                        algoAssignHits_;
      unsigned                                        minNSHits_;
      float                                           tbin_;
      unsigned                                        minTimeYbin_;
      float                                           maxTimeDT_;
      bool                                            filterMVA_;
      unsigned                                        minHitSelect_;
      float                                           minCutMVA_;
      bool                                            splitPhi_;
      bool                                            reqphi_;
      float                                           maxDeltaPhi_;
      int                                             maxNdiff_;
      int                                             diag_;
      std::unique_ptr<ModuleHistToolBase>             diagTool_;
      Data_types                                      data_;


      void findClusters           (const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                   TimeClusterCollection& tccol);
      void findTimePeaks          (const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                   TimePhiCandidateCollection& timeCandidates);
      void assignHits1            (const ComboHitCollection& chcol, const std::vector<unsigned>& chCood,
                                   float timePeakLow, float timePeakHigh, TimePhiCandidate& tc);
      void assignHits2            (const ComboHitCollection& chcol, const std::vector<unsigned>& chCood,
                                   const std::vector<float>& timePeaks, unsigned ipeak, TimePhiCandidate& tc);
      void addCalo                (const art::Handle<CaloClusterCollection>& ccH, std::vector<unsigned>& timeHist, float tmin);
      void addCaloPtr             (const art::Handle<CaloClusterCollection>& ccH,  TimePhiCandidate& tc, float time);
      void findPhiPeaks           (const ComboHitCollection& chcol, TimePhiCandidateCollection& timeCandidates,
                                   TimePhiCandidateCollection& phiCandidates);
      void filterMVACluster       (const ComboHitCollection& chcol, TimePhiCandidateCollection& candidates);
      void calculateMean          (const ComboHitCollection& chcol, TimeCluster& tc);
      void flagDuplicates         (const ComboHitCollection& chcol, TimePhiCandidateCollection& timeCandidates,
                                   TimePhiCandidateCollection& phiCandidates);
      void fillTCcol              (const TimePhiCandidateCollection& candidates,const ComboHitCollection& chcol,
                                         TimeClusterCollection& tccol);
      void fillDiag               (const TimePhiCandidateCollection& timeCandidates, const ComboHitCollection& chcol,
                                   const TimeClusterCollection& tccol);

      inline bool goodHit(const StrawHitFlag& flag) const {return flag.hasAllProperties(hsel_) && !flag.hasAnyProperty(hbkg_);}

  };




  TimeAndPhiClusterFinder::TimeAndPhiClusterFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    iev_(0),
    chToken_             {consumes<ComboHitCollection>      (config().comboHitCollection())     },
    ccToken_             {mayConsume<CaloClusterCollection> (config().caloClusterCollection())  },
    MVATime_             (config().MVATime()),
    hsel_                (config().hsel()),
    hbkg_                (config().hbkg()),
    testflag_            (config().testflag()),
    usecc_               (config().usecc()),
    ccmine_              (config().ccmine()),
    ccweight_            (config().ccweight()),
    algoAssignHits_      (config().algoAssignHits()),
    minNSHits_           (config().minNSHits()),
    tbin_                (config().tbin()),
    minTimeYbin_         (config().minTimeYbin()),
    maxTimeDT_           (config().maxTimeDT()),
    filterMVA_           (config().filterMVA()),
    minHitSelect_        (config().minHitSelect()),
    minCutMVA_           (config().minCutMVA()),
    splitPhi_            (config().splitPhi()),
    reqphi_              (config().reqphi()),
    maxDeltaPhi_         (config().maxDeltaPhi()),
    maxNdiff_            (config().maxNdiff()),
    diag_                (config().diag()),
    diagTool_(),
    data_()
    {
       produces<TimeClusterCollection>();

       if (algoAssignHits_ !=1 && algoAssignHits_ !=2)
           throw cet::exception("CATEGORY")<< "Unrecognized time Algorthm in TimeAndPhiClusterFinder module";

       if (diag_) diagTool_ = art::make_tool<ModuleHistToolBase>(config().diagPlugin," ");
    }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::beginJob()
  {
      MVATime_.initMVA();

      if (diag_){
         art::ServiceHandle<art::TFileService> tfs;
         diagTool_->bookHistograms(tfs);
         MVATime_.showMVA();
      }
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::produce(art::Event & event )
  {
      iev_ = event.id().event();

      art::Handle<CaloClusterCollection> ccH{};
      if (usecc_) ccH = event.getHandle<CaloClusterCollection>(ccToken_);

      const auto& chH = event.getValidHandle(chToken_);
      const auto& chcol(*chH);


      if (diag_) {data_.reset(); data_.event_=&event; data_.chcol_ = &chcol;}

      std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
      findClusters(ccH, chcol, *tccol);

      if (diag_) diagTool_->fillHistograms(&data_);
      event.put(std::move(tccol));
  }


  //----------------------------------------------------------ch----------------------------------------------------
  void TimeAndPhiClusterFinder::findClusters(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                             TimeClusterCollection& tccol)
  {
      // Find time peaks
      std::vector<TimePhiCandidate> timeCandidates;
      timeCandidates.reserve(64);
      findTimePeaks(ccH, chcol, timeCandidates);


      // Refine time peaks and create split phi clusters
      std::vector<TimePhiCandidate> phiCandidates;
      phiCandidates.reserve(64);
      if (splitPhi_) findPhiPeaks(chcol, timeCandidates, phiCandidates);


      // Apply MVA filtering for timeCandidates if needed
      if (filterMVA_) filterMVACluster(chcol,timeCandidates);


      //flag duplicated sequences so that we don't save them
      if (!reqphi_)
        flagDuplicates(chcol,timeCandidates,phiCandidates);


      // Finally create the timeClusters
      tccol.reserve(64);
      if (!reqphi_)
        fillTCcol(timeCandidates,chcol,tccol);
      fillTCcol(phiCandidates,chcol,tccol);

      if (diag_) fillDiag(timeCandidates, chcol, tccol);
  }


  //--------------------------------------------------------------------------------------------------------------
  // Find peaks in time distribution
  //--------------------------------------------------------------------------------------------------------------

  void TimeAndPhiClusterFinder::findTimePeaks(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                              TimePhiCandidateCollection& timeCandidates)
  {
      // Select good hits and sort them by time
      std::vector<unsigned> chGood;
      chGood.reserve(chcol.size());
      for (size_t ich=0; ich<chcol.size();++ich){
          if (!testflag_ || goodHit(chcol[ich].flag())) chGood.emplace_back(ich);
      }
      auto sortFcn = [&chcol](unsigned i1, unsigned i2){return chcol[i1].correctedTime() < chcol[i2].correctedTime();};
      sort(chGood.begin(),chGood.end(),sortFcn);
      if(chGood.empty()) return;

      // Fill the time histogram with hit corrected times (add buffer to hist boundaries). Add calo hits if requested
      float tmax     = chcol[chGood.back()].correctedTime() + 4*maxTimeDT_;
      float tmin     = std::max(0.0f, chcol[chGood.front()].correctedTime() - 4*maxTimeDT_);
      unsigned nbins = unsigned((tmax-tmin)/tbin_);

      std::vector<unsigned> timeHist(nbins,0);
      for (auto ich : chGood){
          unsigned ibin = unsigned((chcol[ich].correctedTime() - tmin)/tbin_);
          timeHist[ibin] += chcol[ich].nStrawHits();
      }
      if (usecc_) addCalo(ccH, timeHist, tmin);

      // Scan for local maxima (algorithm basd on maximum in sub-array with queue)
      std::vector<float> timePeaks;
      unsigned k(2u);
      std::deque<unsigned> q;
      for (unsigned i=0;i<timeHist.size();++i) {
           while (!q.empty() && timeHist[i] > timeHist[q.back()]) q.pop_back();
           while (!q.empty() && q.front()+k < i)                  q.pop_front();
           q.push_back(i);
           if (i<k) continue;

           const auto idx = q.front();
           if (idx+k/2 == i && timeHist[idx] >= minTimeYbin_) {
               float sumBins = timeHist[idx]+timeHist[idx-1]+timeHist[idx+1];
               float average = float(idx*timeHist[idx]+(idx-1)*timeHist[idx-1]+(idx+1)*timeHist[idx+1])/sumBins;
               timePeaks.push_back(tmin + average*tbin_);
           }
      }
      if (timePeaks.empty()) return;

      //Create time cluster candidates
      for (size_t i=0;i<timePeaks.size();++i){
           TimePhiCandidate tc;
           if (algoAssignHits_ == 1) assignHits1(chcol,chGood, timePeaks[i]- maxTimeDT_, timePeaks[i]+maxTimeDT_, tc);
           else                      assignHits2(chcol,chGood, timePeaks, i, tc);

           if (tc.nsh_ < minNSHits_) continue;

           if (usecc_) addCaloPtr(ccH, tc, timePeaks[i]);
           timeCandidates.emplace_back(std::move(tc));
      }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Assign hits to timeCandidates to maximize efficiency - hits can be assigned to several timeCandidates.
  // Note: hits must be time ordered
  void TimeAndPhiClusterFinder::assignHits1(const ComboHitCollection& chcol, const std::vector<unsigned>& chGood,
                                            float timePeakLow, float timePeakHigh, TimePhiCandidate& tc)
  {
      for (const auto& ich : chGood){
          float time = chcol[ich].correctedTime();
          if (time < timePeakLow) continue;
          if (time > timePeakHigh) break;
          tc.strawIdx_.emplace_back(ich);
          tc.nsh_ += chcol.at(ich).nStrawHits();
      }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Assign hits to time clusters to maximize purity - hits are assigned to closest timeCandidates
  // use the fact that TimeClusters are time ordered.  Note: hits must be time ordered
  void TimeAndPhiClusterFinder::assignHits2(const ComboHitCollection& chcol, const std::vector<unsigned>& chGood,
                                            const std::vector<float>& timePeaks, unsigned ipeak, TimePhiCandidate& tc)
  {
      for (const auto& ich : chGood){
          float time = chcol[ich].correctedTime();
          if (time < timePeaks[ipeak] - maxTimeDT_) continue;
          if (time > timePeaks[ipeak] + maxTimeDT_) break;

          // Check if the previous peak or next peak is a better match
          float dt = abs(chcol[ich].correctedTime()-timePeaks[ipeak]);
          float dt0 = (ipeak>0) ? abs(time-timePeaks[ipeak-1]) : 1e6;
          float dt1 = (ipeak+1<timePeaks.size()) ? abs(time-timePeaks[ipeak+1]) : 1e6;
          if (dt > dt0 || dt > dt1 || dt > maxTimeDT_) continue;

          tc.strawIdx_.emplace_back(ich);
          tc.nsh_ += chcol.at(ich).nStrawHits();
      }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::addCalo(const art::Handle<CaloClusterCollection>& ccH, std::vector<unsigned>& timeHist, float tmin){

      const CaloClusterCollection& cccol = *ccH.product();
      if (cccol.empty()) return;

      for (size_t icalo=0;icalo < cccol.size();++icalo){
          unsigned ibin = unsigned((cccol[icalo].time() - tmin)/tbin_);
          if (cccol[icalo].energyDep() > ccmine_ && ibin<timeHist.size()) timeHist[ibin] += ccweight_;
      }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::addCaloPtr(const art::Handle<CaloClusterCollection>& ccH, TimePhiCandidate& tc, float t0){

      const CaloClusterCollection& cccol = *ccH.product();
      if (cccol.empty()) return;

      auto caloFcn  = [t0](const CaloCluster& a, const CaloCluster& b) {return abs(a.time()-t0)<abs(b.time()-t0);};
      unsigned icalo = std::distance(cccol.begin(),std::min_element(cccol.begin(),cccol.end(),caloFcn));
      if (abs(cccol.at(icalo).time()-t0) < maxTimeDT_) tc.caloCluster_ = art::Ptr<CaloCluster>(ccH,icalo);
   }




  //--------------------------------------------------------------------------------------------------------------
  //Filtering and split clusters in phi
  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::findPhiPeaks(const ComboHitCollection& chcol, TimePhiCandidateCollection& timeCandidates,
                                             TimePhiCandidateCollection& phiCandidates)
  {
       // Cache phi value of comboHits for efficiency, this is recalculated every time the fcn is called
       std::vector<float> chPhi;
       chPhi.reserve(chcol.size());
       for (const auto& ch : chcol) chPhi.emplace_back(polyAtan2(ch.pos().y(),ch.pos().x()));

       for (auto& tc : timeCandidates){

           // Sort hits in ascending azimuthal angle
           auto& hits = tc.strawIdx_;
           sort(hits.begin(),hits.end(),[&chPhi](const int a, const int b){return chPhi[a]<chPhi[b];});
           if (hits.size()<minNSHits_) continue;

           // Start by fast forwarding to first large gap in phi.
           // note: (idx+vector_size)%vector_size maps negative and overflow indices -> valid vector indices
           size_t idx(0),vsize(hits.size());
           while (idx < vsize){
               float deltaPhi = chPhi[hits[idx]] - chPhi[hits[(idx+vsize-1)%vsize]];
               if (deltaPhi < 0) deltaPhi += 6.2832;
               if (deltaPhi > maxDeltaPhi_) break;
               ++idx;
           }

           // Find sequences of contiguous phi hits - first (last) point separated by an phi angle greater than threshold
           // with previous (next) point - flag them for removal if their number is below a threshold
           // Only save split clusters smaller than the original cluster to avoid duplicates
           size_t idxMax(idx+vsize), nStrawHits(0);
           std::vector<StrawHitIndex> buffer;
           while (idx < idxMax){
               buffer.emplace_back(hits[idx%vsize]);
               nStrawHits += chcol[hits[idx%vsize]].nStrawHits();

               float deltaPhi = chPhi[hits[(idx+1)%vsize]]-chPhi[hits[idx%vsize]];
               if (deltaPhi < 0) deltaPhi += 6.2832;

               if (deltaPhi< maxDeltaPhi_ && idx+1!=idxMax ){++idx;continue;}
               if (nStrawHits >= minNSHits_) {
                   phiCandidates.emplace_back(TimePhiCandidate(nStrawHits,buffer,tc.caloCluster_));
               }

               buffer.clear();
               nStrawHits = 0;
               ++idx;
           }
       }
   }

  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::filterMVACluster(const ComboHitCollection& chcol, TimePhiCandidateCollection& candidates)
  {
      TimePhiCandidateCollection tempCand;
      std:: vector<float> parsMVA(3,0.0);

      for (auto& cand : candidates){
          auto& hits = cand.strawIdx_;

          // calculate mean values of the cluster
          float mean_r(0),mean_t(0),mean_x(0),mean_y(0),sweight(0);
          for (const auto& hit : hits) {
              const ComboHit& ch = chcol[hit];
              mean_r  += sqrt(ch.pos().x()*ch.pos().x()+ch.pos().y()*ch.pos().y())*ch.nStrawHits();
              mean_x  += ch.pos().x()*ch.nStrawHits();
              mean_y  += ch.pos().y()*ch.nStrawHits();
              mean_t  += ch.correctedTime()*ch.nStrawHits();
              sweight += ch.nStrawHits();
          }
          mean_r /= sweight;
          mean_t /= sweight;
          mean_x /= sweight;
          mean_y /= sweight;
          float mean_p = polyAtan2(mean_y,mean_x);

          if (sweight > minHitSelect_) tempCand.emplace_back(cand);

          // loop ovr hits, calculate MVA output and flag bad hits
          for (auto& hit : hits){
              const ComboHit& ch = chcol[hit];

              float dphi = polyAtan2(ch.pos().y(),ch.pos().x()) - mean_p;
              if (dphi > 3.14159)  dphi -= 6.2832;
              if (dphi < -3.14159) dphi += 6.2832;

              parsMVA[0] = sqrt(ch.pos().x()*ch.pos().x()+ch.pos().y()*ch.pos().y())-mean_r;
              parsMVA[1] = dphi;
              parsMVA[2] = ch.correctedTime()-mean_t;

              //potential speed improvement: if dt, dr and dphi are all "small", then no need to run the MVA.

              Float_t outMVA =  MVATime_.evalMVA(parsMVA);
              if (outMVA < minCutMVA_) {cand.nsh_ -= ch.nStrawHits(); hit = chcol.size()+1;}
          }

          hits.erase(std::remove_if(hits.begin(), hits.end(), [&chcol](auto ich){return ich>chcol.size();}),hits.end());
      }
      std::move(tempCand.begin(), tempCand.end(), std::back_inserter(candidates));
  }



  //--------------------------------------------------------------------------------------------------------------
  // Caluclate cluster variables, find duplicates and fill collections
  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::calculateMean(const ComboHitCollection& chcol, TimeCluster& tc)
  {
      if (tc._strawHitIdxs.empty()){tc._pos = XYZVectorF(0,0,0);tc._t0._t0=0;tc._t0._t0err=0;tc._nsh=0; return;};

      tc._nsh = 0;
      float tacc(0),tacc2(0),xacc(0),yacc(0),zacc(0),weight(0);
      for (const auto& ish :tc._strawHitIdxs){
          const ComboHit& ch = chcol[ish];
          tc._nsh += ch.nStrawHits();

          float htime = chcol[ish].correctedTime();
          float hwt   = ch.nStrawHits();

          weight += hwt;
          tacc   += htime*hwt;
          tacc2  += htime*htime*hwt;
          xacc   += ch.pos().x()*hwt;
          yacc   += ch.pos().y()*hwt;
          zacc   += ch.pos().z()*hwt;
      }

      tacc/=weight;
      tacc2/=weight;
      xacc/=weight;
      yacc/=weight;
      zacc/=weight;

      tc._t0._t0    = tacc;
      tc._t0._t0err = sqrtf(tacc2-tacc*tacc);
      tc._pos       = XYZVectorF(xacc, yacc, zacc);
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::flagDuplicates(const ComboHitCollection& chcol, TimePhiCandidateCollection& timeCandidates,
                                               TimePhiCandidateCollection& phiCandidates)
  {
      for (auto& tc : timeCandidates) sort(tc.strawIdx_.begin(),tc.strawIdx_.end());
      for (auto& pc : phiCandidates)  sort(pc.strawIdx_.begin(),pc.strawIdx_.end());

      for (auto& tc : timeCandidates){
          for (auto& pc : phiCandidates){
              if (pc.nsh_ == 0) continue;
              if (tc.strawIdx_.size()+maxNdiff_ < pc.strawIdx_.size() ||
                  pc.strawIdx_.size()+maxNdiff_ < tc.strawIdx_.size()) continue;

              int nMatch(0),i1(0),i2(0),isize1(tc.strawIdx_.size()),isize2(pc.strawIdx_.size());
              while (i1 < isize1 && i2 < isize2) {
                  if      (tc.strawIdx_[i1] < pc.strawIdx_[i2]) ++i1;
                  else if (tc.strawIdx_[i1] > pc.strawIdx_[i2]) ++i2;
                  else    {++nMatch; ++i1; ++i2;}
              }

              if (abs(nMatch-isize1)<= maxNdiff_ || abs(nMatch-isize2)<= maxNdiff_){
                  pc.nsh_=0;
                  pc.strawIdx_.clear();
              }
          }
      }
  }

  //----------------------------------------------------------ch----------------------------------------------------
  void TimeAndPhiClusterFinder::fillTCcol(const TimePhiCandidateCollection& candidates, const ComboHitCollection& chcol,
                                          TimeClusterCollection& tccol)
  {
      for (const auto& cand : candidates){
          if (cand.nsh_ < minNSHits_) continue;
          TimeCluster tclu;
          tclu._caloCluster  = cand.caloCluster_;
          tclu._nsh          = cand.nsh_;
          tclu._strawHitIdxs = std::move(cand.strawIdx_);

          calculateMean(chcol, tclu);
          tccol.emplace_back(std::move(tclu));
      }
  }



  //--------------------------------------------------------------------------------------------------------------
  // Diagnosis
  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::fillDiag(const TimePhiCandidateCollection& timeCandidates, const ComboHitCollection& chcol,
                                         const TimeClusterCollection& tccol )
  {
      data_.iev_ = iev_;

      data_.Nch_ = chcol.size();
      for (unsigned ich=0; ich<chcol.size();++ich)
      {
         int   selFlag = (!testflag_ || goodHit(chcol[ich].flag())) ? 1 : 0;
         //const StrawHitFlag ener("EnergySelection");
         //if (chcol[ich].flag().hasAllProperties(ener) && !chcol[ich].flag().hasAnyProperty(hbkg_)) selFlag = 2;

         data_.chSel_[ich]  = selFlag;
         data_.chTime_[ich] = chcol[ich].correctedTime();

         data_.chX_[ich]    = chcol[ich].pos().x();
         data_.chY_[ich]    = chcol[ich].pos().y();
         data_.chZ_[ich]    = chcol[ich].pos().z();
         data_.chRad_[ich]  = sqrt(data_.chX_[ich]*data_.chX_[ich]+data_.chY_[ich]*data_.chY_[ich]);
         data_.chPhi_[ich]  = chcol[ich].pos().phi();
         data_.chNhit_[ich] = chcol[ich].nStrawHits();
      }


      int iclu1(0),ih1(0);
      for (const auto& tc : timeCandidates)
      {
        for (auto& ish : tc.strawIdx_)
        {
          data_.nclu1_[ih1]  = iclu1;
          data_.hitIdx1_[ih1]= ish;
          ++ih1;
        }
        ++iclu1;
      }
      data_.nhit1_=ih1;


      int iclu2(0),ih2(0);
      for (const auto& tc : tccol)
      {
        for (auto& ish : tc._strawHitIdxs)
        {
          data_.nclu2_[ih2]    = iclu2;
          data_.hitIdx2_[ih2]  = ish;
          ++ih2;
        }
        ++iclu2;
      }
      data_.nhit2_=ih2;
   }
}


DEFINE_ART_MODULE(mu2e::TimeAndPhiClusterFinder)

