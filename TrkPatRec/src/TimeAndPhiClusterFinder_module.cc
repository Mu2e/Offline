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
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrkPatRec/inc/TimeAndPhiClusterFinder_types.hh"


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
        fhicl::Atom<art::InputTag>              strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection {Name") };
        fhicl::Atom<art::InputTag>              caloClusterCollection  {Name("CaloClusterCollection"),  Comment("Calo cluster collection {Name") };
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
        fhicl::Atom<bool>                       refinePhi              {Name("RefinePhi"),              Comment("Refine time cluster with phi info") };
        fhicl::Atom<float>                      maxDeltaPhiRef         {Name("MaxDeltaPhiRef"),         Comment("Max delta phi between consecutive hits for refine") };
        fhicl::Atom<float>                      maxNhitRef             {Name("MaxNhitRef"),             Comment("Max number of hits for refiningtime clusters") };
        fhicl::Atom<bool>                       splitPhi               {Name("SplitPhi"),               Comment("Split time cluster with phi info") };
        fhicl::Atom<float>                      maxDeltaPhi            {Name("MaxDeltaPhi"),            Comment("Max delta phi between consecutive hits in cluster") };
        fhicl::Atom<bool>                       filterTime             {Name("FilterTime"),             Comment("Final time filtering of clusters") };
        fhicl::Atom<float>                      deltaTFilter           {Name("DeltaTimeFilter"),        Comment("Max time difference between mean time and hit time in cluster") };
        fhicl::Atom<int>                        diag                   {Name("Diag"),                   Comment("Diag level"), 0 };
        fhicl::Table<Config_types>              diagPlugin             {Name("DiagPlugin"),             Comment("Diag Plugin config")};
      };


      explicit TimeAndPhiClusterFinder(const art::EDProducer::Table<Config>& config);
      ~TimeAndPhiClusterFinder() = default;

      void     beginJob() override;
      void     produce(art::Event& e) override;


    private:
      int                                             iev_;
      const art::ProductToken<ComboHitCollection>     chToken_;
      const art::ProductToken<StrawHitFlagCollection> shfToken_;
      const art::ProductToken<CaloClusterCollection>  ccToken_;
      StrawHitFlag                                    hsel_, hbkg_;
      bool                                            testflag_;
      bool                                            usecc_;
      float                                           ccmine_,ccweight_;
      unsigned                                        algoAssignHits_;
      unsigned                                        minNSHits_;
      float                                           tbin_;
      unsigned                                        minTimeYbin_;
      float                                           maxTimeDT_;
      bool                                            refinePhi_;
      float                                           maxDeltaPhiRef_;
      unsigned                                        maxNhitRef_;
      bool                                            splitPhi_;
      float                                           maxDeltaPhi_;
      bool                                            filterTime_;
      float                                           deltaTFilter;
      int                                             diag_;
      std::unique_ptr<ModuleHistToolBase>             diagTool_;
      Data_types                                      data_;


      void findClusters           (const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                   const StrawHitFlagCollection& shfcol,          TimeClusterCollection& tccol);
      void findTimePeaks          (const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                   const StrawHitFlagCollection& shfcol,          TimeClusterCollection& tccol);
      void assignHits1            (const ComboHitCollection& chcol, const std::vector<unsigned>& chCood,
                                   const float timeDown, const float TimeUp, TimeCluster& tc);
      void assignHits2            (const ComboHitCollection& chcol, const std::vector<unsigned>& chCood,
                                   const std::vector<float>& timePeaks, unsigned ipeak, TimeCluster& tc);
      void addCalo                (const art::Handle<CaloClusterCollection>& ccH, std::vector<unsigned>& timeHist, float tmin);
      void addCaloPtr             (const art::Handle<CaloClusterCollection>& ccH,  TimeCluster& tc, float time);
      void findPhiPeaks           (const ComboHitCollection& chcol, TimeClusterCollection& tccol, TimeClusterCollection& tccolPhi);
      void refinePhi              (const ComboHitCollection& chcol, const std::vector<float>& chPhi, TimeCluster& tc);
      void splitPhi               (const ComboHitCollection& chcol, const std::vector<float>& chPhi, TimeCluster& tc,
                                   TimeClusterCollection& tccolPhi);
      void filterTime             (const ComboHitCollection& chcol, TimeCluster& pc);
      void calculateMean          (const ComboHitCollection& chcol, TimeCluster& tc);
//      void fillDiag               (const TimePhiCandidateCollection& timeCandidates, const ComboHitCollection& chcol,
//                                   const StrawHitFlagCollection& shfcol,             const TimeClusterCollection& tccol);

      inline bool goodHit(const StrawHitFlag& flag) const {return flag.hasAllProperties(hsel_) && !flag.hasAnyProperty(hbkg_);}

  };




  TimeAndPhiClusterFinder::TimeAndPhiClusterFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    iev_(0),
    chToken_             {consumes<ComboHitCollection>      (config().comboHitCollection())     },
    shfToken_            {mayConsume<StrawHitFlagCollection>(config().strawHitFlagCollection()) },
    ccToken_             {mayConsume<CaloClusterCollection> (config().caloClusterCollection())  },
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
    refinePhi_           (config().refinePhi()),
    maxDeltaPhiRef_      (config().maxDeltaPhiRef()),
    maxNhitRef_          (config().maxNhitRef()),
    splitPhi_            (config().splitPhi()),
    maxDeltaPhi_         (config().maxDeltaPhi()),
    filterTime_          (config().filterTime()),
    deltaTFilter         (config().deltaTFilter()),
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
      if (diag_){
         art::ServiceHandle<art::TFileService> tfs;
         diagTool_->bookHistograms(tfs);
      }
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::produce(art::Event & event )
  {
      iev_ = event.id().event();

      art::Handle<CaloClusterCollection> ccH{};
      if (usecc_){ccH = event.getHandle<CaloClusterCollection>(ccToken_);}

      const auto& shfH = event.getValidHandle(shfToken_);
      const auto& chH = event.getValidHandle(chToken_);
      const auto& shfcol(*shfH);
      const auto& chcol(*chH);

      if (testflag_ && shfcol.size() != chcol.size())
        throw cet::exception("RECO")<<"TimeAndPhiClusterFinder: inconsistent flag collection length " << std::endl;

      if (diag_){data_.reset(); data_.event_=&event; data_.chcol_ = &chcol;}

      std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
      findClusters(ccH, chcol, shfcol, *tccol);

      if (diag_) diagTool_->fillHistograms(&data_);
      event.put(std::move(tccol));
  }


  //----------------------------------------------------------ch----------------------------------------------------
  void TimeAndPhiClusterFinder::findClusters(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                             const StrawHitFlagCollection& shfcol,          TimeClusterCollection& tccol)
  {
      // Find time peaks
      findTimePeaks(ccH, chcol, shfcol, tccol);

      // Refine time peaks and create split phi clusters
      TimeClusterCollection tccolPhi;
      if (refinePhi_ || splitPhi_) findPhiPeaks(chcol, tccol, tccolPhi);

      // Refine time of the split clusters
      if (filterTime_) {for (auto& tcp : tccolPhi) filterTime(chcol, tcp);}

      // If filtering was done, check if the candidates still have the required number of hits
      //auto predTC = [this](const TimeCluster& a){return a.nsh_ < minNSHits_;};
      //if (refinePhi_)  tccol.erase(std::remove_if(tccol.begin(),tccol.end(),predTC),tccol.end()) ;
      //if (filterTime_) tccolPhi.erase( std::remove_if(tccolPhi.begin(),tccolPhi.end(),predTC),tccolPhi.end()) ;

      std::move(tccolPhi.begin(), tccolPhi.end(), std::back_inserter(tccol));
      for (auto& tc : tccol) calculateMean(chcol, tc);

      //if (diag_) fillDiag(timeCandidates, chcol, shfcol, tccol);
  }


  //--------------------------------------------------------------------------------------------------------------
  // Find peaks in time distribution
  //--------------------------------------------------------------------------------------------------------------

  void TimeAndPhiClusterFinder::findTimePeaks(const art::Handle<CaloClusterCollection>& ccH, const ComboHitCollection& chcol,
                                              const StrawHitFlagCollection& shfcol, TimeClusterCollection& tccol)
  {
      // Select good hits and sort them by time
      std::vector<unsigned> chGood;
      chGood.reserve(chcol.size());
      for (size_t ich=0; ich<chcol.size();++ich){
          if (!testflag_ || goodHit(shfcol[ich])) chGood.emplace_back(ich);
      }
      auto predGood = [&chcol](unsigned i1, unsigned i2){return chcol[i1].correctedTime() < chcol[i2].correctedTime();};
      sort(chGood.begin(),chGood.end(),predGood);


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
      std::vector<float> timePeaks, timePeaksUp, timePeaksDown;
      unsigned k(2u);
      std::deque<unsigned> q;
      for (unsigned i=0;i<timeHist.size();++i) {
           while (!q.empty() && timeHist[i] > timeHist[q.back()]) q.pop_back();
           while (!q.empty() && q.front()+k < i)                  q.pop_front();
           q.push_back(i);

           if (i<k) continue;

           unsigned idx = q.front();
           if (idx+k/2 == i && timeHist[idx] >= minTimeYbin_){

               //Scan to find the boundaries of the timecluster (until there are zero hits)
               unsigned iup(idx);   while (iup < timeHist.size() && timeHist[iup]>0) ++iup;
               unsigned idown(idx); while (idown>0 && timeHist[idown]>0) --idown;

               //Average over two neighboring bins to get the central time and cap timeUp/Down
               float sumBins  = timeHist[idx]+timeHist[idx-1]+timeHist[idx+1];
               float average  = float(idx*timeHist[idx]+(idx-1)*timeHist[idx-1]+(idx+1)*timeHist[idx+1]) / sumBins;
               float time     = tmin + average*tbin_;
               float timeDown = tmin + idown*tbin_ + tbin_;
               float timeUp   = tmin + iup*tbin_;

               if (time-timeDown > maxTimeDT_) timeDown = time-maxTimeDT_;
               if (timeUp-time > maxTimeDT_)   timeUp   = time+maxTimeDT_;

               timePeaks.push_back(time);
               timePeaksDown.push_back(timeDown);
               timePeaksUp.push_back(timeUp);
           }
      }
      if (timePeaks.empty()) return;

      //Create time cluster candidates
      for (size_t i=0;i<timePeaks.size();++i){
          TimeCluster tc;
          if (algoAssignHits_ == 1) assignHits1(chcol,chGood, timePeaksDown[i],timePeaksUp[i],tc);
          else                      assignHits2(chcol,chGood, timePeaks, i, tc);

          if (tc._nsh < minNSHits_) continue;

          if (usecc_) addCaloPtr(ccH, tc, timePeaks[i]);
          tccol.emplace_back(std::move(tc));
      }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Assign hits to timeCandidates to maximize efficiency - hits can be assigned to several timeCandidates
  void TimeAndPhiClusterFinder::assignHits1(const ComboHitCollection& chcol, const std::vector<unsigned>& chGood,
                                            const float timeDown, const float timeUp, TimeCluster& tc)
  {
      for (const auto& ich : chGood){
          if (chcol[ich].correctedTime() < timeDown) continue;
          if (chcol[ich].correctedTime() > timeUp)   break;
          tc._strawHitIdxs.emplace_back(ich);
          tc._nsh += chcol.at(ich).nStrawHits();
      }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Assign hits to time clusters to maximize purity - hits are assigned to closest timeCandidates
  // use the fact that TimeClusters are time ordered
  void TimeAndPhiClusterFinder::assignHits2(const ComboHitCollection& chcol, const std::vector<unsigned>& chGood,
                                            const std::vector<float>& timePeaks, unsigned ipeak, TimeCluster& tc)
  {
      for (const auto& ich : chGood){
          if (chcol[ich].correctedTime() < timePeaks[ipeak] - maxTimeDT_) continue;
          if (chcol[ich].correctedTime() > timePeaks[ipeak] + maxTimeDT_) break;

          // Check if the previous peak or next peak is a better match
          float dt = abs(chcol[ich].correctedTime()-timePeaks[ipeak]);
          float dt0 = (ipeak>0) ? abs(chcol[ich].correctedTime()-timePeaks[ipeak-1]) : 1e6;
          float dt1 = (ipeak+1<timePeaks.size()) ? abs(chcol[ich].correctedTime()-timePeaks[ipeak+1]) : 1e6;
          if (dt > dt0 || dt > dt1 || dt > maxTimeDT_) continue;

          tc._strawHitIdxs.emplace_back(ich);
          tc._nsh += chcol.at(ich).nStrawHits();
      }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::addCalo(const art::Handle<CaloClusterCollection>& ccH, std::vector<unsigned>& timeHist, float tmin){

      const CaloClusterCollection& cccol = *ccH.product();
      if (cccol.empty()) return;

      for (size_t icalo=0;icalo < cccol.size();++icalo){
          unsigned ibin = unsigned((cccol[icalo].time() - tmin)/tbin_);
          if (cccol[icalo].energyDep() > ccmine_) timeHist[ibin] += ccweight_;
      }
  }

  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::addCaloPtr(const art::Handle<CaloClusterCollection>& ccH, TimeCluster& tc, float t0){

      const CaloClusterCollection& cccol = *ccH.product();
      if (cccol.empty()) return;

      auto predCalo = [t0](const CaloCluster& a, const CaloCluster& b) {return abs(a.time()-t0)<abs(b.time()-t0);};
      int icalo     = std::distance(cccol.begin(),std::min_element(cccol.begin(),cccol.end(),predCalo));
      if (abs(cccol.at(icalo).time()-t0) < maxTimeDT_) tc._caloCluster = art::Ptr<CaloCluster>(ccH,icalo);
   }




  //--------------------------------------------------------------------------------------------------------------
  // Refine the timeCandidate  and split clusters in phi
  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::findPhiPeaks(const ComboHitCollection& chcol, TimeClusterCollection& tccol, TimeClusterCollection& tccolPhi)
  {
       // Cache phi value of comboHits for efficiency, this is recalculated every time the fcn is called
       std::vector<float> chphi;
       chphi.reserve(chcol.size());
       for (const auto& ch : chcol) chphi.emplace_back(ch.phi());

       // Refine: remove isolared, small cluster of hits
       // Split: break time clusters in sub-clusters in azimuthal plane
       for (auto& tc : tccol){

           // Sort hits in ascending azimuthal angle
           auto& hits = tc._strawHitIdxs;
           std::sort(hits.begin(),hits.end(),[&chphi](const int a, const int b){return chphi[a]<chphi[b];});

           if (refinePhi_) refinePhi(chcol, chphi, tc);
           if (splitPhi_)  splitPhi( chcol, chphi, tc, tccolPhi);
       }
   }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::refinePhi(const ComboHitCollection& chcol, const std::vector<float>& chPhi, TimeCluster& tc)
  {
        auto& hits = tc._strawHitIdxs;

        // Start by fast forwarding to first large gap in phi.
        // note: (idx+vector_size)%vector_size maps negative and overflow indices -> valid vector indices  -- neat!
        size_t idx(0),vsize(hits.size());
        while (idx < vsize){
            float deltaPhi = chPhi[hits[idx]] - chPhi[hits[(idx+vsize-1)%vsize]];
            if (deltaPhi < 0) deltaPhi += 6.2832; //deal with discontinuities pi -> -pi
            if (deltaPhi > maxDeltaPhiRef_) break;
            ++idx;
        }

        // Find sequences of contiguous phi hits - first (last) point separated by an phi angle greater than threshold
        // with previous (next) point - flag them for removal if their number is below a threshold
        size_t idxMax(idx+vsize);
        std::vector<StrawHitIndex> flaggedHits;
        while (idx < idxMax){
            flaggedHits.emplace_back(idx%vsize);

            float deltaPhi = chPhi[hits[(idx+1)%vsize]]-chPhi[hits[idx%vsize]];
            if (deltaPhi < 0) deltaPhi += 6.2832;

            if (deltaPhi< maxDeltaPhiRef_ && idx+1!=idxMax ){++idx;continue;}

            if (flaggedHits.size() <= maxNhitRef_) {
                for (auto iflag : flaggedHits) {
                    tc._nsh    -= chcol[hits[iflag]].nStrawHits();
                    hits[iflag] = chPhi.size()+1;
                }
            }
            flaggedHits.clear();
            ++idx;
        }
        hits.erase(std::remove_if(hits.begin(), hits.end(),[&chPhi](unsigned i) {return i>chPhi.size();}), hits.end());
   }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::splitPhi(const ComboHitCollection& chcol, const std::vector<float>& chPhi,
                                         TimeCluster& tc, TimeClusterCollection& tccolPhi)
  {
        const auto& hits = tc._strawHitIdxs;
        if (hits.size()<minNSHits_) return;

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
            if (nStrawHits >= minNSHits_ && buffer.size() < hits.size()){

                TimeCluster pc;
                pc._strawHitIdxs = buffer;
                pc._nsh = nStrawHits;
                pc._caloCluster = tc._caloCluster;

                tccolPhi.emplace_back(std::move(pc));
            }
            buffer.clear();
            nStrawHits = 0;
            ++idx;
        }
   }


   //--------------------------------------------------------------------------------------------------------------
   // Calulate mean time of cluster and remove hits too far away (keep correct number of strawHits in cluster)
   void TimeAndPhiClusterFinder::filterTime(const ComboHitCollection& chcol, TimeCluster& tc)
   {
       auto& hits = tc._strawHitIdxs;

       float tmean(0),weight(0);
       for (const auto& ich : hits){
           const ComboHit& ch = chcol[ich];
           weight += ch.nStrawHits();
           tmean  += ch.correctedTime()*ch.nStrawHits();
       }
       tmean /= weight;

       for (auto& hit : hits){
           if (std::abs(chcol[hit].correctedTime()-tmean) < deltaTFilter) continue;
           tc._nsh -= chcol[hit].nStrawHits();
           hit = chcol.size()+1;
       }
       hits.erase(std::remove_if(hits.begin(), hits.end(), [&chcol](auto ich){return ich>chcol.size();}),hits.end());
   }



  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::calculateMean(const ComboHitCollection& chcol, TimeCluster& tc)
  {
      if (tc._strawHitIdxs.empty()){tc._pos = XYZVectorF(0,0,0);tc._t0._t0=0; tc._t0._t0err =0; tc._nsh = 0; return;};

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
  // Diagnosis
  //--------------------------------------------------------------------------------------------------------------
  /*
  void TimeAndPhiClusterFinder::fillDiag(const TimePhiCandidateCollection& timeCandidates, const ComboHitCollection& chcol,
                                         const StrawHitFlagCollection& shfcol,             const TimeClusterCollection& tccol )
  {
      data_.iev_ = iev_;

      data_.Nch_ = chcol.size();
      for (unsigned ich=0; ich<chcol.size();++ich)
      {
         int   selFlag = (!testflag_ || goodHit(shfcol[ich])) ? 1 : 0;
         const StrawHitFlag ener("EnergySelection");
         if (shfcol[ich].hasAllProperties(ener) && !shfcol[ich].hasAnyProperty(hbkg_)) selFlag = 2;

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
   */
}


DEFINE_ART_MODULE(mu2e::TimeAndPhiClusterFinder);

