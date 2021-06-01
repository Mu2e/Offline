//
// Tracker time / phi cluster finder
//
// Original author B. Echenard
//

#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "GeneralUtilities/inc/Angles.hh"
#include "Mu2eUtilities/inc/polyAtan2.hh"
#include "Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Mu2eUtilities/inc/MVATools.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "TrkReco/inc/TrkUtilities.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "TrkPatRec/inc/TimeAndPhiClusterFinder_types.hh"

#include <algorithm>
#include <numeric>


namespace
{ 
   struct TimeCandidate
   {
      TimeCandidate()  {strawIdx_.reserve(32);} 

      float                       t0_       = 0;
      unsigned                    nsh_      = 0;
      std::vector<StrawHitIndex>  strawIdx_ = {};
      art::Ptr<mu2e::CaloCluster> caloCluster_;     
   };  

   typedef std::vector<TimeCandidate> TimeCandidateCollection;
}




namespace mu2e {
  
  using TimeAndPhiClusterFinderTypes::Data_t;

  class TimeAndPhiClusterFinder : public art::EDProducer 
  {     
     public:
       
       struct Config
       {
           using Name    = fhicl::Name;
           using Comment = fhicl::Comment;
           fhicl::Atom<art::InputTag>              comboHitCollection     {Name("ComboHitCollection"),     Comment("ComboHit collection {Name") };
           fhicl::Atom<art::InputTag>              strawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("StrawHitFlag collection {Name") };
           fhicl::Atom<art::InputTag>              caloClusterCollection  {Name("CaloClusterCollection"),  Comment("Calo cluster collection {Name") };
           fhicl::Table<MVATools::Config>          tcMVA                  {Name("TimeClusterMVA"),         Comment("MVA config for TimeAndPhi clutering") }; 
           fhicl::Sequence<std::string>            hsel                   {Name("HitSelectionBits"),       Comment("HitSelectionBits") }; 
           fhicl::Sequence<std::string>            hbkg                   {Name("HitBackgroundBits"),      Comment("HitBackgroundBits") }; 
           fhicl::Atom<bool>                       testflag               {Name("TestFlag"),               Comment("Test hit flags") }; 
           fhicl::Atom<float>                      tmin                   {Name("Tmin"),                   Comment("Time histogram start") }; 
           fhicl::Atom<float>                      tmax                   {Name("Tmax"),                   Comment("Time histogram end") }; 
           fhicl::Atom<float>                      tbin                   {Name("Tbin"),                   Comment("Time histogram bin width") }; 
           fhicl::Atom<float>                      pbin                   {Name("Pbin"),                   Comment("Phi search algorithm bin width") }; 
           fhicl::Atom<unsigned>                   minNSHits              {Name("MinNSHits"),              Comment("Minimum number of hits for cluster") }; 
           fhicl::Atom<bool>                       usecc                  {Name("UseCaloCluster"),         Comment("Use calorimeter cluster") }; 
           fhicl::Atom<float>                      ccmine                 {Name("CaloClusterMinE"),        Comment("Minimum energy for calorimeter cluster") }; 
           fhicl::Atom<float>                      ccweight               {Name("CaloClusterWeight"),      Comment("Weight of cluster in tracker hits") }; 
           fhicl::Table<TrkTimeCalculator::Config> ttcalc                 {Name("T0Calculator"),           Comment("TimeTracker calculator config") }; 
           fhicl::Atom<float>                      pitch                  {Name("AveragePitch"),           Comment("Average track pitch") }; 
           fhicl::Atom<unsigned>                   timeStrategy           {Name("TimeStrategy"),           Comment("Strategy to search for time peaks") }; 
           fhicl::Atom<unsigned>                   windowTime             {Name("WindowTime"),             Comment("Window width for time scan") }; 
           fhicl::Atom<unsigned>                   minTimeYbin            {Name("MinTimeYbin"),            Comment("Minimum number of bins to start recording max for scanning algo") }; 
           fhicl::Atom<float>                      maxTimeDT              {Name("MaxTimeDT"),              Comment("Max time difference for hits in cluster") }; 
           fhicl::Atom<float>                      maxTimeDTCal           {Name("MaxTimeDTCal"),           Comment("Max calo time difference for hits in cluster") }; 
           fhicl::Atom<unsigned>                   minHitSplit            {Name("MinHitSplit"),            Comment("Minimum numnber of hits to split phi cluster") }; 
           fhicl::Atom<float>                      MVAScoreCut            {Name("MVAScoreCut"),            Comment("min MVA score to remove hit") }; 
           fhicl::Atom<bool>                       recover                {Name("Recover"),                Comment("Apply hit recovery algorithm") }; 
           fhicl::Atom<float>                      dphiMaxReco            {Name("DphiMaxReco"),            Comment("Max phi difference for recovering hits") }; 
           fhicl::Atom<int>                        printfreq              {Name("Printfreq"),              Comment("Print frequency"), 100 }; 
           fhicl::Atom<int>                        diag                   {Name("Diag"),                   Comment("Diag level"), 0 }; 
           fhicl::Table<TimeAndPhiClusterFinderTypes::Config> diagPlugin  {Name("DiagPlugin"),             Comment("Diag Plugin config")}; 
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
       const StrawHitFlagCollection*                   shfcol_;
       const ComboHitCollection*                       chcol_;
       MVATools                                        tcMVA_; 
       StrawHitFlag                                    hsel_, hbkg_;
       bool		                               testflag_;
       float                                           tmin_, tmax_, tbin_, pbin_;
       unsigned                                        minNSHits_;
       bool                                            usecc_;
       float                                           ccmine_,ccweight_;
       TrkTimeCalculator                               ttcalc_;
       float		                               pitch_; 
       unsigned                                        timeStrategy_;
       unsigned                                        windowTime_;
       unsigned                                        minTimeYbin_;
       float                                           maxTimeDT_;  
       float                                           maxTimeDTCal_;  
       unsigned                                        minHitSplit_;
       float                                           MVAScoreCut_;
       bool                                            recover_;
       float                                           dphiMaxReco_;
       int                                             printfreq_, diag_;      
       std::unique_ptr<ModuleHistToolBase>             diagTool_;
       Data_t                                          data_;


       void findClusters           (const art::Handle<CaloClusterCollection>& ccH, TimeClusterCollection& tccol);
       void findCaloSeeds          (const art::Handle<CaloClusterCollection>& ccH, TimeCandidateCollection& tccol);
       void findTimePeaks          (const art::Handle<CaloClusterCollection>& ccH, TimeCandidateCollection& seeds,const std::vector<float>& comboHitTime);
       void timePeakScan           (const std::vector<unsigned>& timeHist, std::vector<unsigned>& timePeaks);
       void timePeakBlind          (const std::vector<unsigned>& timeHist, std::vector<unsigned>& timePeaks, const unsigned nbins);
       void findPhiPeaks           (TimeCandidate& tc, TimeClusterCollection& tccol); 
       void associateScan          (const std::vector<unsigned>& phiHist, std::vector<unsigned>& istart, std::vector<unsigned>& iend);
       void filterClusterSequential(TimeCluster& tc, const std::vector<float>& comboHitTime);
       void recoverHits            (TimeCluster& tc, const std::vector<float>& comboHitTime, const std::vector<std::vector<unsigned>>& hitIndex);       
       void calculateMean          (TimeCluster& tc, const std::vector<float>& comboHitTime);       
       void fillDiagPreClean       (const TimeCandidateCollection& timeCandidates, const TimeClusterCollection& tccol);
       void fillDiagPostClean      (const TimeClusterCollection& tccol);
       void fillDiagMVA            (unsigned iworseHit, const std::vector<float>& mva);

       inline bool goodHit(const StrawHitFlag& flag) const {return flag.hasAllProperties(hsel_) && !flag.hasAnyProperty(hbkg_);}
       
  };

  TimeAndPhiClusterFinder::TimeAndPhiClusterFinder(const art::EDProducer::Table<Config>& config) :
     art::EDProducer{config},
     iev_(0),
     chToken_      {consumes<ComboHitCollection>      (config().comboHitCollection()) },
     shfToken_     {mayConsume<StrawHitFlagCollection>(config().strawHitFlagCollection()) },
     ccToken_      {mayConsume<CaloClusterCollection> (config().caloClusterCollection()) },
     tcMVA_        (config().tcMVA()),
     hsel_         (config().hsel()),
     hbkg_         (config().hbkg()),
     testflag_     (config().testflag()),
     tmin_         (config().tmin()),
     tmax_         (config().tmax()),
     tbin_         (config().tbin()),
     pbin_         (config().pbin()),
     minNSHits_    (config().minNSHits()),
     usecc_        (config().usecc()),
     ccmine_       (config().ccmine()),
     ccweight_     (config().ccweight()),
     ttcalc_       (config().ttcalc()),
     pitch_        (config().pitch()),
     timeStrategy_ (config().timeStrategy()),
     windowTime_   (config().windowTime()),
     minTimeYbin_  (config().minTimeYbin()),
     maxTimeDT_    (config().maxTimeDT()),
     maxTimeDTCal_ (config().maxTimeDTCal()),
     minHitSplit_  (config().minHitSplit()),
     MVAScoreCut_  (config().MVAScoreCut()),
     recover_      (config().recover()),
     dphiMaxReco_  (config().dphiMaxReco()),
     printfreq_    (config().printfreq()),
     diag_         (config().diag()),
     diagTool_(),
     data_()
  {
       produces<TimeClusterCollection>();
       
       unsigned nbins = (unsigned)rint((tmax_-tmin_)/tbin_);
       tmax_ = tmin_+ nbins*tbin_;
       
       //if (diag_>0) diagTool_ = art::make_tool<ModuleHistToolBase>(config.get_PSet().get<fhicl::ParameterSet>("diagPlugin"));       
       if (diag_>0) diagTool_ = art::make_tool<ModuleHistToolBase>(config().diagPlugin," ");
       
       if (diag_>1)
       {
          std::cout<<"[TimeAndPhiClusterFinder] Request MVA diagnosis, setting MVA MVAScoreCut to 0.0"<<std::endl;
          MVAScoreCut_ = 0.0;
       }
  }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::beginJob()
  {
      tcMVA_.initMVA();
      
      if (diag_ > 0)
      {
          art::ServiceHandle<art::TFileService> tfs;
          diagTool_->bookHistograms(tfs);
      }
   }


  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::produce(art::Event & event )
  {

      if (diag_ > 0 && (iev_%printfreq_)==0) std::cout<<"TimeAndPhiClusterFinder: event="<<iev_<<std::endl;

      const auto& chH = event.getValidHandle(chToken_);
      chcol_ = chH.product();

      art::Handle<CaloClusterCollection> ccH{}; 
      if (usecc_) {ccH = event.getHandle<CaloClusterCollection>(ccToken_);}

      if (testflag_)
      {
          const auto& shfH = event.getValidHandle(shfToken_);
          shfcol_ = shfH.product();
          if (shfcol_->size() != chcol_->size())
	     throw cet::exception("RECO")<<"TimeAndPhiClusterFinder: inconsistent flag collection length " << std::endl;
      }

      if (diag_>0) {data_.reset(); data_.event_=&event; data_.chcol_ = chcol_;}
      
      std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
      findClusters(ccH,*tccol);

      if (diag_>0) diagTool_->fillHistograms(&data_);

      event.put(std::move(tccol));
      ++iev_;
  }


  

  //--------------------------------------------------------------------------------------------------------------
  void TimeAndPhiClusterFinder::findClusters(const art::Handle<CaloClusterCollection>& ccH, TimeClusterCollection& tccol)
  {            
      
      // Cache good ComboHit and their recalculated time. Map hits per time slices 
      // (the 0.02 corresponds to 50 bins, an internal parameter that isn't part of the interface)
      
      unsigned hitIndexSize = unsigned((tmax_-tmin_)*0.02f);
      std::vector<std::vector<unsigned>> hitIndex(hitIndexSize,std::vector<unsigned>());
      for (auto& vec : hitIndex) vec.reserve(128);
           
      std::vector<float> comboHitTime(chcol_->size(),0.0);
      for (unsigned ich=0; ich<chcol_->size();++ich)
      {
          if (testflag_ && !goodHit((*shfcol_)[ich])) continue;
          float time        = ttcalc_.comboHitTime((*chcol_)[ich],pitch_);                
          comboHitTime[ich] = time;
          
          unsigned idx = unsigned((time-tmin_)*0.02f);
          if (idx <=hitIndexSize) hitIndex[idx].push_back(ich);          
      }


      // Start by finding time peaks
      std::vector<TimeCandidate> timeCandidates;
      timeCandidates.reserve(32);
      findTimePeaks(ccH, timeCandidates, comboHitTime);     

      // Then find the phi clusters from those time clusters and init the time cluster
      for (auto& tc : timeCandidates) findPhiPeaks(tc,tccol);
      std::remove_if(tccol.begin(),tccol.end(),[this](const auto& tc){return tc._nsh<minNSHits_;});            
      
      for (auto& tc : tccol) calculateMean(tc,comboHitTime);
      if (diag_==1 || diag_==3) fillDiagPreClean(timeCandidates, tccol);
    

      // Finally refine the clusters and remove those with too few hits
      auto it = tccol.begin();
      while (it !=  tccol.end()) 
      {            
          filterClusterSequential(*it,comboHitTime);          
          if (recover_)   recoverHits(*it,comboHitTime, hitIndex);

          if (it->nStrawHits() < minNSHits_) it = tccol.erase(it); 
          else          ++it;          
      }

      if (diag_==1 || diag_==3) fillDiagPostClean(tccol);      
  }

  


  //--------------------------------------------------------------------------------------------------------------
  // Find peaks in time distribution
  //--------------------------------------------------------------------------------------------------------------
  
  //--------------------------------------------------------------------------------------------------------------
  // Find peaks in time distribution
  void TimeAndPhiClusterFinder::findTimePeaks(const art::Handle<CaloClusterCollection>& ccH, TimeCandidateCollection& timeCandidates,
                                              const std::vector<float>& comboHitTime)
  {            
      std::vector<unsigned> timePeaks;
      const unsigned nbins = unsigned((tmax_-tmin_)/tbin_);
      std::vector<unsigned> timeHist(nbins,0);
     
      TimeCandidateCollection timeCaloCandidates;
      if (usecc_) findCaloSeeds(ccH, timeCaloCandidates); 
      
      for (unsigned ich=0; ich<comboHitTime.size();++ich)
      {
          if (comboHitTime[ich]<1.0) continue;
          unsigned ibin = unsigned((comboHitTime[ich]-tmin_)/tbin_);
          if (ibin <= nbins) timeHist[ibin] += (*chcol_)[ich].nStrawHits();
      }
      
      //extract the time peak position
      if (timeStrategy_==1) timePeakScan(timeHist,timePeaks);
      else                  timePeakBlind(timeHist,timePeaks,nbins);

      if (timePeaks.empty()) return;
      
      // Refine the time estimate for each peak and create time cluster candidate
      for (const auto& ipeak : timePeaks)
      {          
          float nsh(0),t0(0);       
          for (unsigned ibin = ipeak-windowTime_;ibin <= std::min(nbins,ipeak+windowTime_); ++ibin)
          {
             float time = tmin_ + tbin_*(ibin + 0.5);
             t0  += timeHist[ibin]*time;
             nsh += timeHist[ibin];          
          }
          t0/=nsh;
                    
          TimeCandidate tc;
	  tc.t0_ = t0; 
          timeCandidates.push_back(tc);            
      } 
             
      // Assign hit to cluster, and remove clusters with small number of hits (might be associated to calo cluster and non-calo cluster)
      for (unsigned ich=0;ich<comboHitTime.size();++ich)
      {
          float time0 = comboHitTime[ich];
          if (time0<1.0) continue;
                    
          // look at the calorimeter clusters if needed
          if (usecc_ && !timeCaloCandidates.empty())
          {
              auto besttcc = std::min_element(timeCaloCandidates.begin(),timeCaloCandidates.end(),
                                              [time0](const auto& a, const auto& b)
                                              {return std::abs(a.t0_-time0) < std::abs(b.t0_-time0);});                        
              if (std::abs(besttcc->t0_-time0) < maxTimeDTCal_) {besttcc->strawIdx_.emplace_back(ich); besttcc->nsh_ += chcol_->at(ich).nStrawHits();}         
          }
                            
          auto besttc = std::min_element(timeCandidates.begin(),timeCandidates.end(),
                                         [time0](const auto& a, const auto& b)
                                         {return std::abs(a.t0_-time0) < std::abs(b.t0_-time0);});                        
          if (std::abs(besttc->t0_-time0) < maxTimeDT_) {besttc->strawIdx_.emplace_back(ich); besttc->nsh_ += chcol_->at(ich).nStrawHits();}
      }
                 
      // Final filtering
      if (usecc_) std::move(timeCaloCandidates.begin(),timeCaloCandidates.end(),std::back_inserter(timeCandidates));
      std::remove_if(timeCandidates.begin(),timeCandidates.end(),[this](const auto& a){return a.nsh_<minNSHits_;});
   }
         
   //--------------------------------------------------------------------------------------------------------------
   void TimeAndPhiClusterFinder::timePeakScan(const std::vector<unsigned>& timeHist, std::vector<unsigned>& timePeaks)
   {   
       unsigned k(2u*windowTime_);
       std::deque<unsigned> q;      
       for (unsigned i=0;i<timeHist.size();++i)
       {
           while (!q.empty() && timeHist[i] > timeHist[q.back()]) q.pop_back();
           while (!q.empty() && q.front()+k < i)                  q.pop_front();
           q.push_back(i);                            

           if (i<k && timeHist[i]>timeHist[i+1] && timeHist[i]>minTimeYbin_) timePeaks.push_back(i);
           if (q.front() > windowTime_ && q.front()+windowTime_==i && timeHist[q.front()]>minTimeYbin_)
           {
              //check if the next bin has the same value and break the tie by looking at neighbors 
              unsigned idx = q.front();
              if (idx+2 < timeHist.size() && timeHist[idx]==timeHist[idx+1] && timeHist[idx-1]>timeHist[idx+2]) ++idx;
              timePeaks.push_back(idx);
           }
       }     
       std::sort(timePeaks.begin(),timePeaks.end(),[&timeHist](unsigned i, unsigned j){return timeHist[i]>timeHist[j];});   
   }

   //--------------------------------------------------------------------------------------------------------------
   void TimeAndPhiClusterFinder::timePeakBlind(const std::vector<unsigned>& timeHist, std::vector<unsigned>& timePeaks, const unsigned nbins )
   {   
       std::vector<bool> used(nbins,false); 
       std::vector<unsigned> indices(nbins);
       std::iota (indices.begin(), indices.end(), 0);
       std::sort(indices.begin(),indices.end(),
                [&timeHist](unsigned i, unsigned j)
                {if (timeHist[i]!=timeHist[j]) return timeHist[i]>timeHist[j]; return i<j;});

       for (const auto& idx : indices)
       {
           if (used[idx]) continue;
           timePeaks.push_back(idx);
           used[idx] = true;
           if (idx >0 ) used[idx-1] = true;
           if (idx+1<nbins) used[idx+1] = true; 
           if (timeHist[idx]<=minTimeYbin_) break;  
       }
   }
   
   //-------------------------------------------------------------------------------------------------------------------
   void TimeAndPhiClusterFinder::findCaloSeeds(const art::Handle<CaloClusterCollection>& ccH, TimeCandidateCollection& timeCandidates)
   {
       const CaloClusterCollection& cccol = *ccH.product();
       for (size_t icalo=0; icalo < cccol.size(); ++icalo)
       {
           const auto& calo = cccol.at(icalo);
           if (calo.energyDep() > ccmine_)
           {
	      TimeCandidate tc;
	      tc.t0_          = ttcalc_.caloClusterTime(calo,pitch_);
              tc.caloCluster_ = art::Ptr<CaloCluster>(ccH,icalo);
	      timeCandidates.push_back(tc);

              if (diag_==1 || diag_==3) data_.calTime_[data_.Ncal_++] = tc.t0_ ;
          }
       }
   }
  
  
  
  
  
  
  
  //--------------------------------------------------------------------------------------------------------------
  // FIND PEAKS IN PHI DISTRIBUTION
  //--------------------------------------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------------------------------------
  // Find peaks in phi coordinates for a given time cluster. This is a full TimeCluster candidate
  void TimeAndPhiClusterFinder::findPhiPeaks(TimeCandidate& tc, TimeClusterCollection& tccol) 
  {
       const static unsigned len = unsigned(2*M_PI/pbin_)+1;
       std::vector<unsigned> phiHist(len,0);

       for (auto ish : tc.strawIdx_) 
       {
            const ComboHit& ch = chcol_->at(ish);
            float phi          = ch.pos().phi()+M_PI;
            unsigned ibin      = unsigned(phi/pbin_);
            phiHist[ibin]     += ch.nStrawHits();
       }
       
       std::vector<unsigned> istart,iend;
       associateScan(phiHist,istart,iend);

       for (unsigned ic=0;ic<istart.size();++ic)
       {
           TimeCluster tclu;
           tclu._t0 = TrkT0(tc.t0_,tbin_/0.5);
           tclu._caloCluster = tc.caloCluster_;
           tclu._strawHitIdxs.reserve(16);
           
           for (const auto sIdx :tc.strawIdx_ )
           {
               const ComboHit& ch = chcol_->at(sIdx);
               float phi          = ch.pos().phi()+M_PI;
               unsigned ibin      = unsigned(phi/pbin_);
               if (ibin < istart[ic]) ibin +=len;        
                                     
               if (ibin >= istart[ic] && ibin<iend[ic])
               {
                   tclu._nsh +=ch.nStrawHits(); 
                   tclu._strawHitIdxs.emplace_back(sIdx);
               }
           } 
           if (tclu._nsh > minNSHits_) tccol.push_back(tclu);
       }        
   }
    
   //---------------------------------------------------------------------------------------------------------------------  
   // Associate hits to a given phi cluster by scanning the phi vector
   void TimeAndPhiClusterFinder::associateScan(const std::vector<unsigned>& phiHist, std::vector<unsigned>& istart, std::vector<unsigned>& iend)
   {
      const unsigned len(phiHist.size());
      
      //find the first end of sequence (two consecutive 0) and create the buffer vector
      unsigned offset(0);
      while (offset+2<len && (phiHist[offset]+phiHist[offset+1])>0) ++offset;

      std::vector<unsigned> phiHistWork(len+offset+2);
      std::copy(phiHist.begin(),phiHist.end(), phiHistWork.begin());
      std::copy(phiHist.begin(),phiHist.begin()+offset+2,phiHistWork.end()-offset-2);

      
      // search for isolated sequences surreounded by separator sequence (0 0)
      for (unsigned is=std::max(2u,offset); is+2<phiHistWork.size();++is)
      {
          //move until you see the pass the separator sequence and find the end
          if (phiHistWork[is]==0 || (phiHistWork[is-1]+phiHistWork[is-2])>0) continue;

          unsigned ie(is+1),maxVal(phiHist[is%len]),content(0);
          while (ie%len != is && (phiHist[ie%len]+phiHist[(ie+1)%len])>0)
          {
             content += phiHistWork[ie%len];
             maxVal = std::max(maxVal,phiHistWork[ie%len]);
             ++ie;
          }
                                        
          if (ie-is<3) continue;
                             
          // try to break two large clusters close to each other separated by small distance
          if (ie-is > 6 && content > 2*minHitSplit_)
          {
              unsigned minSum(99),imin(99),leftSum(0),rightSum(0);
              for (unsigned i=is+2;i<ie-2;++i)
              {
                  unsigned sum = phiHist[(i-1)%len]+phiHist[i%len]+phiHist[(i+1)%len];
                  if (sum < minSum) {minSum=sum; imin=i;}          
              }
              for (unsigned i=is;i<imin;++i) leftSum += phiHist[i%len];
              for (unsigned i=imin;i<ie;++i) rightSum += phiHist[i%len];

              if ((leftSum>minHitSplit_ && rightSum > minHitSplit_) && 2*phiHist[imin%len] < maxVal)
              {                 
                  istart.push_back(is);    iend.push_back(imin);
                  istart.push_back(imin);  iend.push_back(ie);
                  continue;
              }
           }
 
           // save slightly larger sequences to recollect stragglers
           istart.push_back(is>3? is-3 : 0);
           iend.push_back(ie+3);
      }
      
      
      // if there is only one huge sequence, try looser separator
      if (istart.empty())
      {
          unsigned is0(0);
          for (unsigned i=2;i<phiHist.size()-2;++i)
          {             
              if (phiHist[(i-1)%len]+phiHist[(i+1)%len]+phiHist[i%len]>1) continue;
              if (phiHist[i]==1) ++i;
              if (i-is0>1) { istart.push_back(is0);iend.push_back(i);}
              is0=i;
          }
          if (phiHist.size()-is0>1) {istart.push_back(is0);iend.push_back(phiHist.size());}
      }
     
      //though luck, put all hits inside a single sequence and pray it works! Should be rare though
      if (istart.empty()) {istart.push_back(0); iend.push_back(len);}
                 
      return;   
   }





  //--------------------------------------------------------------------------------------------------------------
  // Filtering and recovering hits
  //--------------------------------------------------------------------------------------------------------------
 
  //--------------------------------------------------------------------------------------------------------------
  // sequential removal of bad hits
  void TimeAndPhiClusterFinder::filterClusterSequential(TimeCluster& tc, const std::vector<float>& comboHitTime)
  {     
     std::vector<float> mva(4,0.0), mvaWorst;
     while (tc._strawHitIdxs.size() >= minNSHits_)
     {
         float cluPhiMean = polyAtan2( tc._pos.y(), tc._pos.x());
         float cluRadMean = sqrtf(tc._pos.y()*tc._pos.y()+tc._pos.x()*tc._pos.x());
         float cluTime    = tc._t0._t0; 

         float scoreWorse(0);
         auto it(tc._strawHitIdxs.begin()), worseHit(tc._strawHitIdxs.end());
         while (it != tc._strawHitIdxs.end())
         {
             const ComboHit& ch = chcol_->at(*it);
	     float phi          = ch.pos().phi();
             float rad          = sqrtf(ch.pos().x()*ch.pos().x()+ch.pos().y()*ch.pos().y());          

             mva[0] = M_PI - abs(abs(phi-cluPhiMean) - M_PI); 
             mva[1] = abs(comboHitTime[*it] - cluTime);
             mva[2] = abs(rad-cluRadMean);
             mva[3] = rad;

             float score = tcMVA_.evalMVA(mva);
             if (score > scoreWorse) {worseHit=it; scoreWorse=score; if (diag_>1) mvaWorst=mva;} 
             ++it;
         }

         if (scoreWorse < MVAScoreCut_) break;
         if (diag_>1) fillDiagMVA(*worseHit, mvaWorst);
         
         tc._strawHitIdxs.erase(worseHit);                              
         calculateMean(tc,comboHitTime); //nearly as fast as recalculating time/pos by removing one hit but much cleaner!
     }
     if (diag_>1) ++data_.icluMVA_;
  }
  
  //-----------------------------------------------------------------------------------------------
  // Hit recovery selection
  void TimeAndPhiClusterFinder::recoverHits(TimeCluster& tc, const std::vector<float>& comboHitTime, 
                                            const std::vector<std::vector<unsigned>>& hitIndex)
  {           
      float cluPhiMean = polyAtan2(tc._pos.y(), tc._pos.x());
      float cluRadMean = sqrtf(tc._pos.y()*tc._pos.y()+tc._pos.x()*tc._pos.x());
      float cluTime    = tc._t0._t0;
      
      bool recover(false);
      std::vector<float> mva(4,0.0);
      int idx = int((tc._t0._t0-tmin_)*0.02f);     
      for (int i=std::max(0,idx-1);i<idx+1;++i)
      {
          for (auto ich : hitIndex[i])
          {
              if (std::abs(comboHitTime[ich]-cluTime) > 40.0f) continue;
              if (std::find(tc._strawHitIdxs.begin(),tc._strawHitIdxs.end(),ich) != tc._strawHitIdxs.end()) continue;

              const ComboHit& ch = (*chcol_)[ich];
              float phi          = ch.pos().phi();
              float dphi         = M_PI - abs(abs(phi-cluPhiMean) - M_PI); 
              if (dphi > dphiMaxReco_) continue;

              float dtime = abs(comboHitTime[ich] - cluTime);
              float rad   = sqrtf(ch.pos().x()*ch.pos().x()+ch.pos().y()*ch.pos().y());          

              mva[0] = dphi; 
              mva[1] = dtime;
              mva[2] = abs(rad-cluRadMean);
              mva[3] = rad;

              float score = tcMVA_.evalMVA(mva);            
              if (score < MVAScoreCut_) {tc._strawHitIdxs.emplace_back(ich);recover=true;}
          }
      }
      
      if (recover) calculateMean(tc,comboHitTime);
  }
 
  
  
 

  //--------------------------------------------------------------------------------------------------------------
  // calculate the mean position / time of cluster (as good as median in this case)
  void TimeAndPhiClusterFinder::calculateMean(TimeCluster& tc, const std::vector<float>& comboHitTime) 
  {        
      if (tc._strawHitIdxs.empty()) {tc._pos = XYZVec(0,0,0);tc._t0._t0=0; tc._t0._t0err =0; tc._nsh = 0; return;};
      
      tc._nsh = 0;            
      float tacc(0),tacc2(0),xacc(0),yacc(0),zacc(0),weight(0);
      for (auto ish :tc._strawHitIdxs) 
      {
          const ComboHit& ch = chcol_->at(ish);
          tc._nsh += ch.nStrawHits();

          float htime = comboHitTime[ish];
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
      tc._pos       = XYZVec(xacc, yacc, zacc);          
  }


  


  //--------------------------------------------------------------------------------------------------------------
  // Diagnosis
  //--------------------------------------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------------------------------------------------  
  // Diagnosis routines
  void TimeAndPhiClusterFinder::fillDiagPreClean(const TimeCandidateCollection& timeCandidates,const TimeClusterCollection& tccol )
  {
     data_.iev_ = iev_;
     
     data_.Nch_=chcol_->size();
     for (unsigned ich=0; ich<chcol_->size();++ich)
     {
         data_.chTime_[ich]= 0;
         data_.chSel_[ich] = 0;
         if (testflag_ && !goodHit(shfcol_->at(ich))) continue;
         float time         = ttcalc_.comboHitTime(chcol_->at(ich),pitch_);                
         data_.chTime_[ich]= time;
         data_.chSel_[ich] = 1;
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
           data_.clu2Time_[ih2] = tc._t0._t0;
           data_.clu2Phi_[ih2]  = polyAtan2(tc._pos.y(),tc._pos.x());
           data_.clu2Rad_[ih2]  = sqrtf(tc._pos.y()*tc._pos.y()+tc._pos.x()*tc._pos.x());
           ++ih2; 
        }                                                      
        ++iclu2;            
     }
     data_.nhit2_=ih2;
  }

  void TimeAndPhiClusterFinder::fillDiagPostClean(const TimeClusterCollection& tccol)
  {
      int iclu3(0), ih3(0);         
      for (const auto& tc : tccol) 
      {
          for (auto ish : tc._strawHitIdxs)
          {
             data_.nclu3_[ih3]  = iclu3; 
             data_.hitIdx3_[ih3]= ish; 
             ++ih3; 
          }                                                      
          ++iclu3;    
      }
      data_.nhit3_=ih3;
   }
  
   void TimeAndPhiClusterFinder::fillDiagMVA(unsigned iworseHit, const std::vector<float>& mva)
   {      
       data_.hitIdxMVA_[data_.nhitMVA_] = iworseHit;
       data_.MVAvar1_[data_.nhitMVA_]   = mva[0];
       data_.MVAvar2_[data_.nhitMVA_]   = mva[1];
       data_.MVAvar3_[data_.nhitMVA_]   = mva[2];
       data_.MVAvar4_[data_.nhitMVA_]   = mva[3];
       data_.MVAvar5_[data_.nhitMVA_]   = mva[4];
       data_.ncluMVA_[data_.nhitMVA_]   = data_.icluMVA_;
       ++data_.nhitMVA_;
   }
  
}


DEFINE_ART_MODULE(mu2e::TimeAndPhiClusterFinder);

