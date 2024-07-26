#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/TrkPatRec/inc/RobustMultiHelixFinder_types.hh"

#include <numeric>
#include <stdint.h>


namespace {

  struct CandHelix
  {
    using artCptr = art::Ptr<mu2e::CaloCluster>;

    CandHelix() : nStrawHits_(0),x_(0),y_(0),r_(0),fita_zp_(0),fitb_zp_(0),fita_zt_(0),fitb_zt_(0),caloPtr_{} {};

    CandHelix(int nsh, float x, float y, float r, const std::vector<size_t>& hitIdxs, artCptr& caloPtr) :
       nStrawHits_(nsh), x_(x), y_(y),r_(r), fita_zp_(0),fitb_zp_(0),fita_zt_(0),fitb_zt_(0),
       hitIndices_(hitIdxs),caloPtr_(caloPtr)
    {hitIndices_.reserve(32);};

    void remove(unsigned it, unsigned weight) {hitIndices_.erase(hitIndices_.begin()+it); nStrawHits_ -= weight;}

    unsigned              nStrawHits_;
    float                 x_,y_,r_,fita_zp_,fitb_zp_,fita_zt_,fitb_zt_;
    std::vector<size_t>   hitIndices_;  //combo hit indices for the helix candidate
    artCptr caloPtr_;
  };

  struct LSFitter
  {
    LSFitter() : sn_(0),sx_(0),sx2_(0),sy_(0),sxy_(0) {};

    float fa()  {return fabs(sn_*sx2_-sx_*sx_)>1e-6 ? (sn_*sxy_-sx_*sy_) /(sn_*sx2_-sx_*sx_) : 0.0;}
    float fb()  {return fabs(sn_*sx2_-sx_*sx_)>1e-6 ? (sy_*sx2_-sx_*sxy_)/(sn_*sx2_-sx_*sx_) : 0.0;}
    void  add(float x, float y, float w=1.0) {sx_+=x*w; sx2_+=x*x*w; sy_+=y*w; sxy_+=x*y*w; sn_+=w;}

    float sn_,sx_,sx2_,sy_,sxy_;
  };
}



namespace mu2e {

  class RobustMultiHelixFinder : public art::EDProducer
  {
    public:
      using mapHelix        = std::map<Helicity,std::unique_ptr<HelixSeedCollection>>;
      using strawHitIndices = std::vector<StrawHitIndex>;
      using Config_types    = RobustMultiHelixFinderTypes::Config ;
      using Data_types      = RobustMultiHelixFinderTypes::Data_t;

      enum circleFitter {NoChoice, HyperFit, ChiSquared};

      struct Config
      {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<art::InputTag> comboHitCollection    {Name("ComboHitCollection"),     Comment("ComboHit collection {Name")    };
        fhicl::Atom<art::InputTag> timeClusterCollection {Name("TimeClusterCollection"),  Comment("TimeCluster collection {Name") };
        fhicl::Sequence<int>       helicities            {Name("Helicities"),             Comment("Helicity values") };
        fhicl::Atom<float>         clusteringPhiBin      {Name("ClusteringPhiBin"),       Comment("Compton clustering phi bin ") };
        fhicl::Atom<float>         clusteringMinBin      {Name("ClusteringMinBin"),       Comment("Minimum bin content for a Compton cluster") };
        fhicl::Atom<float>         minDRCircle           {Name("MinDRCircle"),            Comment("Maximum radial distance between circle and hit") };
        fhicl::Atom<float>         minRadCircle          {Name("MinRadCircle"),           Comment("Minimum circle radius") };
        fhicl::Atom<float>         maxRadCircle          {Name("MaxRadCircle"),           Comment("Maximum circle radius") };
        fhicl::Atom<float>         minDXY2Circle         {Name("MinDXY2Circle"),          Comment("Minimum XY distance between hits for circle fit") };
        fhicl::Atom<bool>          targetCon             {Name("TargetCon"),              Comment("Require track to be produced in target") };
        fhicl::Atom<float>         targetRadius          {Name("TargetRadius"),           Comment("Target radius") };
        fhicl::Atom<float>         ccMinEnergy           {Name("CaloClusterMinE"),        Comment("Minimum calo cluster energy") };
        fhicl::Atom<int>           ccWeight              {Name("CaloClusterWeight"),      Comment("Calo cluster weight ") };
        fhicl::Atom<float>         minDPDZSlope          {Name("MinDPDZSlope"),           Comment("Minimum dz/dphi slope for track") };
        fhicl::Atom<float>         maxDPDZSlope          {Name("MaxDPDZSlope"),           Comment("Maximum dz/dphi slope for track") };
        fhicl::Atom<float>         DPDZStep              {Name("DPDZStep"),               Comment("Step in dz/dphi scan for dz/dp fit init ") };
        fhicl::Atom<float>         maxDPhiHelInit        {Name("MaxDPhiHelInit"),         Comment("Maximum phi difference between hit and dz/dphi line init fit") };
        fhicl::Atom<float>         maxDPhiHelFit         {Name("MaxDPhiHelFit"),          Comment("Maximum phi difference between hit and dz/dphi line") };
        fhicl::Atom<float>         maxDTHelFit           {Name("MaxDTHelFit"),            Comment("Maximum time difference between hit time and z-time fit") };
        fhicl::Atom<float>         maxChi2Hit            {Name("MaxChi2Hit"),             Comment("Maximum chi2 for a hit to be associated to the helix") };
        fhicl::Atom<float>         minDZTrk              {Name("MinDZTrk"),               Comment("Minimum z span of a trk") };
        fhicl::Atom<std::string>   fitCircleStr          {Name("FitCircleStrategy"),      Comment("Fit Circle algorhithm HyperFit or ChiSquared") };
        fhicl::Atom<unsigned>      minStrawHits          {Name("MinStrawHits"),           Comment("Minimum number of Straw hits for a helix candidate") };
        fhicl::Atom<unsigned>      nMaxTrkIter           {Name("NMaxTrkIter"),            Comment("Number of track finding iterations ") };
        fhicl::Atom<int>           diagLevel             {Name("DiagLevel"),              Comment("Diag level"), 0 };
        fhicl::Table<Config_types> diagPlugin            {Name("DiagPlugin"),             Comment("Diag Plugin config")};
        fhicl::Atom<float>         maxEDepAvg            {Name("maxEDepAvg"),             Comment("Maximum EDep average")};
      };
      explicit RobustMultiHelixFinder(const art::EDProducer::Table<Config>& config);
      virtual void produce(art::Event& event);
      virtual void beginJob();


    private:
      const art::ProductToken<ComboHitCollection>    chToken_;
      const art::ProductToken<TimeClusterCollection> tcToken_;
      float                                          clusteringPhiBin_;
      float                                          clusteringMinBin_;
      float                                          minDR2Circle_;
      float                                          minRadCircle_;
      float                                          maxRadCircle_;
      float                                          minDXY2Circle_;
      bool                                           targetCon_;
      float                                          targetRadius_;
      float                                          ccMinEnergy_;
      int                                            ccWeight_;
      float                                          minDPDZSlope_;
      float                                          maxDPDZSlope_;
      float                                          DPDZStep_;
      float                                          MaxDPhiHelInit_;
      float                                          maxDPhiHelFit_;
      float                                          maxDTHelFit_;
      float                                          maxChi2Hit_;
      float                                          minDZTrk_;
      circleFitter                                   fitCircleStrategy_;
      unsigned                                       minStrawHits_;
      unsigned                                       nMaxTrkIter_;
      const Calorimeter*                             cal_;
      int                                            iev_;
      int                                            diag_;
      std::unique_ptr<ModuleHistToolBase>            diagTool_;
      float                                          maxEDepAvg_;
      Data_types                                     data_;
      std::vector<Helicity>                          hels_;


      void              findAllHelices        (art::Event& event, mapHelix& helcols, const art::ValidHandle<TimeClusterCollection>& tcH,
                                               const art::ValidHandle<ComboHitCollection>& chH);
      void              findHelicesInTC       (mapHelix& helcols, const art::Ptr<TimeCluster>& tcArtPtr, const TimeCluster& tc, const ComboHitCollection& chcol);
      void              findHelicesInHits     (const ComboHitCollection& chcol, strawHitIndices& hitsToProcess, const art::Ptr<CaloCluster>& caloPtr,
                                               CandHelix& poshelix, CandHelix& neghelix);
      CandHelix         findCircleCandidate   (const ComboHitCollection& chcol, const strawHitIndices& hits, const art::Ptr<CaloCluster>& caloPtr);
      void              findHelixCandidate    (CandHelix& circle,       const ComboHitCollection& chcol, const strawHitIndices& hitsToProcess, Helicity helicity);
      std::vector<bool> flagCompton           (const ComboHitCollection& chcol, const std::vector<StrawHitIndex>& hitsToProcess);
      unsigned          init_dzdp             (CandHelix& circle,       const ComboHitCollection& chcol, Helicity helicity);
      void              fit_dzdp              (CandHelix& circle,       const ComboHitCollection& chcol, Helicity helicity);
      void              fit_dzdt              (CandHelix& circle,       const ComboHitCollection& chcol);
      int               filterZPhi            (CandHelix& circle,       const ComboHitCollection& chcol, float maxDphi);
      int               filterZT              (CandHelix& circle,       const ComboHitCollection& chcol, float maxDt);
      void              fitCircleAlg          (const CandHelix& circle, const ComboHitCollection& chcol, float& centerX, float& centerY, float& radius);
      void              fitCircleChi2         (const CandHelix& circle, const ComboHitCollection& chcol, float& centerX, float& centerY, float& radius);
      float             weight                (const ComboHit& ch,      const CandHelix& circle);
      float             chi2XYHelix           (const ComboHit& ch,      const CandHelix& circle);
      float             chi2XYCircle          (const ComboHit& ch,      float centerX, float centerY, float radius);
      float             chi2XYHelix           (const ComboHit& ch,      float centerX, float centerY, float radius, float dzdp ,float fzb);
      float             timeUncertainty       (const ComboHitCollection& chcol, const CandHelix& helix);
      void              filterDuplicateHelices(HelixSeedCollection& helices);
      void              printHelix            (const HelixSeed& helix);
      void              fillDiag              (mapHelix& helcols, const ComboHitCollection& chcol);
  };


  RobustMultiHelixFinder::RobustMultiHelixFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    chToken_       {consumes<ComboHitCollection>   (config().comboHitCollection())},
    tcToken_       {consumes<TimeClusterCollection>(config().timeClusterCollection())},
    clusteringPhiBin_ (config().clusteringPhiBin()),
    clusteringMinBin_ (config().clusteringMinBin()),
    minDR2Circle_     (config().minDRCircle()*config().minDRCircle()),
    minRadCircle_     (config().minRadCircle()),
    maxRadCircle_     (config().maxRadCircle()),
    minDXY2Circle_    (config().minDXY2Circle()),
    targetCon_        (config().targetCon()),
    targetRadius_     (config().targetRadius()),
    ccMinEnergy_      (config().ccMinEnergy()),
    ccWeight_         (config().ccWeight()),
    minDPDZSlope_     (config().minDPDZSlope()),
    maxDPDZSlope_     (config().maxDPDZSlope()),
    DPDZStep_         (config().DPDZStep()),
    MaxDPhiHelInit_   (config().maxDPhiHelInit()),
    maxDPhiHelFit_    (config().maxDPhiHelFit()),
    maxDTHelFit_      (config().maxDTHelFit()),
    maxChi2Hit_       (config().maxChi2Hit()),
    minDZTrk_         (config().minDZTrk()),
    fitCircleStrategy_(circleFitter::NoChoice),
    minStrawHits_     (std::max(config().minStrawHits(),3u)),
    nMaxTrkIter_      (config().nMaxTrkIter()),
    iev_              (0),
    diag_             (config().diagLevel()),
    diagTool_(),
    maxEDepAvg_       (config().maxEDepAvg()),
    data_()
  {
    std::vector<int> helvals = config().helicities();
    for (const auto& hv : helvals) {
      Helicity hel(hv);
      hels_.emplace_back(hel);
      produces<HelixSeedCollection>(Helicity::name(hel));
    }
    if (diag_) diagTool_ = art::make_tool<ModuleHistToolBase>(config().diagPlugin," ");

    if      (config().fitCircleStr()=="HyperFit")   fitCircleStrategy_ = circleFitter::HyperFit;
    else if (config().fitCircleStr()=="ChiSquared") fitCircleStrategy_ = circleFitter::ChiSquared;
    else    throw cet::exception("CATEGORY")<< "RobustMultiHelixFinder: unrecognixed FitCirclestrategy specified";
  }


  //--------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::beginJob(){
    if (diag_){
       art::ServiceHandle<art::TFileService> tfs;
       diagTool_->bookHistograms(tfs);
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  void RobustMultiHelixFinder::produce(art::Event& event )
  {
    if (diag_>0) std::cout<<"Event "<<event.id().event()<<std::endl;
    iev_ = event.id().event();

    const auto& tcH = event.getValidHandle(tcToken_);
    const auto& chH = event.getValidHandle(chToken_);

    mapHelix helcols;
    for (const auto& hel : hels_) helcols[hel] = std::unique_ptr<HelixSeedCollection>(new HelixSeedCollection());

    findAllHelices(event, helcols,tcH, chH);

    for (const auto& hel : hels_) event.put(std::move(helcols[hel]),Helicity::name(hel));
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Find all helices in an event (hopefully), some timeCluster implementations have overlapping content, so filter duplicates
  void RobustMultiHelixFinder::findAllHelices(art::Event& event, mapHelix& helcols, const art::ValidHandle<TimeClusterCollection>& tcH,
                                              const art::ValidHandle<ComboHitCollection>& chH)
  {
    cal_ = &(*GeomHandle<Calorimeter>());

    const TimeClusterCollection& tccol(*tcH);
    const ComboHitCollection&    chcol(*chH);

    if (diag_) {data_.reset(); data_.event_=&event; data_.chcol_ = &chcol;}

    for (size_t index=0;index<tccol.size();++index) {
      const auto tcArtPtr = art::Ptr<TimeCluster>(tcH,index);
      const auto& tc = tccol[index];
      findHelicesInTC(helcols,tcArtPtr,tc,chcol);
    }

    filterDuplicateHelices(*helcols[Helicity::poshel]);
    filterDuplicateHelices(*helcols[Helicity::neghel]);

    if (diag_) {
      fillDiag(helcols,chcol);
      diagTool_->fillHistograms(&data_);
      if (diag_>2) for (const auto& hel : *helcols[Helicity::poshel]) printHelix(hel);
      if (diag_>2) for (const auto& hel : *helcols[Helicity::neghel]) printHelix(hel);
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Filter helices which have the same set of hits in common
  void RobustMultiHelixFinder::filterDuplicateHelices(HelixSeedCollection& helices)
  {
    for (auto first = helices.begin(); first != helices.end(); ++first){
      for (auto second = std::next(first); second != helices.end(); ++second){

        size_t result(0);
        auto helix1 = first->hits();
        auto helix2 = second->hits();
        auto iter1 = helix1.begin();
        auto iter2 = helix2.begin();
        while (iter1!=helix1.end() && iter2!=helix2.end()) {
           if (iter1->index(0)<iter2->index(0)) ++iter1;
           else if (iter2->index(0)<iter1->index(0)) ++iter2;
           else {++result; ++iter1; ++iter2;}
        }

        if (result != helix1.size() && result != helix2.size()) continue;

        if (first->hits().size()>second->hits().size()) second->_hhits.clear();
        else first->_hhits.clear();
      }
    }

    auto predEmpty = [](const HelixSeed& hel) {return hel.hits().empty();};
    helices.erase(std::remove_if(helices.begin(),helices.end(),predEmpty),helices.end());
  }

  //---------------------------------------------------------------------------------------------------------------------------
  // Find all helices in a timeCluster for both helicities. Can find up to nMaxTrkIter helices for each helicity
  void RobustMultiHelixFinder::findHelicesInTC(mapHelix& helcols, const art::Ptr<TimeCluster>& tcArtPtr, const TimeCluster& tc,
                                               const ComboHitCollection& chcol)
  {
    std::vector<StrawHitIndex> usedHits;
    for (size_t it=0;it<nMaxTrkIter_;++it){

      size_t strawHitsToGo(0);
      strawHitIndices hitsToProcess;
      hitsToProcess.reserve(64);
      for (const auto& ich : tc.hits()) {
        if (find(usedHits.begin(),usedHits.end(),ich) != usedHits.end()) continue;
        hitsToProcess.emplace_back(ich);
        strawHitsToGo += chcol[ich].nStrawHits();
      }
      if (strawHitsToGo < minStrawHits_) break;

      // search for helices for each helicity (do both as it may reduce fake rate)
      CandHelix poshelix, neghelix;
      findHelicesInHits(chcol, hitsToProcess, tc.caloCluster(), poshelix, neghelix);


      // pick the best helix
      const Helicity bestHelicity = poshelix.nStrawHits_ >= neghelix.nStrawHits_ ? Helicity::poshel: Helicity::neghel;
      CandHelix& bestHelix        = poshelix.nStrawHits_ >= neghelix.nStrawHits_ ? poshelix : neghelix;
      for (const auto& ich : bestHelix.hitIndices_) usedHits.push_back(ich);
      if (bestHelix.nStrawHits_ < minStrawHits_) continue;

      // other criteria to reject fake helices
      float deltaZ = chcol[bestHelix.hitIndices_.back()].pos().z() - chcol[bestHelix.hitIndices_.front()].pos().z();
      if (deltaZ < minDZTrk_)  continue;


      //calculate the remaining quantities to fill the RobustHelix object and add the helixSeed to the list
      float chi2dXY(0), chi2dZPhi(0);
      for (const auto& ich : bestHelix.hitIndices_) {
         chi2dXY   += chi2XYCircle(chcol[ich],bestHelix.x_,bestHelix.y_,bestHelix.r_)*chcol[ich].nStrawHits();
         chi2dZPhi += chi2XYHelix(chcol[ich],bestHelix)*chcol[ich].nStrawHits();
      }
      chi2dXY   /= bestHelix.nStrawHits_;
      chi2dZPhi /= bestHelix.nStrawHits_;
//Dirty hack to save particle propagation direction, will be gone when we have updated the data products
chi2dXY = bestHelix.fita_zt_;

      float Rcent  = sqrt(bestHelix.x_*bestHelix.x_+bestHelix.y_*bestHelix.y_);
      float Fcent  = polyAtan2(bestHelix.y_,bestHelix.x_);
      float fz0    = -bestHelix.fitb_zp_/bestHelix.fita_zp_;
      float t0     = bestHelix.fitb_zt_;
      float t0err  = timeUncertainty(chcol,bestHelix);

      HelixSeed hseed;
      hseed._helix = RobustHelix(Rcent,Fcent,bestHelix.r_,bestHelix.fita_zp_,fz0);
      hseed._hhits.setParent(chcol.parent());
      hseed._helix._helicity  = bestHelicity;
      hseed._helix._chi2dXY   = chi2dXY;
      hseed._helix._chi2dZPhi = chi2dZPhi;
      hseed._t0 = TrkT0(t0,t0err);
      for (const auto& ich : bestHelix.hitIndices_) hseed._hhits.emplace_back(chcol[ich]);
      hseed._status.merge(TrkFitFlag::MPRHelix);
      hseed._status.merge(TrkFitFlag::helixOK);
      hseed._timeCluster = tcArtPtr;
      hseed._eDepAvg = hseed._hhits.eDepAvg();
      if (hseed._eDepAvg > maxEDepAvg_) continue;



      helcols[bestHelicity]->push_back(std::move(hseed));
    }
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Find the best helix in a set of hits
  void RobustMultiHelixFinder::findHelicesInHits(const ComboHitCollection& chcol, strawHitIndices& hitsToProcess,
                                                 const art::Ptr<CaloCluster>& caloPtr, CandHelix& poshelix, CandHelix& neghelix)
  {
    auto pred = [&chcol](const auto& i, const auto& j) {return chcol[i].strawId().uniquePanel()<chcol[j].strawId().uniquePanel();};
    sort(hitsToProcess.begin(),hitsToProcess.end(),pred);

    CandHelix circle = findCircleCandidate(chcol, hitsToProcess, caloPtr);
    if (circle.nStrawHits_ < minStrawHits_) return;

    poshelix = circle;
    findHelixCandidate(poshelix, chcol, hitsToProcess, Helicity::poshel);

    neghelix = circle;
    findHelixCandidate(neghelix, chcol, hitsToProcess, Helicity::neghel);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Loop over triplets and find the circle with the maximum hits using the distance between hit and circle as (fast) metric
  // or chi2. Exclude compton hits in triplet search, but include them in circle assignment. Add calo cluster if present
  CandHelix RobustMultiHelixFinder::findCircleCandidate(const ComboHitCollection& chcol, const strawHitIndices& hits,
                                                        const art::Ptr<CaloCluster>& caloPtr)
  {
    CandHelix bestCircle;
    float sumChisqBest(1e6);
    std::vector<size_t> stubby;
    stubby.reserve(32);

    std::vector<bool> hitComptonFlag = flagCompton(chcol,hits);

    for (size_t i=0; i+2<hits.size();++i){
      if (hitComptonFlag[hits[i]]) continue;
      const auto& chi = chcol[hits[i]];
      float x1        = chi.pos().x();
      float y1        = chi.pos().y();
      float rad2x1y1  = x1*x1+y1*y1;

      for (size_t j=i+1; j+1<hits.size(); ++j){
        if (hitComptonFlag[hits[j]]) continue;
        const auto& chj    = chcol[hits[j]];
        float x2           = chj.pos().x()-x1;
        float y2           = chj.pos().y()-y1;
        float x12x12y12y12 = x2*x2+y2*y2;
        if (x12x12y12y12<minDXY2Circle_) continue;
        float rad2x2y2 = chj.pos().x()*chj.pos().x()+chj.pos().y()*chj.pos().y();

        for (size_t k=j+1; k<hits.size(); ++k){
          if (hitComptonFlag[hits[k]]) continue;
          const auto& chk    = chcol[hits[k]];
          float x3           = chk.pos().x()-x1;
          float y3           = chk.pos().y()-y1;
          float x13x13y13y13 = x3*x3+y3*y3;
          if (x13x13y13y13<minDXY2Circle_) continue;

          float x4           = chk.pos().x()-x2;
          float y4           = chk.pos().y()-y2;
          float x23x23y23y23 = x4*x4+y4*y4;
          if (x23x23y23y23<minDXY2Circle_) continue;
          float rad2x3y3 = chk.pos().x()*chk.pos().x()+chk.pos().y()*chk.pos().y();

          float denominator = 2*(x2*y3 - x3*y2);
          float numeratorX  = x12x12y12y12*y3 - x13x13y13y13*y2;
          float numeratorY  = x13x13y13y13*x2 - x12x12y12y12*x3;

          float centerX = numeratorX/denominator+x1;
          float centerY = numeratorY/denominator+y1;
          float radius  = sqrt((x1-centerX)*(x1-centerX) + (y1-centerY)*(y1-centerY));
          float radCenter2 = centerX*centerX+centerY*centerY;

          if (radius < minRadCircle_ || radius > maxRadCircle_ ) continue;
          if (radCenter2 > rad2x1y1 || radCenter2 > rad2x2y2 || radCenter2 > rad2x3y3) continue;

          // Optional target consistency criteria
          if (targetCon_ && abs(sqrt(radCenter2)-radius)>targetRadius_) continue;

          // Select hits with good deltaR and pick best deltaR if several hits have same unique panel id
          // could alternatively use the chi2
          stubby.clear();
          unsigned id_prev(chcol[hits[0]].strawId().uniquePanel()), ilast(chcol.size()), nsh(0);
          float    val_prev(1e6),sumChisq(0);

          for (size_t i=0;i<hits.size();++i){
             const auto& ch = chcol[hits[i]];
             float dx       = ch.pos().x() - centerX;
             float dy       = ch.pos().y() - centerY;
             float deltaR   = sqrt(dx*dx+dy*dy)-radius;
             float chisq  = deltaR*deltaR;
             //float chisq    = chi2XYCircle(ch,centerX,centerY,radius);

             if (chisq > minDR2Circle_) continue;
             //if (chisq > maxChi2Hit_) continue;

             if (ch.strawId().uniquePanel()!=id_prev && ilast<chcol.size() ) {
                 stubby.push_back(hits[ilast]);
                 nsh      += chcol[hits[ilast]].nStrawHits();
                 val_prev  = chisq;
                 sumChisq += val_prev;
                 id_prev   = ch.strawId().uniquePanel();
                 ilast     = i;
             }
             else if (chisq<val_prev) {val_prev = chisq; ilast=i;}
          }
          if (ilast<chcol.size()) {stubby.push_back(hits[ilast]);nsh += chcol[hits[ilast]].nStrawHits(); sumChisq+=val_prev;}
          sumChisq /= nsh;

          art::Ptr<CaloCluster> thisCaloPtr{};
          if (caloPtr && caloPtr->energyDep()>ccMinEnergy_) {
             float dx      = caloPtr->cog3Vector().x() - centerX;
             float dy      = caloPtr->cog3Vector().y() - centerY;
             float deltaR  = sqrt(dx*dx+dy*dy)-radius;
             float deltaR2 = deltaR*deltaR;
             if (deltaR2 < minDR2Circle_) nsh += ccWeight_;
             thisCaloPtr = caloPtr;
          }

          if (nsh < bestCircle.nStrawHits_ || (nsh==bestCircle.nStrawHits_ && sumChisq > sumChisqBest)) continue;
          bestCircle = CandHelix(nsh,centerX,centerY,radius,stubby,thisCaloPtr);
          sumChisqBest = sumChisq;
        }
      }
    }

    //Recalculate the number of Straw Hits to exclude the calo cluster contribution
    bestCircle.nStrawHits_ = 0;
    for (const auto& ich : bestCircle.hitIndices_) bestCircle.nStrawHits_ += chcol[ich].nStrawHits();
    if (bestCircle.caloPtr_) bestCircle.nStrawHits_ += ccWeight_;

    return bestCircle;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Simple algorithms to find peaks in phi distribution and flag hits in largest peak if above threshiold
  std::vector<bool> RobustMultiHelixFinder::flagCompton(const ComboHitCollection& chcol, const std::vector<StrawHitIndex>& hitsToProcess)
  {
     std::vector<bool> comptonFlags(chcol.size(),false);

     std::vector<float> phiVec;
     for (const auto& ch : chcol) {
        float phi = polyAtan2(ch.pos().y(),ch.pos().x());
        if (phi<0) phi += 6.28318;
        phiVec.push_back(phi);
     }

     int nBins = int(6.28318/clusteringPhiBin_)+1;
     std::vector<int> phiCluster(nBins,0);
     for (const auto ich : hitsToProcess) {
       int ibin = int(phiVec[ich]/clusteringPhiBin_);
       if (ibin<0) ibin=0;
       phiCluster[ibin] += chcol[ich].nStrawHits();
     }

     //compute the rolling sum over three bins and find the maximum.
     //use cyclic buffer trick to go past the end of the veccctor to correctly llok at boundaries
     int sumMax(0),imax(0);
     for (int i=0;i<nBins+3;++i){
       int sum = phiCluster[(i+nBins-1)%nBins] + phiCluster[i%nBins] + phiCluster[(i+1)%nBins];
       if (sum > sumMax) {sumMax=sum;imax=i;}
     }

     if (sumMax<clusteringMinBin_) return comptonFlags;

     float phiMax = (imax+0.5)*clusteringPhiBin_;
     for (const auto& ich : hitsToProcess){
       float dphi = phiVec[ich]-phiMax;
       if (dphi > 6.2831) dphi -= 6.2831;
       if (dphi < -6.2831) dphi += 6.2831;
       if (abs(dphi) < 1.5*clusteringPhiBin_) comptonFlags[ich]=true;
     }

     return comptonFlags;
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Find helix candidate from circle by fitting hits in z-phi plane wiht a robust search for init parameters. Then time-z fit
  // Finally refit circle with algebraic or chi2 fit and filter bad chi2
  void RobustMultiHelixFinder::findHelixCandidate(CandHelix& helix, const ComboHitCollection& chcol,
                                                  const strawHitIndices& hitsToProcess, Helicity helicity)
  {
    //perform dz/dphi fit and filter hits based on the fit result
    init_dzdp(helix,chcol,helicity);
    fit_dzdp(helix,chcol,helicity);
    if (abs(helix.fita_zp_)<1e-3) {helix.hitIndices_.clear(); helix.nStrawHits_=0; return;}
    filterZPhi(helix,chcol,maxDPhiHelFit_);

    //Perform a dz/dt fit and filter hits based on the fit result
    fit_dzdt(helix,chcol);
    filterZT(helix,chcol,maxDTHelFit_);

    //refit helix with algebraic or chisq fit and remove worst chi2 hit if above threshold
    //repeat until no hit is above threshold or we have too few hits
    while (helix.nStrawHits_>minStrawHits_){
      float helixFitX(0),helixFitY(0),helixFitR(0);
      if (fitCircleStrategy_==circleFitter::HyperFit) fitCircleAlg (helix,chcol,helixFitX,helixFitY,helixFitR);
      else                                            fitCircleChi2(helix,chcol,helixFitX,helixFitY,helixFitR);
      if (helixFitR < minRadCircle_ || helixFitR > maxRadCircle_) break;
      helix.x_ = helixFitX;
      helix.y_ = helixFitY;
      helix.r_ = helixFitR;

      unsigned iWorst(-1);
      float chi2Worst(0),nshWorst(0);
      for (size_t i=0;i<helix.hitIndices_.size();++i) {
        unsigned ich = helix.hitIndices_[i];
        float chi2   = chi2XYCircle(chcol[ich],helix.x_,helix.y_,helix.r_);
        if (chi2>chi2Worst) {chi2Worst = chi2; iWorst = i; nshWorst=chcol[ich].nStrawHits();}
      }
      if (chi2Worst < maxChi2Hit_) break;
      helix.remove(iWorst, nshWorst);
    }

    //final dz/dphi and dz/dt refit
    fit_dzdp(helix,chcol,helicity);
    fit_dzdt(helix,chcol);

    helix.nStrawHits_=0;
    for (const auto& ich : helix.hitIndices_) helix.nStrawHits_ += chcol[ich].nStrawHits();
    if (helix.caloPtr_) helix.nStrawHits_ += ccWeight_;

    return;
  }




  //---------------------------------------------------------------------------------------------------------------------------
  // Initialize the dz/dphi linear fit by performing a grid search for the line parameters using calo cluster if present
  // Loop factor is determined dynamically and must increase by 0 or 1 between consecutive z-ordered hits. Add calo cluster if present
  unsigned RobustMultiHelixFinder::init_dzdp(CandHelix& circle, const ComboHitCollection& chcol, Helicity helicity)
  {
    if (circle.hitIndices_.empty()) return 0;

    float caloPhi(-999),caloZ(-9999);
    if (circle.caloPtr_ && circle.caloPtr_->energyDep()>ccMinEnergy_){
      const auto& caloCluster    = *(circle.caloPtr_);
      const auto  caloPosInMu2e  = cal_->geomUtil().diskFFToMu2e(caloCluster.diskID(),caloCluster.cog3Vector());
      const auto  caloPosInTrk   = cal_->geomUtil().mu2eToTracker(caloPosInMu2e);
      caloZ                      = caloPosInTrk.z();
      caloPhi                    = polyAtan2(caloPosInTrk.y()-circle.y_,caloPosInTrk.x()-circle.x_);
    }

    const auto& chits = circle.hitIndices_;
    std::vector<float> p_vec(chcol.size(),-999.9);
    for (const auto& ich : chits) p_vec[ich] = polyAtan2(chcol[ich].pos().y() - circle.y_, chcol[ich].pos().x() - circle.x_);

    unsigned nhitsMax(0);
    float    sumDeltaMax(1e6);
    float    fa_init = helicity == Helicity::poshel ? minDPDZSlope_ : -maxDPDZSlope_;
    float    fa_end  = helicity == Helicity::poshel ? maxDPDZSlope_ : -minDPDZSlope_;
    int      nsteps  = int(std::round((fa_end-fa_init)/DPDZStep_));
    for (int i=0;i<=nsteps;++i){
      float fa = fa_init + i*DPDZStep_;
      for (unsigned j=0;j<chits.size();j+=2){
         float fb = chcol[chits[j]].pos().z()-fa*p_vec[chits[j]];

         unsigned nhits(0);
         float    sumDelta(0);
         int      nloop_prev(round((p_vec[chits[0]]-(chcol[chits[0]].pos().z()-fb)/fa)/6.283185));
         for (const auto& ich : chits){
            int   nloop = round((p_vec[ich]-(chcol[ich].pos().z()-fb)/fa)/6.283185);
            float p2    = p_vec[ich] - nloop*6.283185;
            float delta = abs((chcol[ich].pos().z()-fb)/fa-p2);

            //abort if we skip more than one loop from one hit to the next
            if (nloop>nloop_prev+1) {nhits=0; break;}
            if (delta > MaxDPhiHelInit_) continue;

            nloop_prev = nloop;
            nhits     += chcol[ich].nStrawHits();
            sumDelta  += abs(delta);
         }

         if (circle.caloPtr_){
            int   nloop = round((caloPhi-(caloZ-fb)/fa)/6.283185);
            float delta = abs((caloZ-fb)/fa-caloPhi + nloop*6.283185);
            if (delta < MaxDPhiHelInit_) nhits += ccWeight_;
         }

         if (nhits<nhitsMax || (nhits==nhitsMax && sumDelta < sumDeltaMax)) continue;
         circle.fita_zp_ = fa;
         circle.fitb_zp_ = fb;
         nhitsMax        = nhits;
         sumDeltaMax     = sumDelta;
       }
    }
    return nhitsMax;
  }

  //---------------------------------------------------------------------------------------------------------------------------
  // Perform dz/dphi linear fit, add calo cluster if present
  void RobustMultiHelixFinder::fit_dzdp(CandHelix& circle, const ComboHitCollection& chcol, Helicity helicity)
  {
    const auto& chits = circle.hitIndices_;

    LSFitter zphiFitter;
    for (const auto& ich : chits){
       float phiC    =  polyAtan2(chcol[ich].pos().y() - circle.y_, chcol[ich].pos().x() - circle.x_);
       int   nloop   = round((phiC-(chcol[ich].pos().z()-circle.fitb_zp_)/circle.fita_zp_)/6.283185);
       float phiLoop = phiC - nloop*6.283185;
       float delta   = abs((chcol[ich].pos().z()-circle.fitb_zp_)/circle.fita_zp_- phiLoop);
       if (delta < MaxDPhiHelInit_) zphiFitter.add(phiLoop,chcol[ich].pos().z(),chcol[ich].nStrawHits());
    }

    if (circle.caloPtr_ && circle.caloPtr_->energyDep()>ccMinEnergy_){
      const auto& caloCluster    = *(circle.caloPtr_);
      const auto  caloPosInMu2e  = cal_->geomUtil().diskFFToMu2e(caloCluster.diskID(),caloCluster.cog3Vector());
      const auto  caloPosInTrk   = cal_->geomUtil().mu2eToTracker(caloPosInMu2e);
      float caloPhi              = polyAtan2(caloPosInTrk.y()-circle.y_,caloPosInTrk.x()-circle.x_);
      int   nloop   = round((caloPhi-(caloPosInTrk.z()-circle.fitb_zp_)/circle.fita_zp_)/6.283185);
      float phiLoop = caloPhi - nloop*6.283185;
      float delta = abs((caloPosInTrk.z()-circle.fitb_zp_)/circle.fita_zp_-phiLoop);
      if (delta < MaxDPhiHelInit_) zphiFitter.add(phiLoop,caloPosInTrk.z(),ccWeight_);
    }

    float fa = zphiFitter.fa();
    float fb = zphiFitter.fb();
    if ((helicity==Helicity::poshel && fa<1e-3) || (helicity==Helicity::neghel  && fa>-1e-3)) return;

    circle.fita_zp_ = fa;
    circle.fitb_zp_ = fb - int(fb/fa/6.283185)*6.283185*fa; //adjust the intercept to have minimum at z=0
    return;
  }

  //---------------------------------------------------------------------------------------------------------------------------
  // Perform dz/dt linear fit including calo cluster if present
  void RobustMultiHelixFinder::fit_dzdt(CandHelix& circle, const ComboHitCollection& chcol)
  {
    LSFitter tzFitter;
    for (const auto& ich : circle.hitIndices_) tzFitter.add(chcol[ich].pos().z(),chcol[ich].correctedTime());

    if (circle.caloPtr_ && circle.caloPtr_->energyDep()>ccMinEnergy_){
      const auto& caloCluster    = *(circle.caloPtr_);
      const auto  posInMu2e      = cal_->geomUtil().diskFFToMu2e(caloCluster.diskID(),caloCluster.cog3Vector());
      const auto  caloClusterPos = cal_->geomUtil().mu2eToTracker(posInMu2e);
      tzFitter.add(caloClusterPos.z(),caloCluster.time(),ccWeight_);
    }

    circle.fita_zt_= tzFitter.fa();
    circle.fitb_zt_= tzFitter.fb();
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Filter hits in z-phi plane. Need to go backwards since erase invalidates forward iterator
  int RobustMultiHelixFinder::filterZPhi(CandHelix& circle, const ComboHitCollection& chcol, float maxDphi)
  {
    int nRemoved(0);
    auto it = circle.hitIndices_.size();
    while (it>0){
       --it;
       unsigned ich = circle.hitIndices_[it];
       float phi    = polyAtan2(chcol[ich].pos().y() - circle.y_,chcol[ich].pos().x() - circle.x_);
       int   n      = round((phi-(chcol[ich].pos().z()-circle.fitb_zp_)/circle.fita_zp_)/6.293185);
       float delta  = abs((chcol[ich].pos().z()-circle.fitb_zp_)/circle.fita_zp_-phi+n*6.293185);
       if (delta < maxDphi) continue;

       circle.remove(it,  chcol[ich].nStrawHits());
       ++nRemoved;
    }
    return nRemoved;
  }

  //---------------------------------------------------------------------------------------------------------------------------
  // Filter hits in z-t plane. Need to go backwards since erase invalidates forward iterator
  int RobustMultiHelixFinder::filterZT(CandHelix& circle, const ComboHitCollection& chcol, float maxDt)
  {
    int nRemoved(0);
    auto it = circle.hitIndices_.size();
    while (it>0){
       --it;
       unsigned ich = circle.hitIndices_[it];
       float delta  = abs(circle.fita_zt_*chcol[ich].pos().z()+circle.fitb_zt_ - chcol[ich].correctedTime());
       if (delta < maxDt) continue;

       circle.remove(it, chcol[ich].nStrawHits());
       ++nRemoved;
    }
    return nRemoved;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Algebraic weighted circle fit using weights from hit resolution (adapted from HyperFit from N. Chernov)
  void RobustMultiHelixFinder::fitCircleAlg(const CandHelix& circle, const ComboHitCollection& chcol, float& centerX,
                                            float& centerY, float& radius)
  {
    std::vector<float> hitWeight(chcol.size(),0);
    for (const auto& ich : circle.hitIndices_) hitWeight[ich] = weight(chcol[ich],circle);//*chcol[ich].nStrawHits();


    float sumWeights(0),xmean(0),ymean(0);
    for (const auto& ich : circle.hitIndices_) {
      xmean      += chcol[ich].pos().x()*hitWeight[ich];
      ymean      += chcol[ich].pos().y()*hitWeight[ich];
      sumWeights += hitWeight[ich];
    }
    xmean /= sumWeights;
    ymean /= sumWeights;

    float Mxx(0),Mxy(0),Myy(0),Mxz(0),Myz(0),Mzz(0);
    for (const auto& ich : circle.hitIndices_){
        float xx = chcol[ich].pos().x()-xmean;
        float yy = chcol[ich].pos().y()-ymean;
        float zz = xx*xx+yy*yy;
        float w  = hitWeight[ich];
        Mxy += xx*yy*w;
        Mxx += xx*xx*w;
        Myy += yy*yy*w;
        Mxz += xx*zz*w;
        Myz += yy*zz*w;
        Mzz += zz*zz*w;
    }
    Mxx /= sumWeights;
    Myy /= sumWeights;
    Mxy /= sumWeights;
    Mxz /= sumWeights;
    Myz /= sumWeights;
    Mzz /= sumWeights;

    float Mz     = Mxx + Myy;
    float Cov_xy = Mxx*Myy - Mxy*Mxy;
    float Var_z  = Mzz - Mz*Mz;
    float A2     = 4.0*Cov_xy - 3.0*Mz*Mz - Mzz;
    float A1     = Var_z*Mz + 4.0*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    float A0     = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    float A22    = A2+A2;

    unsigned iterMAX = 10;
    float x(0.0),y(A0);
    for (unsigned iter=0; iter<iterMAX; ++iter)  {
        float Dy = A1 + x*(A22 + 16.*x*x);
        float xnew = x - y/Dy;
        if (abs(xnew-x) < 1e-4 || abs(x) > 1e10) break;
        float ynew = A0 + xnew*(A1 + xnew*(A2 + 4.0*xnew*xnew));
        if (abs(ynew)> abs(y)) break;
        x = xnew;  y = ynew;
    }

    float DET = x*x - x*Mz + Cov_xy;
    float centerXMean = (Mxz*(Myy - x) - Myz*Mxy)/DET/2.0;
    float centerYMean = (Myz*(Mxx - x) - Mxz*Mxy)/DET/2.0;

    centerX = centerXMean + xmean;
    centerY = centerYMean + ymean;
    radius  = sqrt(centerXMean*centerXMean + centerYMean*centerYMean + Mz - x - x);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Circle fit minimizing chi2 = (delta R2)^2
  void RobustMultiHelixFinder::fitCircleChi2(const CandHelix& circle, const ComboHitCollection& chcol, float& centerX,
                                             float& centerY, float& radius)
  {
    float w_(0),x_(0),x2_(0),x3_(0),y_(0),y2_(0),y3_(0),xy_(0),x2y_(0),xy2_(0);

    for (const auto& ich : circle.hitIndices_) {
      float w = weight(chcol[ich],circle);
      float x = chcol[ich].pos().x();
      float y = chcol[ich].pos().y();
       x_   += x*w;
       x2_  += x*x*w;
       x3_  += x*x*x*w;
       y_   += y*w;
       y2_  += y*y*w;
       y3_  += y*y*y*w;
       xy_  += x*y*w;
       x2y_ += x*x*y*w;
       xy2_ += x*y*y*w;
       w_   += w;
    }

    float d11 = w_*xy_  - x_*y_;
    float d20 = w_*x2_  - x_*x_;
    float d30 = w_*x3_  - x2_*x_;
    float d21 = w_*x2y_ - x2_*y_;
    float d02 = w_*y2_  - y_*y_;
    float d03 = w_*y3_  - y2_*y_;
    float d12 = w_*xy2_ - x_*y2_;
    float den = 2*(d20*d02-d11*d11);
    centerX  = ((d30+d12)*d02-(d03+d21)*d11)/den;
    centerY  = ((d03+d21)*d20-(d30+d12)*d11)/den;
    float c  = (x2_+y2_-2*centerX*x_-2*centerY*y_)/w_;
    radius   = sqrt(c+centerX*centerX+centerY*centerY);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate weight based on circle fit and straw resolution
  float RobustMultiHelixFinder::weight(const ComboHit& ch, const CandHelix& circle)
  {
    float dx     = ch.pos().x() - circle.x_;
    float dy     = ch.pos().y() - circle.y_;
    float rwdot  = (ch.uDir2D().x()*dx + ch.uDir2D().y()*dy);
    float costh2 = rwdot*rwdot/(dx*dx+dy*dy);
    float werr   = ch.posRes(StrawHitPosition::wire);
    float terr   = ch.posRes(StrawHitPosition::trans);
    float res2   = werr*werr*costh2 + terr*terr*(1.0-costh2);
    return 1.0/res2;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate chi2 with full helix hypothesis
  float RobustMultiHelixFinder::chi2XYHelix(const ComboHit& ch, const CandHelix& helix)
  {
    return chi2XYHelix(ch,helix.x_,helix.y_,helix.r_,helix.fita_zp_,helix.fitb_zp_);
  }

  float RobustMultiHelixFinder::chi2XYHelix(const ComboHit& ch,float centerX,float centerY,float radius,float dzdp,float fzb)
  {
    float wx     = ch.uDir2D().x();
    float wy     = ch.uDir2D().y();
    float phiAtZ = (ch.pos().z()-fzb) / dzdp;
    float dhx    = ch.pos().x() - (centerX + radius*cos(phiAtZ));
    float dhy    = ch.pos().y() - (centerY + radius*sin(phiAtZ));
    float dtrans = fabs(-wy*dhx + wx*dhy);
    float dwire  = fabs(wx*dhx + wy*dhy);
    float werr   = ch.posRes(StrawHitPosition::wire);
    float terr   = ch.posRes(StrawHitPosition::trans);
    float wres2  = werr*werr;
    float tres2  = terr*terr;
    float chisq  = dwire*dwire/wres2 + dtrans*dtrans/tres2;
    return chisq;

    /*
    This is a fast version of this
    const XYZVectorF& wdir = hit->uDir2D();
    XYZVectorF wtdir = zaxis.Cross(wdir);   // transverse direction to the wire

    XYZVectorF hpos = hit->pos(); // this sets the z position to the hit z
    helix.position(hpos);                // this computes the helix expectation at that z
    XYZVectorF dh = hit->pos() - hpos;   // this is the vector between them
    float dtrans = fabs(dh.Dot(wtdir)); // transverse projection
    float dwire = fabs(dh.Dot(wdir));   // projection along wire direction

    // compute the total resolution including hit and helix parameters first along the wire
    float wres2 = std::pow(hit->posRes(StrawHitPosition::wire),(int)2)
    // transverse to the wires
    float wtres2 = std::pow(hit->posRes(StrawHitPosition::trans),(int)2)
    float chisq = dwire*dwire/wres2 + dtrans*dtrans/wtres2;
    */
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Calculate chi2 based on circle fit
  float RobustMultiHelixFinder::chi2XYCircle(const ComboHit& ch, float centerX, float centerY, float radius)
  {
    float dx     = ch.pos().x() - centerX;
    float dy     = ch.pos().y() - centerY;
    float r      = sqrt(dx*dx+dy*dy);
    float dr     = r - radius;
    float rwdot  = (ch.uDir2D().x()*dx + ch.uDir2D().y()*dy)/r;
    float rwdot2 = rwdot*rwdot;
    float werr   = ch.posRes(StrawHitPosition::wire);
    float terr   = ch.posRes(StrawHitPosition::trans);
    float rres2  = werr*werr*rwdot2 + terr*terr*(1.0-rwdot2);
    float chisq  = dr*dr/rres2;
    return chisq;
  }



  //---------------------------------------------------------------------------------------------------------------------------
  // Estimate time uncertainty - variance of hit with simple nHit weight
  float RobustMultiHelixFinder::timeUncertainty(const ComboHitCollection& chcol, const CandHelix& helix)
  {
    float st(0),st2(0),sw(0);
    for (const auto& ich : helix.hitIndices_){
      st  += chcol[ich].correctedTime()*chcol[ich].nStrawHits();
      st2 += chcol[ich].correctedTime()*chcol[ich].correctedTime()*chcol[ich].nStrawHits();
      sw  += chcol[ich].nStrawHits();
    }
    st /=sw;
    st2/=sw;
    return sqrt(st2-st*st);
  }


  //---------------------------------------------------------------------------------------------------------------------------
  // Print me!
  void RobustMultiHelixFinder::printHelix(const HelixSeed& helix)
  {
    std::cout<<"Helicity   "<<Helicity::name(helix.helix().helicity())<<std::endl;
    std::cout<<"Radius     "<<helix.helix().radius()<<std::endl;
    std::cout<<"Rcent      "<<helix.helix().rcent()<<std::endl;
    std::cout<<"Fcent      "<<helix.helix().fcent()<<std::endl;
    std::cout<<"Lambda     "<<helix.helix().lambda()<<std::endl;
    std::cout<<"fz0        "<<helix.helix().fz0()<<std::endl;
    std::cout<<"Time       "<<helix.t0().t0()<<"  +- "<<helix.t0().t0Err()<<std::endl;
    std::cout<<"Hits       "<<helix.hits().size()<<std::endl;
    std::cout<<"SH indx(0) ";
    for (const auto& ch : helix.hits()) std::cout<<ch.index(0)<<" ";
    std::cout<<std::endl;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Diagnosis
  void RobustMultiHelixFinder::fillDiag(mapHelix& helcols, const ComboHitCollection& chcol)
  {
    data_.iev_ = iev_;

    data_.Nch_ = chcol.size();
    for (unsigned ich=0; ich<chcol.size();++ich)
    {
       data_.chTime_[ich] = chcol[ich].correctedTime();
       data_.chX_[ich]    = chcol[ich].pos().x();
       data_.chY_[ich]    = chcol[ich].pos().y();
       data_.chZ_[ich]    = chcol[ich].pos().z();
       data_.chRad_[ich]  = sqrt(data_.chX_[ich]*data_.chX_[ich]+data_.chY_[ich]*data_.chY_[ich]);
       data_.chPhi_[ich]  = chcol[ich].pos().phi();
       data_.chNhit_[ich] = chcol[ich].nStrawHits();
       data_.chUId_[ich]  = chcol[ich].strawId().uniquePanel();
       data_.chTerr_[ich] = chcol[ich].posRes(StrawHitPosition::trans);
       data_.chWerr_[ich] = chcol[ich].posRes(StrawHitPosition::wire);
       data_.chWDX_[ich]  = chcol[ich].uDir2D().x();
       data_.chWDY_[ich]  = chcol[ich].uDir2D().y();

       std::vector<StrawHitIndex> strawHitIdxs;
       chcol.fillStrawHitIndices(ich,strawHitIdxs);
       data_.chStrawId_.push_back(strawHitIdxs);
    }

    for (const auto& helix : *helcols[Helicity::poshel]){
      std::vector<int> hhits;
      for (auto hi : helix.hits()){
        for (size_t j=0;j<chcol.size();++j){
          if (hi.index(0)==chcol[j].index(0)) {hhits.push_back(j); break;}
        }
      }
      data_.helhel_[data_.Nhel_]  = 1;
      data_.helrad_[data_.Nhel_]  = helix.helix().radius();
      data_.helrcen_[data_.Nhel_] = helix.helix().rcent();
      data_.helfcen_[data_.Nhel_] = helix.helix().fcent();
      data_.hellam_[data_.Nhel_]  = helix.helix().lambda();
      data_.helfz0_[data_.Nhel_]  = helix.helix().fz0();
      data_.helchi2_[data_.Nhel_] = helix.helix().chi2dZPhi();
      data_.heldzdt_[data_.Nhel_] = helix.helix().chi2dXY();
      data_.helnhi_[data_.Nhel_]  = hhits.size();
      data_.helhits_.push_back(hhits);
      ++data_.Nhel_;
   }

   for (const auto& helix : *helcols[Helicity::neghel]){
      std::vector<int> hhits;
      for (auto hi : helix.hits()){
        for (size_t j=0;j<chcol.size();++j){
          if (hi.index(0)==chcol[j].index(0)) {hhits.push_back(j); break;}
        }
      }
      data_.helhel_[data_.Nhel_]  = -1;
      data_.helrad_[data_.Nhel_]  = helix.helix().radius();
      data_.helrcen_[data_.Nhel_] = helix.helix().rcent();
      data_.helfcen_[data_.Nhel_] = helix.helix().fcent();
      data_.hellam_[data_.Nhel_]  = helix.helix().lambda();
      data_.helfz0_[data_.Nhel_]  = helix.helix().fz0();
      data_.heldzdt_[data_.Nhel_] = helix.helix().chi2dXY();
      data_.helnhi_[data_.Nhel_]  = hhits.size();
      data_.helchi2_[data_.Nhel_] = helix.helix().chi2dZPhi();
      data_.helhits_.push_back(hhits);
      ++data_.Nhel_;
    }
  }


}

DEFINE_ART_MODULE(mu2e::RobustMultiHelixFinder);
