///////////////////////////////////////////////////////////////////////////////
// TZClusterFinder
// M. Stortini
///////////////////////////////////////////////////////////////////////////////

#include "Offline/CalPatRec/inc/TZClusterFinder_types.hh"

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "art/Utilities/make_tool.h"

#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include <vector>
#include <TH1.h>
#include <TF1.h>
#include <TROOT.h>
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"

using CLHEP::HepVector;
using CLHEP::HepSymMatrix;
using CLHEP::Hep3Vector;

namespace mu2e {

  using namespace TZClusterFinderTypes;

  class TZClusterFinder: public art::EDProducer {

  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>               diagLevel        {Name("diagLevel"        ), Comment("turn tool on or off"         ) };
      fhicl::Atom<int>               debugLevel       {Name("debugLevel"       ), Comment("turn on/off debug"           ) };
      fhicl::Atom<int>               printFrequency   {Name("printFrequency"   ), Comment("print frequency"             ) };
      fhicl::Atom<int>               runDisplay       {Name("runDisplay"       ), Comment("will plot t vs z"            ) };
      fhicl::Atom<int>               saveProtCand     {Name("saveProtCand"     ), Comment("save proton candidate TCs"   ) };
      fhicl::Atom<int>               useCaloClusters  {Name("useCaloClusters"  ), Comment("use calo cluster option"     ) };
      fhicl::Atom<int>               recoverCaloClusters {Name("recoverCaloClusters"), Comment("recover TCs using CCs"  ) };
      fhicl::Atom<art::InputTag>     ComboHitCollectionLabel {Name("ComboHitCollectionLabel"), Comment("chColl"         ) };
      fhicl::Atom<art::InputTag>     ComboHitCollectionLabel2 {Name("ComboHitCollectionLabel2"), Comment("for MC tool"  ) };
      fhicl::Atom<art::InputTag>     TimeClusterCollectionLabel {Name("TimeClusterCollectionLabel"), Comment("tcColl"   ) };
      fhicl::Atom<art::InputTag>     StrawHitFlagCollection {Name("StrawHitFlagCollection"), Comment("shfColl"          ) };
      fhicl::Atom<art::InputTag>     caloClusterModuleLabel {Name("caloClusterModuleLabel"), Comment("CC module label"  ) };
      fhicl::Sequence<std::string>   HitBackgroundBits {Name("HitBackgroundBits"      ), Comment("background bits"      ) };
      fhicl::Atom<int>               chunkSep         {Name("chunkSep"         ), Comment("max num of planes for chunk" ) };
      fhicl::Atom<double>            chunkWindow      {Name("chunkWindow"      ), Comment("time window in ns"           ) };
      fhicl::Atom<int>               chunkThresh      {Name("chunkThresh"      ), Comment("number of combo hits"        ) };
      fhicl::Atom<double>            chi2combineThresh {Name("chi2combineThresh"), Comment("thresh for combining chks"  ) };
      fhicl::Atom<double>            combineWindow    {Name("combineWindow"    ), Comment("time window in ns"           ) };
      fhicl::Atom<double>            maxCombineSep    {Name("maxCombineSep"    ), Comment("z range in mm"               ) };
      fhicl::Atom<int>               chunkFitThresh   {Name("chunkFitThresh"   ), Comment("number of combo hits"        ) };
      fhicl::Atom<double>            recoverWindow    {Name("recoverWindow"    ), Comment("time window in ns"           ) };
      fhicl::Atom<int>               clusterThresh    {Name("clusterThresh"    ), Comment("number of straw hits"        ) };
      fhicl::Atom<int>               doComptonClean   {Name("doComptonClean"   ), Comment("whether or not to do clean"  ) };
      fhicl::Atom<double>            maxIntersectSigma {Name("maxIntersectSigma"), Comment("for compton clean"          ) };
      fhicl::Atom<double>            maxApproachSigma {Name("maxApproachSigma"), Comment("for compton clean"            ) };
      fhicl::Atom<double>            maxApproachSigmaTrans {Name("maxApproachSigmaTrans"), Comment("for compton clean"  ) };
      fhicl::Atom<int>               comptonThresh    {Name("comptonThresh"    ), Comment("number of combo hits"        ) };
      fhicl::Atom<int>               doIsoClean       {Name("doIsoClean"    ), Comment("whether or not to do clean"     ) };
      fhicl::Atom<double>            isoRad           {Name("isoRad"), Comment("radial search distance (mm) for clean"  ) };
      fhicl::Atom<int>               doPhiClean       {Name("doPhiClean"    ), Comment("whether or not to do clean"     ) };
      fhicl::Atom<double>            phiCleanPhi      {Name("phiCleanPhi"), Comment("in degrees"                        ) };
      fhicl::Atom<int>               phiCleanThresh   {Name("phiCleanThresh"    ), Comment("number of combo hits"       ) };
      fhicl::Atom<int>               minCaloSize      {Name("minCaloSize"    ), Comment("number of crystals"            ) };
      fhicl::Atom<double>            minCaloEnergy    {Name("minCaloEnergy"), Comment("in MeV"                          ) };
      fhicl::Atom<double>            caloDtMax        {Name("caloDtMax"), Comment("search time window (ns)"             ) };
      fhicl::Atom<double>            caloTimeOffset   {Name("caloTimeOffset"), Comment("in ns"                          ) };

      fhicl::Table<TZClusterFinderTypes::Config> diagPlugin{Name("diagPlugin"      ), Comment("Diag plugin") };
    };

  private:

    //-----------------------------------------------------------------------------
    // data members
    //-----------------------------------------------------------------------------
    unsigned         _iev;
    int              _diagLevel;
    int              _debugLevel;
    int              _printfreq;
    int              _runDisplay;
    int              _saveProtCand;
    int              _useCaloClusters;
    int              _recoverCaloClusters;

    //-----------------------------------------------------------------------------
    // event object labels
    //-----------------------------------------------------------------------------
    art::InputTag   _chLabel ;
    art::InputTag   _chLabel2;
    art::InputTag   _tcLabel ;
    art::InputTag   _shfLabel;
    art::InputTag   _ccmLabel;
    StrawHitFlag    _hbkg;

    //-----------------------------------------------------------------------------
    // cluster search parameters
    //-----------------------------------------------------------------------------
    int      _chunkSep; // number of planes we allow chunks to be combined in
    double   _chunkWindow; // time window we allow hits to live within to be chunked together
    int      _chunkThresh; // how many hits need to be in a chunk to be saved
    double   _chi2combineThresh; // value below which chi2dof is considered valid for combining chunks
    double   _combineWindow; // time window in which chunks combining may be considered
    double   _maxCombineSep; // max z separation to consider combining chunks
    int      _chunkFitThresh; // how many hits chunk must have to do fit to before recovering hits
    double   _recoverWindow; // time hit must be within chunk fit to be added to chunk
    int      _clusterThresh; // number of combo hits needed to save cluster found
    int      _doComptonClean; // whether or not to clean up cluster with compton search
    double   _maxIntersectSigma; // max allowable distance from hit centers to intersection point for seed pair in compton clean
    double   _maxApproachSigma; // max allowable distance from hit center to closest approach to seed intersection in compton clean
    double   _maxApproachSigmaTrans; // like _maxApproachSigma but transverse direction
    int      _comptonThresh; // how many hits need to be consistent with a compton to clean up cluster
    int      _doIsoClean; // whether or not to clean up cluster by removing hits that are too isolated
    double   _isoRad; // size of circle around test hit that must have no other hits in it for test hit to be considered isolated
    int      _doPhiClean; // whether or not to do phi clean up to chunks
    double   _phiCleanPhi; // phi value to use during phi clean up
    int      _phiCleanThresh; // number of hits below which hit is considered too isolated in phi
    int      _minCaloSize; // number of crystals for calo cluster to be considered
    double   _minCaloEnergy; // minimum energy for calo cluster to be considered
    double   _caloDtMax; // max time from time cluster for calo cluster to be associated with time cluster
    double   _caloTimeOffset; // time offset for calo clusters

    //-----------------------------------------------------------------------------
    // diagnostics
    //-----------------------------------------------------------------------------
    Data_t                                 _data;
    art::Handle<CaloClusterCollection>     _ccHandle;
    facilitateVars                         _f;
    std::unique_ptr<ModuleHistToolBase>    _hmanager;
    TCanvas* c1;

    //-----------------------------------------------------------------------------
    // functions
    //-----------------------------------------------------------------------------
  public:

    explicit TZClusterFinder(const art::EDProducer::Table<Config>& config);
    virtual ~TZClusterFinder();

    virtual void beginJob ();
    virtual void beginRun (art::Run&);
    virtual void produce (art::Event& e);
    virtual void endJob   ();

    //-----------------------------------------------------------------------------
    // helper functions
    //-----------------------------------------------------------------------------
    bool findData               (const art::Event& e);
    void cHitsFill              ();
    bool bkgHit                 (const StrawHitFlag& flag) const;
    void plotTZ                 ();
    void setSeed                (int& seedPln, size_t& seedPlnHit);
    void getValidStartIndex     (int& testPln, int& seedPln, size_t& seedPlnHit);
    void setTestHit             (int& testPln, size_t& testPlnHit);
    void testTestHit            (int& testPln, size_t& testPlnHit);
    void flagUsedHits           ();
    void chunkHits              ();
    void combineChunks          ();
    void recoverHits            ();
    void countProtons           (IntensityInfoTimeCluster& outIITC);
    void cleanCompton           ();
    void isoClean               ();
    void phiClean               ();
    void checkCaloClusters      ();
    void findClusters           (TimeClusterCollection& OutSeeds);
  };

  //-----------------------------------------------------------------------------
  // module constructor
  //-----------------------------------------------------------------------------
  TZClusterFinder::TZClusterFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _diagLevel              (config().diagLevel()                               ),
    _debugLevel             (config().debugLevel()                              ),
    _printfreq              (config().printFrequency()                          ),
    _runDisplay             (config().runDisplay()                              ),
    _saveProtCand           (config().saveProtCand()                            ),
    _useCaloClusters        (config().useCaloClusters()                         ),
    _recoverCaloClusters    (config().recoverCaloClusters()                     ),
    _chLabel                (config().ComboHitCollectionLabel()                 ),
    _chLabel2               (config().ComboHitCollectionLabel2()                ),
    _tcLabel                (config().TimeClusterCollectionLabel()              ),
    _shfLabel               (config().StrawHitFlagCollection()                  ),
    _ccmLabel               (config().caloClusterModuleLabel()                  ),
    _hbkg                   (config().HitBackgroundBits()                       ),
    _chunkSep               (config().chunkSep()                                ),
    _chunkWindow            (config().chunkWindow()                             ),
    _chunkThresh            (config().chunkThresh()                             ),
    _chi2combineThresh      (config().chi2combineThresh()                       ),
    _combineWindow          (config().combineWindow()                           ),
    _maxCombineSep          (config().maxCombineSep()                           ),
    _chunkFitThresh         (config().chunkFitThresh()                          ),
    _recoverWindow          (config().recoverWindow()                           ),
    _clusterThresh          (config().clusterThresh()                           ),
    _doComptonClean         (config().doComptonClean()                          ),
    _maxIntersectSigma      (config().maxIntersectSigma()                       ),
    _maxApproachSigma       (config().maxApproachSigma()                        ),
    _maxApproachSigmaTrans  (config().maxApproachSigmaTrans()                   ),
    _comptonThresh          (config().comptonThresh()                           ),
    _doIsoClean             (config().doIsoClean()                              ),
    _isoRad                 (config().isoRad()                                  ),
    _doPhiClean             (config().doPhiClean()                              ),
    _phiCleanPhi            (config().phiCleanPhi()                             ),
    _phiCleanThresh         (config().phiCleanThresh()                          ),
    _minCaloSize            (config().minCaloSize()                             ),
    _minCaloEnergy          (config().minCaloEnergy()                           ),
    _caloDtMax              (config().caloDtMax()                               ),
    _caloTimeOffset         (config().caloTimeOffset()                          )
    {

      consumes<ComboHitCollection>(_chLabel);
      consumes<CaloClusterCollection>(_ccmLabel);
      produces<TimeClusterCollection>();
      produces<IntensityInfoTimeCluster>();


      if (_runDisplay == 1) { c1 = new TCanvas("c1", "t vs. z", 900, 900); }

      if (_debugLevel != 0) _printfreq = 1;

      if (_diagLevel  != 0) _hmanager = art::make_tool<ModuleHistToolBase>(config().diagPlugin, "diagPlugin");
      else                  _hmanager = std::make_unique<ModuleHistToolBase>();

    }

  //-----------------------------------------------------------------------------
  // destructor
  //-----------------------------------------------------------------------------
  TZClusterFinder::~TZClusterFinder() {}

  //-----------------------------------------------------------------------------
  // beginJob
  //-----------------------------------------------------------------------------
  void TZClusterFinder::beginJob(){
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

  //-----------------------------------------------------------------------------
  // beginRun
  //-----------------------------------------------------------------------------
  void TZClusterFinder::beginRun(art::Run& ) {}

  //-----------------------------------------------------------------------------
  // find input things
  //-----------------------------------------------------------------------------
  bool TZClusterFinder::findData(const art::Event& evt) {

    auto chcolH = evt.getValidHandle<ComboHitCollection>(_chLabel);
    if (chcolH.product() != 0){
      _data.chcol = chcolH.product();
    }
    else {
      _data.chcol  = 0;
      std::cout << ">>> ERROR in TZClusterFinder::findData: ComboHitCollection not found." << std::endl;
    }


    if (_diagLevel  != 0) {
      auto chcolH2 = evt.getValidHandle<ComboHitCollection>(_chLabel2);
      if (chcolH2.product() != 0){
        _data.chcol2 = chcolH2.product();
      }
    }

    if (_useCaloClusters == 1) {
      if (evt.getByLabel(_ccmLabel, _ccHandle)) {
        _data.ccCollection = _ccHandle.product();
      }
      else {
        _data.ccCollection = NULL;
      }
    }

    auto shfH = evt.getValidHandle<StrawHitFlagCollection>(_shfLabel);
    if (shfH.product() != 0) { _data.shfcol = shfH.product(); }
    else {
      _data.shfcol = 0;
      std::cout << ">>> ERROR in TZClusterFinder::findData: StrawHitFlagCollection not found." << std::endl;
    }

    if(_data.shfcol->size() != _data.chcol->size())
      throw cet::exception("RECO")<<"TimeClusterFinder: inconsistent flag collection length " << std::endl;

    return (_data.chcol != 0);
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void TZClusterFinder::produce(art::Event& event) {

    const char* oname = "TZClusterFinder::filter";

    // event printout
    _iev = event.id().event();
    if ((_debugLevel > 0) && (_iev%_printfreq) == 0) printf("[%s] : START event number %8i\n", oname,_iev);

    // for when we are in runDisplay mode
    if (_runDisplay == 1) {c1->Clear();}

    _data._event = &event;

    std::unique_ptr<TimeClusterCollection>  tcColl(new TimeClusterCollection);
    std::unique_ptr<IntensityInfoTimeCluster> iiTC(new IntensityInfoTimeCluster);

    _data._tcColl = tcColl.get();
    _data._iiTC = iiTC.get();

    bool ok = findData(event);

    if (ok) { findClusters(*_data._tcColl); }
    else    printf("%s ERROR: No straw hits found in event %i\n",oname,_iev);

    if (_diagLevel > 0) { _hmanager->fillHistograms(&_data);}

    //-----------------------------------------------------------------------------
    // put time cluster collection and protons counted into the event record
    //-----------------------------------------------------------------------------
    event.put(std::move(tcColl));
    event.put(std::move(iiTC));

  }

  //-----------------------------------------------------------------------------
  // endJob
  //-----------------------------------------------------------------------------
  void TZClusterFinder::endJob(){ }

  //-----------------------------------------------------------------------------
  // order combo hits to aid in finding clusters
  //-----------------------------------------------------------------------------
  void TZClusterFinder::cHitsFill() {

    const mu2e::ComboHit* hit;
    // fill cHits indexed by pln, each column being a vector housing cHit info
    for (size_t i=0; i<_data.chcol->size(); i++) {
      if ((*_data.shfcol)[i].hasAnyProperty(StrawHitFlag::energysel)) { if (bkgHit((*_data.shfcol)[i])) {continue;} }
      hit = &_data.chcol->at(i);
      int plnID = hit->strawId().plane();
      cHit comboHit;
      comboHit.hIndex = i;
      comboHit.hTime = hit->correctedTime();
      comboHit.hWeight = 1/(hit->timeVar());
      comboHit.hZpos = hit->pos().z();
      comboHit.nStrawHits = hit->nStrawHits();
      comboHit.hIsUsed = 0;
      _f.cHits[plnID].plnHits.push_back(comboHit);
    }

    // order the plnHits in each pln by descending time
    for (size_t i=0; i<_f.cHits.size(); i++) {
      std::sort(_f.cHits[i].plnHits.begin(), _f.cHits[i].plnHits.end(),
                [](const cHit &a, const cHit &b) -> bool
                { return a.hTime > b.hTime;});
    }
  }


  //-----------------------------------------------------------------------------
  // for testing if hit is flagged as background
  //-----------------------------------------------------------------------------
  bool TZClusterFinder::bkgHit(const StrawHitFlag& flag) const {

    return flag.hasAnyProperty(_hbkg);

  }

  //-----------------------------------------------------------------------------
  // if in runDisplay mode we plot T vs. Z
  //-----------------------------------------------------------------------------
  void TZClusterFinder::plotTZ() {

    if (_runDisplay == 1) {
      const mu2e::ComboHit* hit;
      // set up multigraph that graphs can be added to
      TMultiGraph* mg = new TMultiGraph();
      mg->GetXaxis()->SetLimits(-1600,1600);
      mg->GetYaxis()->SetRangeUser(0,1800);
      mg->GetXaxis()->SetTitle("z (mm)");
      mg->GetYaxis()->SetTitle("t (ns)");
      // create graph for t vs. z
      int dPoints = 0;
      for (size_t i=0; i<_f.cHits.size(); i++) { dPoints = dPoints + _f.cHits[i].plnHits.size(); }
      float zPoints[dPoints];
      float tPoints[dPoints];
      int arrayIndex = 0;
      for (size_t i=0; i<_f.cHits.size(); i++) {
        for (size_t j=0; j<_f.cHits[i].plnHits.size(); j++) {
          zPoints[arrayIndex] = _f.cHits[i].plnHits[j].hZpos;
          tPoints[arrayIndex] = _f.cHits[i].plnHits[j].hTime;
          arrayIndex++;
        }
      }
      if (dPoints > 0 ) {
        TGraph* g1 = new TGraph(dPoints, zPoints, tPoints);
        g1->SetMarkerStyle(20);
        g1->SetMarkerSize(1.1);
        mg->Add(g1);
        TGraph* clusterPlots[_data._nTZClusters];
        TF1* fit = new TF1("fit", "pol1");
        for (int i=0; i<_data._nTZClusters; i++) {
          int dcPoints = (int)_data._tcColl->at(i)._strawHitIdxs.size();
          double_t zcPoints[dcPoints];
          double_t tcPoints[dcPoints];
          for (int j=0; j<dcPoints; j++) {
            int index = (int)_data._tcColl->at(i)._strawHitIdxs[j];
            hit = &_data.chcol->at(index);
            double_t zPosition = hit->pos().z();
            double_t hitTime = hit->correctedTime();
            zcPoints[j] = zPosition;
            tcPoints[j] = hitTime;
          }
          clusterPlots[i] = new TGraph(dcPoints, zcPoints, tcPoints);
          clusterPlots[i]->SetMarkerStyle(20);
          clusterPlots[i]->SetMarkerSize(0.65);
          if (i>7) {
            clusterPlots[i]->SetMarkerColor(i+20);
            fit->SetLineColor(i+20);
          }
          else {
            clusterPlots[i]->SetMarkerColor(i+2);
            fit->SetLineColor(i+2);
          }
          mg->Add(clusterPlots[i]);
          fit->FixParameter(0,_data.lineIntercept[i]);
          fit->FixParameter(1,_data.lineSlope[i]);
          clusterPlots[i]->Fit("fit","B");
        }
        // draw graphs
        c1->cd();
        mg->Draw("AP");
        c1->Update();
      }
    }

  }

  //-----------------------------------------------------------------------------
  // for setting original seed in findClusters logic
  //----------------------------------------------------------------------------
  void TZClusterFinder::setSeed(int& seedPln, size_t& seedPlnHit) {

    _f.clear_chunkInfo();
    _f.holdIndices.clear(); // hold index in cHits where seeds are so you can flag them if cluster is saved
    _f.seedIndice = _f.cHits[seedPln].plnHits[seedPlnHit].hIndex;
    _f.seedTime = _f.cHits[seedPln].plnHits[seedPlnHit].hTime;
    _f.seedWeight = _f.cHits[seedPln].plnHits[seedPlnHit].hWeight;
    _f.seedZpos = _f.cHits[seedPln].plnHits[seedPlnHit].hZpos;
    if ((*_data.shfcol)[_f.seedIndice].hasAnyProperty(StrawHitFlag::energysel)) { _f.seedNRGselection = 1; }
    else { _f.seedNRGselection = 0; }
    _f.indicePair.first = seedPln;
    _f.indicePair.second = seedPlnHit;
    _f.holdIndices.push_back(_f.indicePair);
    _f._chunkInfo.hIndices.push_back(_f.seedIndice);
    _f._chunkInfo.fitter.addPoint(_f.seedZpos, _f.seedTime, _f.seedWeight);
    _f.nHitsInChunk = 1;
    _f.nStrawHitsInChunk =  _f.cHits[seedPln].plnHits[seedPlnHit].nStrawHits;
    _f.totalTime = _f.seedTime;
    _f.totalZpos = _f.seedZpos;

  }


  //-----------------------------------------------------------------------------
  // after setting test plane we need to set where to start our loop over hits
  //-----------------------------------------------------------------------------
  void TZClusterFinder::getValidStartIndex (int& testPln, int& seedPln, size_t& seedPlnHit) {

    if (testPln == seedPln) {
      if (seedPlnHit == _f.cHits[seedPln].plnHits.size()-1) {_f.startIndex=-1;}
      else {_f.startIndex = seedPlnHit+1;}
    }
    else { _f.startIndex = 0; }

  }

  //-----------------------------------------------------------------------------
  // for setting test hit in findClusters logic
  //-----------------------------------------------------------------------------
  void TZClusterFinder::setTestHit(int& testPln, size_t& testPlnHit) {

    _f.testIndice = _f.cHits[testPln].plnHits[testPlnHit].hIndex;
    _f.testTime = _f.cHits[testPln].plnHits[testPlnHit].hTime;
    _f.testWeight = _f.cHits[testPln].plnHits[testPlnHit].hWeight;
    _f.testZpos = _f.cHits[testPln].plnHits[testPlnHit].hZpos;;
    if ((*_data.shfcol)[_f.testIndice].hasAnyProperty(StrawHitFlag::energysel)) { _f.testNRGselection = 1; }
    else { _f.testNRGselection = 0; }

  }

  //-----------------------------------------------------------------------------
  // put in findclusters logic to test if hit is compatible with seed
  //-----------------------------------------------------------------------------
  void TZClusterFinder::testTestHit(int& testPln, size_t& testPlnHit) {

    if (std::abs(_f.testTime - _f.seedTime) < _chunkWindow) {
      _f.indicePair.first = testPln;
      _f.indicePair.second = testPlnHit;
      _f.holdIndices.push_back(_f.indicePair);
      _f._chunkInfo.hIndices.push_back(_f.testIndice);
      _f._chunkInfo.fitter.addPoint(_f.testZpos, _f.testTime, _f.testWeight);
      _f.nHitsInChunk++;
      _f.nStrawHitsInChunk = _f.nStrawHitsInChunk + _f.cHits[testPln].plnHits[testPlnHit].nStrawHits;
      _f.totalTime += _f.testTime;
      _f.totalZpos += _f.testZpos;
    }

  }

  //-----------------------------------------------------------------------------
  // for flagging hits used in saved cluster
  //-----------------------------------------------------------------------------
  void TZClusterFinder::flagUsedHits() {

    // flag hits used so they only are used in one cluster
    for (size_t r=0; r<_f.holdIndices.size(); r++) {
      _f.cHits[_f.holdIndices[r].first].plnHits[_f.holdIndices[r].second].hIsUsed = 1;
    }

  }

  //-----------------------------------------------------------------------------
  // logic for looping over cHits to create chunks of hits
  //-----------------------------------------------------------------------------
  void TZClusterFinder::chunkHits() {

    // first two for loops create seed point
    for (int i=(int)_f.cHits.size()-1; i>=0; i--) {;
      for (size_t j=0; j<_f.cHits[i].plnHits.size(); j++) {
        if ( _f.cHits[i].plnHits[j].hIsUsed != 0 ) {continue;}
        setSeed(i,j);
        // now we find points around the seed
        for (int k=i; k>=0; k--) {
          if ( std::abs(k-i) > _chunkSep ) { break; }
          getValidStartIndex(k,i,j);
          if (_f.startIndex == -1) { continue; }
          size_t start = (size_t)_f.startIndex;
          for (size_t n=start; n<_f.cHits[k].plnHits.size(); n++) {
            if ( _f.cHits[k].plnHits[n].hIsUsed != 0 ) {continue;}
            setTestHit(k,n);
            if (_f.seedNRGselection != _f.testNRGselection) {continue;}
            testTestHit(k,n);
          }
        }
        // if we exceed threshold we save chunk of hits
        if ((int)_f._chunkInfo.hIndices.size()>=_chunkThresh) {
          _f._chunkInfo.avgTime = _f.totalTime/_f.nHitsInChunk;
          _f._chunkInfo.avgZpos = _f.totalZpos/_f.nHitsInChunk;
          _f._chunkInfo.nHits = _f.nHitsInChunk;
          _f._chunkInfo.nStrawHits = _f.nStrawHitsInChunk;
          _f._chunkInfo.nrgSelection = _f.seedNRGselection;
          _f._chunkInfo.nCombines = 0;
          _f._chunkInfo.caloIndex = -1;
          _f.chunks.push_back(_f._chunkInfo);
          flagUsedHits();
        }
      }
    }

  }


  //-----------------------------------------------------------------------------
  // logic for combining chunks of hits together
  //-----------------------------------------------------------------------------
  void TZClusterFinder::combineChunks() {

    _f.biggestChi2combine = 1000.;

    // set seed chunk
    for (size_t i=0; i<_f.chunks.size()-1; i++) {
      // add seed to fitter
      _f.chi2seed = _f.chunks[i].fitter.chi2Dof();
      _f.seedTime = _f.chunks[i].avgTime;
      _f.seedZpos = _f.chunks[i].avgZpos;
      // now set test chunk
      for (size_t k=i+1; k<_f.chunks.size(); k++) {
        if (_f.chunks[i].nrgSelection != _f.chunks[k].nrgSelection) {continue;}
        _f.testTime = _f.chunks[k].avgTime;
        _f.testZpos = _f.chunks[k].avgZpos;
        if (std::abs(_f.seedTime-_f.testTime) >= _combineWindow) {continue;}
        if (_f.chunks[i].nCombines == 0 && _f.chunks[k].nCombines == 0
            && std::abs(_f.seedZpos-_f.testZpos) > _maxCombineSep) {continue;}
        // add test to fitter
        _f.chunks[i].fitter.addSum(_f.chunks[k].fitter);
        _f.chi2combineTest = _f.chunks[i].fitter.chi2Dof();
        if (_f.chi2combineTest < _chi2combineThresh) {
          if (_f.chi2combineTest < _f.biggestChi2combine) {
            _f.biggestChi2combine = _f.chi2combineTest;
            _f.chunkOneIdx = i;
            _f.chunkTwoIdx = k;
          }
        }
        // remove test from fitter
        _f.chunks[i].fitter.removeSum(_f.chunks[k].fitter);
      }
    }

    if (_f.biggestChi2combine < 1000.) {
      for (size_t i=0; i<_f.chunks[_f.chunkTwoIdx].hIndices.size(); i++) {
        _f.chunks[_f.chunkOneIdx].hIndices.push_back(_f.chunks[_f.chunkTwoIdx].hIndices[i]);
      }
      _f.chunks[_f.chunkOneIdx].fitter.addSum(_f.chunks[_f.chunkTwoIdx].fitter);
      _f.chunks[_f.chunkOneIdx].avgTime = _f.chunks[_f.chunkOneIdx].avgTime * _f.chunks[_f.chunkOneIdx].nHits;
      _f.chunks[_f.chunkOneIdx].avgTime += _f.chunks[_f.chunkTwoIdx].avgTime * _f.chunks[_f.chunkTwoIdx].nHits;
      _f.chunks[_f.chunkOneIdx].avgTime /= _f.chunks[_f.chunkOneIdx].nHits + _f.chunks[_f.chunkTwoIdx].nHits;
      _f.chunks[_f.chunkOneIdx].avgZpos = _f.chunks[_f.chunkOneIdx].avgZpos * _f.chunks[_f.chunkOneIdx].nHits;
      _f.chunks[_f.chunkOneIdx].avgZpos += _f.chunks[_f.chunkTwoIdx].avgZpos * _f.chunks[_f.chunkTwoIdx].nHits;
      _f.chunks[_f.chunkOneIdx].avgZpos /= _f.chunks[_f.chunkOneIdx].nHits + _f.chunks[_f.chunkTwoIdx].nHits;
      _f.chunks[_f.chunkOneIdx].nHits = _f.chunks[_f.chunkOneIdx].nHits + _f.chunks[_f.chunkTwoIdx].nHits;
      _f.chunks[_f.chunkOneIdx].nStrawHits = _f.chunks[_f.chunkOneIdx].nStrawHits + _f.chunks[_f.chunkTwoIdx].nStrawHits;
      _f.chunks[_f.chunkOneIdx].nCombines++;
      _f.chunks.erase(_f.chunks.begin()+_f.chunkTwoIdx);
    }

    else { _f.moreCombines = false; }

  }

  //-----------------------------------------------------------------------------
  // for recovering hits
  //-----------------------------------------------------------------------------
  void TZClusterFinder::recoverHits() {

    double smallestDeltaT;
    double deltaTtest;
    int strawhits = 0;
    int chunkIndex = 0;

    for (int i=(int)_f.cHits.size()-1; i>=0; i--) {
      for (size_t j=0; j<_f.cHits[i].plnHits.size(); j++) {
        if (_f.cHits[i].plnHits[j].hIsUsed == 1) {continue;}
        _f.testTime = _f.cHits[i].plnHits[j].hTime;
        _f.testWeight = _f.cHits[i].plnHits[j].hWeight;
        _f.testZpos = _f.cHits[i].plnHits[j].hZpos;
        _f.testIndice = _f.cHits[i].plnHits[j].hIndex;
        strawhits = _f.cHits[i].plnHits[j].nStrawHits;
        if ((*_data.shfcol)[_f.testIndice].hasAnyProperty(StrawHitFlag::energysel)) { _f.testNRGselection = 1; }
        else { _f.testNRGselection = 0; }
        smallestDeltaT = 1000.;
        for (size_t k=0; k<_f.chunks.size(); k++) {
          if ((int)_f.chunks[k].hIndices.size() < _chunkFitThresh) {continue;}
          if (_f.testNRGselection != _f.chunks[k].nrgSelection) {continue;}
          deltaTtest = std::abs(_f.testTime - (_f.chunks[k].fitter.dydx()*_f.testZpos + _f.chunks[k].fitter.y0()));
          if (deltaTtest < smallestDeltaT) {
            smallestDeltaT = deltaTtest;
            chunkIndex = k;
          }
        }
        if (smallestDeltaT <= _recoverWindow) {
          _f.chunks[chunkIndex].hIndices.push_back(_f.testIndice);
          _f.chunks[chunkIndex].nStrawHits = _f.chunks[chunkIndex].nStrawHits + strawhits;
          _f.chunks[chunkIndex].nHits++;
          _f.chunks[chunkIndex].fitter.addPoint(_f.testZpos, _f.testTime, _f.testWeight);
        }
      }
    }

  }

  //-----------------------------------------------------------------------------
  // function for counting protons (prediction not truth)
  //-----------------------------------------------------------------------------
  void TZClusterFinder::countProtons(IntensityInfoTimeCluster& outIITC) {

    unsigned short nProtons = 0;

    for (size_t i=0; i<_f.chunks.size(); i++) {
      if (_f.chunks[i].nrgSelection == 1) {continue;}
      if (_f.chunks[i].nStrawHits < _clusterThresh) {continue;}
      nProtons++;
    }

    outIITC.setNProtonTCs(nProtons);

  }

  //-----------------------------------------------------------------------------
  // function for removing compton consistent hits from CE-like chunks
  //-----------------------------------------------------------------------------
  void TZClusterFinder::cleanCompton() {

    // vectors to specify indices of the chunk vector
    std::vector<int> tempIndices;
    std::vector<int> comptonIndices;

    // combo hit pointer
    const mu2e::ComboHit* hit;

    for (size_t i=0; i<_f.chunks.size(); i++) { // clean comptons in each chunk
      comptonIndices.clear();
      tempIndices.clear();
      if (_f.chunks[i].nrgSelection == 0 || _f.chunks[i].nStrawHits < _clusterThresh) {continue;}
      for (size_t j=0; j<_f.chunks[i].hIndices.size(); j++) { // set first point for seed pair
        hit = &_data.chcol->at(_f.chunks[i].hIndices[j]);
        XYZVectorF h1pos = hit->pos(); // hit position
        const XYZVectorF& h1wdir = hit->wdir(); // wire direction of hit
        float h1wsig = hit->wireRes(); // hit resolution along wire
        //float h1wtsig = hit->transRes(); // hit resolution transverse to wire
        // now set up linear line equation for wire
        float m1 = h1wdir.y()/h1wdir.x();
        float b1 = h1pos.y() - m1*h1pos.x();
        for (size_t k=0; k<_f.chunks[i].hIndices.size(); k++) { // set second point for seed pair
          if (k==j) {continue;}
          tempIndices.clear();
          hit = &_data.chcol->at(_f.chunks[i].hIndices[k]);
          XYZVectorF h2pos = hit->pos();
          const XYZVectorF& h2wdir = hit->wdir();
          float h2wsig = hit->wireRes();
          //float h2wtsig = hit->transRes();
          float m2 = h2wdir.y()/h2wdir.x();
          float b2 = h2pos.y() - m2*h2pos.x();
          if (std::fabs(m1-m2)/m1 < 0.01) {continue;} // dont consider hits whose wires are parallel
          // find where wire error bars intersect
          float xInt = (b2-b1)/(m1-m2);
          float yInt = m1*xInt+b1;
          // continue to new point if intersection is too far from pair points
          float dh1 = std::sqrt((h1pos.x()-xInt)*(h1pos.x()-xInt)+(h1pos.y()-yInt)*(h1pos.y()-yInt))/h1wsig;
          float dh2 = std::sqrt((h2pos.x()-xInt)*(h2pos.x()-xInt)+(h2pos.y()-yInt)*(h2pos.y()-yInt))/h2wsig;
          if (dh1 > _maxIntersectSigma || dh2 > _maxIntersectSigma) {continue;}
          tempIndices.push_back(j);
          tempIndices.push_back(k);
          for (size_t q=0; q<_f.chunks[i].hIndices.size(); q++) { // test points with intersection found
            if (q==j || q==k) {continue;}
            hit = &_data.chcol->at(_f.chunks[i].hIndices[q]);
            XYZVectorF h3pos = hit->pos();
            const XYZVectorF& h3wdir = hit->wdir();
            float h3wsig = hit->wireRes();
            float h3wtsig = hit->transRes();
            float m3 = h3wdir.y()/h3wdir.x();
            float b3 = h3pos.y() - m3*h3pos.x();
            // find closest approach
            float m = -1.0/m3;
            float b = yInt - m*xInt;
            float xClosest = (b-b3)/(m3-m);
            float yClosest = m*xClosest+b;
            float dApproach = std::sqrt((xClosest-xInt)*(xClosest-xInt)+(yClosest-yInt)*(yClosest-yInt))/h3wtsig;
            if (dApproach > _maxApproachSigmaTrans) {continue;}
            float dh3 = std::sqrt((h3pos.x()-xClosest)*(h3pos.x()-xClosest)+(h3pos.y()-yClosest)*(h3pos.y()-yClosest))/h3wsig;
            if (dh3 > _maxApproachSigma) {continue;}
            tempIndices.push_back(q);
          }
          if (tempIndices.size() > comptonIndices.size()) {
            comptonIndices.clear();
            for (size_t n=0; n<tempIndices.size(); n++) {
              comptonIndices.push_back(tempIndices[n]);
            }
          }
        }
      }
      if ((int)comptonIndices.size() > _comptonThresh) {
        for (size_t m=0; m<comptonIndices.size(); m++) {
          int chunkIndex = comptonIndices[m];
          hit = &_data.chcol->at(_f.chunks[i].hIndices[chunkIndex]);
          _f.chunks[i].nStrawHits = _f.chunks[i].nStrawHits - hit->nStrawHits();
          _f.chunks[i].nHits = _f.chunks[i].nHits - 1;
          _f.chunks[i].hIndices[chunkIndex] = -1;
        }
        _f.chunks[i].hIndices.erase(std::remove(_f.chunks[i].hIndices.begin(), _f.chunks[i].hIndices.end(), -1), _f.chunks[i].hIndices.end());
      }
    }

  }

  //-----------------------------------------------------------------------------
  // logic to do isolated hits clean up
  //-----------------------------------------------------------------------------
  void TZClusterFinder::isoClean() {

    // combo hit pointers
    const mu2e::ComboHit* hit;
    const mu2e::ComboHit* test;

    // hit position variables
    double hitX = 0.0;
    double hitY = 0.0;
    double testX = 0.0;
    double testY = 0.0;
    double separation = 0.0;
    int isolated = 1;


    for (size_t i=0; i<_f.chunks.size(); i++) { // loop over chunks
      if (_f.chunks[i].nrgSelection == 0 || _f.chunks[i].nStrawHits < _clusterThresh) {continue;}
      for (size_t j=0; j<_f.chunks[i].hIndices.size(); j++) { // set hit to check if isolated
        hit = &_data.chcol->at(_f.chunks[i].hIndices[j]);
        hitX = hit->pos().x();
        hitY = hit->pos().y();
        isolated = 1;
        for (size_t k=0; k<_f.chunks[i].hIndices.size(); k++) { // look for hits near hit being tested
          if (k==j || _f.chunks[i].hIndices[k] == -1 ) {continue;}
          test =  &_data.chcol->at(_f.chunks[i].hIndices[k]);
          testX = test->pos().x();
          testY = test->pos().y();
          separation = std::sqrt((hitX-testX)*(hitX-testX)+(hitY-testY)*(hitY-testY));
          if (separation <= _isoRad) {
            isolated = 0;
            break;
          }
        }
        if (isolated == 1) {
          _f.chunks[i].hIndices[j] = -1;
          _f.chunks[i].nHits = _f.chunks[i].nHits - 1;
          _f.chunks[i].nStrawHits = _f.chunks[i].nStrawHits - hit->nStrawHits();

        }
      }
      _f.chunks[i].hIndices.erase(std::remove(_f.chunks[i].hIndices.begin(), _f.chunks[i].hIndices.end(), -1), _f.chunks[i].hIndices.end());
    }

  }

  //-----------------------------------------------------------------------------
  // logic to do phi clean up
  //-----------------------------------------------------------------------------
  void TZClusterFinder::phiClean() {

    // combo hit pointers
    const mu2e::ComboHit* hit;
    const mu2e::ComboHit* test;

    // hit position variables
    double hitPhi = 0.0;
    double testPhi = 0.0;
    double deltaPhi = 0.0;

    for (size_t i=0; i<_f.chunks.size(); i++) { // loop over chunks
      if (_f.chunks[i].nrgSelection == 0 || _f.chunks[i].nStrawHits < _clusterThresh) {continue;}
      for (size_t j=0; j<_f.chunks[i].hIndices.size(); j++) { // set hit to check
        hit = &_data.chcol->at(_f.chunks[i].hIndices[j]);
        hitPhi = hit->phi();
        int hitsInRange = 1;
        for (size_t k=0; k<_f.chunks[i].hIndices.size(); k++) { // test other hits
          if (k==j || _f.chunks[i].hIndices[k] == -1) {continue;}
          test =  &_data.chcol->at(_f.chunks[i].hIndices[k]);
          testPhi = test->phi();
          deltaPhi = std::fmod(testPhi-hitPhi, CLHEP::twopi);
          if (deltaPhi>CLHEP::pi) deltaPhi -= CLHEP::twopi;
          if (deltaPhi<-CLHEP::pi) deltaPhi += CLHEP::twopi;
          deltaPhi = std::abs(deltaPhi)*180.0/CLHEP::pi;
          if (deltaPhi <= _phiCleanPhi) {
            hitsInRange++;
          }
          if (hitsInRange >= _phiCleanThresh) {break;}
        }
        if (hitsInRange < _phiCleanThresh) {
          _f.chunks[i].hIndices[j] = -1;
          _f.chunks[i].nHits = _f.chunks[i].nHits - 1;
          _f.chunks[i].nStrawHits = _f.chunks[i].nStrawHits - hit->nStrawHits();
        }
      }
      _f.chunks[i].hIndices.erase(std::remove(_f.chunks[i].hIndices.begin(), _f.chunks[i].hIndices.end(), -1), _f.chunks[i].hIndices.end());
    }

  }


  //-----------------------------------------------------------------------------
  // logic to use calo clusters
  //-----------------------------------------------------------------------------
  void TZClusterFinder::checkCaloClusters() {

    const CaloCluster* cc;
    const ComboHit* hit;
    int ncc = _data.ccCollection->size();
    int nchunks = _f.chunks.size();
    double ccTime = 0.0;
    int addedToTC = 0;

    for (int i=0; i<ncc; i++) {
      cc = &_data.ccCollection->at(i);
      if (cc->energyDep() < _minCaloEnergy || cc->size() < _minCaloSize) {continue;}
      ccTime = cc->time() + _caloTimeOffset;
      addedToTC = 0;
      for (int j=0; j<nchunks; j++) {
        if ((int)_f.chunks[j].hIndices.size() < _chunkFitThresh) {continue;}
        double dT = std::abs((double)_f.chunks[j].fitter.y0() - ccTime);
        if (dT < _caloDtMax) {
          if (_f.chunks[j].caloIndex != -1) {
            _f._chunkInfo = _f.chunks[j];
            _f._chunkInfo.caloIndex = i;
            _f.chunks.push_back(_f._chunkInfo);
          }
          else {
            _f.chunks[j].caloIndex = i;
          }
          addedToTC = 1;
        }
      }
      if (_recoverCaloClusters == 1) {
        if (addedToTC == 0) {
          _f.clear_chunkInfo();
          _f._chunkInfo.nHits = 0;
          _f._chunkInfo.nStrawHits = 0;
          _f._chunkInfo.caloIndex = i;
          for (size_t k=0; k<_data.chcol->size(); k++) {
            if (bkgHit((*_data.shfcol)[k])) {continue;}
            if (!(*_data.shfcol)[k].hasAnyProperty(StrawHitFlag::energysel)) {continue;}
            hit = &_data.chcol->at(k);
            if (std::abs(hit->correctedTime() - ccTime) < _caloDtMax) {
              _f._chunkInfo.hIndices.push_back(k);
              _f._chunkInfo.fitter.addPoint(hit->pos().z(), hit->correctedTime(), 1/(hit->timeVar()));
              _f._chunkInfo.nHits++;
              _f._chunkInfo.nStrawHits = _f._chunkInfo.nStrawHits + hit->nStrawHits();
            }
          }
          if (_f._chunkInfo.nStrawHits >= _clusterThresh) {
            _f.chunks.push_back(_f._chunkInfo);
          }
        }
      }
    }

  }

  //-----------------------------------------------------------------------------
  // logic to find clusters
  //-----------------------------------------------------------------------------
  void TZClusterFinder::findClusters(TimeClusterCollection& TimeClusterColl) {

    // clear relevant data members that carry over from processing of previous event
    _data.clearRelevant();
    _f.clear_cHits();
    _f.clear_chunks();
    _f._clusterInfo._caloCluster = art::Ptr<mu2e::CaloCluster>();

    // fill cHits array to loop over
    cHitsFill();

    // loop over cHits looking for hits that can be chunked together
    chunkHits();

    // look for chunks that can be combined together
    if (_f.chunks.size() != 0) {
      _f.moreCombines = true;
      while(_f.moreCombines) { combineChunks(); }
    }

    // recover hits that were missed
    recoverHits();

    // count number of protons
    countProtons(*_data._iiTC);

    // use calo clusters
    if (_useCaloClusters == 1) { checkCaloClusters(); }

    // do clean ups
    if (_doComptonClean == 1) { cleanCompton(); }
    if (_doIsoClean == 1) { isoClean();}
    if (_doPhiClean == 1) { phiClean();}

    // fill TimeClusterColl and perform final fit to save
    for (size_t i=0; i<_f.chunks.size(); i++) {
      if ((int)_f.chunks[i].nStrawHits < _clusterThresh) {continue;}
      if (_saveProtCand != 1 && _f.chunks[i].nrgSelection == 0) {continue;}
      _f.clear_clusterInfo();
      for (size_t j=0; j<_f.chunks[i].hIndices.size(); j++) {
        _f._clusterInfo._strawHitIdxs.push_back(StrawHitIndex(_f.chunks[i].hIndices[j]));
      }
      _f._clusterInfo._t0 = TrkT0(_f.chunks[i].fitter.y0(), 0.);
      int caloIdx = _f.chunks[i].caloIndex;
      if (caloIdx != -1) {
        _f._clusterInfo._caloCluster = art::Ptr<mu2e::CaloCluster>(_ccHandle, caloIdx);
      }
      else {
        _f._clusterInfo._caloCluster = art::Ptr<mu2e::CaloCluster>();
      }
      TimeClusterColl.push_back(_f._clusterInfo);
      _data.lineSlope.push_back(_f.chunks[i].fitter.dydx());
      _data.lineIntercept.push_back(_f.chunks[i].fitter.y0());
      _data.chi2DOF.push_back(_f.chunks[i].fitter.chi2Dof());
    }


    _data._nTZClusters = (int)TimeClusterColl.size();

    plotTZ();

  }

}

using mu2e::TZClusterFinder;
DEFINE_ART_MODULE(TZClusterFinder);
