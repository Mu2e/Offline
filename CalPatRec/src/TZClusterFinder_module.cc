///////////////////////////////////////////////////////////////////////////////
// TZClusterFinder
// M. Stortini and A. M. Ricci
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

namespace mu2e {

  using namespace TZClusterFinderTypes;

  class TZClusterFinder: public art::EDProducer {

  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>               diagLevel          {Name("diagLevel"          ), Comment("turn tool on or off"                             ) };
      fhicl::Atom<int>               debugLevel         {Name("debugLevel"         ), Comment("turn on/off debug"                               ) };
      fhicl::Atom<int>               printFrequency     {Name("printFrequency"     ), Comment("print frequency"                                 ) };
      fhicl::Atom<int>               runDisplay         {Name("runDisplay"         ), Comment("will plot t vs z"                                ) };
      fhicl::Atom<int>               useCCs             {Name("useCCs"             ), Comment("add CCs to TCs"                                  ) };
      fhicl::Atom<int>               recoverCCs         {Name("recoverCCs"         ), Comment("recover TCs using CCs"                           ) };
      fhicl::Atom<art::InputTag>     chCollLabel        {Name("chCollLabel"        ), Comment("combo hit collection label"                      ) };
      fhicl::Atom<art::InputTag>     chCollLabel2       {Name("chCollLabel2"       ), Comment("for MC tool"                                     ) };
      fhicl::Atom<art::InputTag>     tcCollLabel        {Name("tcCollLabel"        ), Comment("time cluster coll label"                         ) };
      fhicl::Atom<art::InputTag>     ccCollLabel        {Name("ccCollLabel"        ), Comment("Calo Cluster coll label"                         ) };
      fhicl::Sequence<std::string>   hitBkgBits         {Name("hitBkgBits"         ), Comment("background bits"                                 ) };
      fhicl::Atom<int>               radSelect          {Name("radSelect"          ), Comment("whether or not to radial cut"                    ) };
      fhicl::Atom<int>               chunkSep           {Name("chunkSep"           ), Comment("max # of planes for chunk"                       ) };
      fhicl::Atom<double>            chunkWindow        {Name("chunkWindow"        ), Comment("time window in ns"                               ) };
      fhicl::Atom<int>               chunkThresh        {Name("chunkThresh"        ), Comment("number of straw hits"                            ) };
      fhicl::Atom<double>            combineWindow      {Name("combineWindow"      ), Comment("time window in ns"                               ) };
      fhicl::Atom<double>            maxCombineSep      {Name("maxCombineSep"      ), Comment("z range in mm"                                   ) };
      fhicl::Atom<int>               chunkFitThresh     {Name("chunkFitThresh"     ), Comment("number of combo hits"                            ) };
      fhicl::Atom<double>            recoverWindow      {Name("recoverWindow"      ), Comment("time window in ns"                               ) };
      fhicl::Atom<int>               clusterThresh      {Name("clusterThresh"      ), Comment("number of straw hits"                            ) };
      fhicl::Atom<int>               minCaloSize        {Name("minCaloSize"        ), Comment("number of crystals"                              ) };
      fhicl::Atom<double>            minCaloEnergy      {Name("minCaloEnergy"      ), Comment("in MeV"                                          ) };
      fhicl::Atom<double>            caloDtMax          {Name("caloDtMax"          ), Comment("search time window (ns)"                         ) };
      fhicl::Atom<double>            caloTimeOffset     {Name("caloTimeOffset"     ), Comment("in ns"                                           ) };
      fhicl::Atom<int>               doRefine           {Name("doRefine"           ), Comment("filter out bad TCs at end"                       ) };
      fhicl::Atom<size_t>            minSHsInCluster    {Name("minSHsInCluster"    ), Comment("min number of straw hits in time cluster"        ) };
      fhicl::Atom<size_t>            minCHsInCluster    {Name("minCHsInCluster"    ), Comment("min number of combo hits in time cluster"        ) };
      fhicl::Atom<size_t>            nPairs             {Name("nPairs"             ), Comment("min number of pairs of neighboring stations"     ) };
      fhicl::Atom<size_t>            minCHsInCubes      {Name("minCHsInCubes"      ), Comment("min number of combo hits in the cubes dR-dPhi-dZ") };
      fhicl::Atom<double>            maxDeltaRInStn     {Name("maxDeltaRInStn"     ), Comment("max deltaR (mm) in the station"                  ) };
      fhicl::Atom<double>            maxDeltaRBtwStns   {Name("maxDeltaRBtwStns"   ), Comment("max deltaR (mm) between stations"                ) };
      fhicl::Atom<double>            maxDeltaPhiInStn   {Name("maxDeltaPhiInStn"   ), Comment("max deltaPhi in the station"                     ) };
      fhicl::Atom<double>            maxDeltaPhiBtwStns {Name("maxDeltaPhiBtwStns" ), Comment("max deltaPhi between stations"                   ) };
      fhicl::Atom<bool>              thirdStation       {Name("thirdStation"       ), Comment("check the stations near the pairs"               ) };

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
    int              _useCaloClusters;
    int              _recoverCaloClusters;

    //-----------------------------------------------------------------------------
    // event object labels
    //-----------------------------------------------------------------------------
    art::InputTag   _chLabel ;
    art::InputTag   _chLabel2;
    art::InputTag   _tcLabel ;
    art::InputTag   _ccLabel;
    StrawHitFlag    _hbkg;
    int             _radSelect;

    //-----------------------------------------------------------------------------
    // cluster search parameters
    //-----------------------------------------------------------------------------
    int      _chunkSep; // number of planes we allow chunks to be combined in
    double   _chunkWindow; // time window we allow hits to live within to be chunked together
    int      _chunkThresh; // how many hits need to be in a chunk to be saved
    double   _combineWindow; // time window in which chunks combining may be considered
    double   _maxCombineSep; // max z separation to consider combining chunks
    int      _chunkFitThresh; // how many hits chunk must have to do fit to before recovering hits
    double   _recoverWindow; // time hit must be within chunk fit to be added to chunk
    int      _clusterThresh; // number of straw hits needed to save cluster found
    int      _minCaloSize; // number of crystals for calo cluster to be considered
    double   _minCaloEnergy; // minimum energy for calo cluster to be considered
    double   _caloDtMax; // max time from time cluster for calo cluster to be associated with time cluster
    double   _caloTimeOffset; // time offset for calo clusters
    int      _doRefine; // if set to 1 then pattern recognition is used to filter out bad TC candidates
    int      _minSHsInCluster; // minimum number of straw hits in TC to pass Refine selection
    int      _minCHsInCluster; // minimum number of combo hits in TC to pass Refine selection
    size_t   _nPairs; // minimum number of pairs of neighboring stations
    size_t   _minCHs; // minimum number of combo hits in the cubes dR-dPhi-dZ
    double   _deltaR1; // maximum deltaR (mm) in the station
    double   _deltaR2; // maximum deltaR (mm) between stations
    double   _deltaPhi_1; // maximum deltaPhi in the station
    double   _deltaPhi_2; // maximum deltaPhi between stations
    bool     _thirdStation; // check the stations near the pairs

    //-----------------------------------------------------------------------------
    // diagnostics
    //-----------------------------------------------------------------------------
    Data_t                                 _data;
    art::Handle<CaloClusterCollection>     _ccHandle;
    facilitateVars                         _f;
    std::unique_ptr<ModuleHistToolBase>    _hmanager;
    TCanvas*                               _c1;

    //-----------------------------------------------------------------------------
    // functions
    //-----------------------------------------------------------------------------

  public:

    explicit TZClusterFinder(const art::EDProducer::Table<Config>& config);
    virtual ~TZClusterFinder();

    virtual void beginJob ();
    virtual void beginRun (art::Run&);
    virtual void produce  (art::Event& e);
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
    void checkCaloClusters      ();
    void refineChunks           ();
    void findClusters           (TimeClusterCollection& OutSeeds);
    void fillClusterHits        (size_t tc);
    bool checkPopulation        (size_t tc);
    bool checkPattern           (size_t tc);
    bool checkStation           (size_t stn);
    bool checkPair              (size_t stn1, size_t stn2);
    bool checkThirdStation      (size_t stn1, size_t stn2);
    void clusterInfo            (size_t tc);
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
    _useCaloClusters        (config().useCCs()                                  ),
    _recoverCaloClusters    (config().recoverCCs()                              ),
    _chLabel                (config().chCollLabel()                             ),
    _chLabel2               (config().chCollLabel2()                            ),
    _tcLabel                (config().tcCollLabel()                             ),
    _ccLabel                (config().ccCollLabel()                             ),
    _hbkg                   (config().hitBkgBits()                              ),
    _radSelect              (config().radSelect()                               ),
    _chunkSep               (config().chunkSep()                                ),
    _chunkWindow            (config().chunkWindow()                             ),
    _chunkThresh            (config().chunkThresh()                             ),
    _combineWindow          (config().combineWindow()                           ),
    _maxCombineSep          (config().maxCombineSep()                           ),
    _chunkFitThresh         (config().chunkFitThresh()                          ),
    _recoverWindow          (config().recoverWindow()                           ),
    _clusterThresh          (config().clusterThresh()                           ),
    _minCaloSize            (config().minCaloSize()                             ),
    _minCaloEnergy          (config().minCaloEnergy()                           ),
    _caloDtMax              (config().caloDtMax()                               ),
    _caloTimeOffset         (config().caloTimeOffset()                          ),
    _doRefine               (config().doRefine()                                ),
    _minSHsInCluster        (config().minSHsInCluster()                         ),
    _minCHsInCluster        (config().minCHsInCluster()                         ),
    _nPairs                 (config().nPairs()                                  ),
    _minCHs                 (config().minCHsInCubes()                           ),
    _deltaR1                (config().maxDeltaRInStn()                          ),
    _deltaR2                (config().maxDeltaRBtwStns()                        ),
    _deltaPhi_1             (config().maxDeltaPhiInStn()                        ),
    _deltaPhi_2             (config().maxDeltaPhiBtwStns()                      ),
    _thirdStation           (config().thirdStation()                            )
    {

      consumes<ComboHitCollection>(_chLabel);
      consumes<CaloClusterCollection>(_ccLabel);
      produces<TimeClusterCollection>();
      produces<IntensityInfoTimeCluster>();


      if (_runDisplay == 1) { _c1 = new TCanvas("_c1", "t vs. z", 900, 900); }

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
      _data._chColl = chcolH.product();
    }
    else {
      _data._chColl  = 0;
      std::cout << ">>> ERROR in TZClusterFinder::findData: ComboHitCollection not found." << std::endl;
    }


    if (_diagLevel  != 0) {
      auto chcolH2 = evt.getValidHandle<ComboHitCollection>(_chLabel2);
      if (chcolH2.product() != 0){
        _data._chColl2 = chcolH2.product();
      }
    }

    if (_useCaloClusters == 1) {
      if (evt.getByLabel(_ccLabel, _ccHandle)) {
        _data._ccColl = _ccHandle.product();
      }
      else {
        _data._ccColl = NULL;
      }
    }

    return (_data._chColl != 0);
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
    if (_runDisplay == 1) { _c1->Clear(); }

    _data._event = &event;

    std::unique_ptr<IntensityInfoTimeCluster> iiTC(new IntensityInfoTimeCluster);
    std::unique_ptr<TimeClusterCollection>    tcColl(new TimeClusterCollection);

    _data._tcColl = tcColl.get();
    _data._iiTC = iiTC.get();

    bool ok = findData(event);

    if (ok) { findClusters(*_data._tcColl); }
    else    { printf("%s ERROR: No straw hits found in event %i\n",oname,_iev); }

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
    for (size_t i=0; i<_data._chColl->size(); i++) {
      const StrawHitFlag flag = _data._chColl->at(i).flag();
      if (!flag.hasAnyProperty(StrawHitFlag::radsel) && _radSelect == 1) {continue;}
      if (flag.hasAnyProperty(StrawHitFlag::energysel)) { if (bkgHit(flag)) {continue;} }
      hit = &_data._chColl->at(i);
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
  // order combo hits of the TC candidates to aid the selection
  //-----------------------------------------------------------------------------
  void TZClusterFinder::fillClusterHits(size_t tc) {

    _f.clear_chStn();

    cHit comboHit;

    // fill chStn indexed by stn, each column being a vector housing cHit info
    for(size_t i=0; i<_f.chunks[tc].hIndices.size(); i++) {
      int hitIndex = _f.chunks[tc].hIndices[i];
      int stnID = _data._chColl->at(hitIndex).strawId().station();
      comboHit.center = false;
      comboHit.hIndex = hitIndex;
      comboHit.nStrawHits = _data._chColl->at(hitIndex).nStrawHits();
      comboHit.x = _data._chColl->at(hitIndex).pos().x();
      comboHit.y = _data._chColl->at(hitIndex).pos().y();
      comboHit.z = _data._chColl->at(hitIndex).pos().z();
      comboHit.phi = _data._chColl->at(hitIndex).phi();
      _f.chStn[stnID].stnHits.push_back(comboHit);
    }

    // sort the vector in ascending order of the z-coordinate
    for(size_t i=0; i<_f.chStn.size(); i++) {
      std::sort(_f.chStn[i].stnHits.begin(), _f.chStn[i].stnHits.end(),
                [](const cHit& a, const cHit& b) -> bool { return a.z < b.z; } );
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
    std::vector<float> zPoints(dPoints);
    std::vector<float> tPoints(dPoints);
    int arrayIndex = 0;
    for (size_t i=0; i<_f.cHits.size(); i++) {
      for (size_t j=0; j<_f.cHits[i].plnHits.size(); j++) {
        zPoints[arrayIndex] = _f.cHits[i].plnHits[j].hZpos;
        tPoints[arrayIndex] = _f.cHits[i].plnHits[j].hTime;
        arrayIndex++;
      }
    }
    if (dPoints > 0 ) {
      TGraph* g1 = new TGraph(dPoints, zPoints.data(), tPoints.data());
      g1->SetMarkerStyle(20);
      g1->SetMarkerSize(1.1);
      mg->Add(g1);
      std::vector<TGraph*> clusterPlots(_data._nTZClusters);
      TF1* fit = new TF1("fit", "pol1");
      for (int i=0; i<_data._nTZClusters; i++) {
        int dcPoints = (int)_data._tcColl->at(i)._strawHitIdxs.size();
        std::vector<double_t> zcPoints(dcPoints);
        std::vector<double_t> tcPoints(dcPoints);
        for (int j=0; j<dcPoints; j++) {
          int index = (int)_data._tcColl->at(i)._strawHitIdxs[j];
          hit = &_data._chColl->at(index);
          double_t zPosition = hit->pos().z();
          double_t hitTime = hit->correctedTime();
          zcPoints[j] = zPosition;
          tcPoints[j] = hitTime;
        }
        clusterPlots[i] = new TGraph(dcPoints, zcPoints.data(), tcPoints.data());
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
      _c1->cd();
      mg->Draw("AP");
      _c1->Update();
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
    const StrawHitFlag flag = _data._chColl->at(_f.seedIndice).flag();
    if (flag.hasAnyProperty(StrawHitFlag::energysel)) { _f.seedNRGselection = 1; }
    else { _f.seedNRGselection = 0; }
    _f._indicePair.first = seedPln;
    _f._indicePair.second = seedPlnHit;
    _f.holdIndices.push_back(_f._indicePair);
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
    _f.testZpos = _f.cHits[testPln].plnHits[testPlnHit].hZpos;
    const StrawHitFlag flag = _data._chColl->at(_f.testIndice).flag();
    if (flag.hasAnyProperty(StrawHitFlag::energysel)) { _f.testNRGselection = 1; }
    else { _f.testNRGselection = 0; }

  }

  //-----------------------------------------------------------------------------
  // put in findClusters logic to test if hit is compatible with seed
  //-----------------------------------------------------------------------------
  void TZClusterFinder::testTestHit(int& testPln, size_t& testPlnHit) {

    if (std::abs(_f.testTime - _f.seedTime) < _chunkWindow) {
      _f._indicePair.first = testPln;
      _f._indicePair.second = testPlnHit;
      _f.holdIndices.push_back(_f._indicePair);
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

    size_t nProtonChunks = 0;
    size_t nCeLikeChunks = 0;

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
          _f._chunkInfo.goodCluster = true;
          _f.chunks.push_back(_f._chunkInfo);
          if (_f._chunkInfo.nrgSelection == 0) { nProtonChunks++;}
          else { nCeLikeChunks++; }
          flagUsedHits();
        }
      }
    }

    // separate chunks by nrgSelection, CE first then protons, then order by ascending time
    if (nProtonChunks != 0 && nCeLikeChunks != 0) {
      std::sort(_f.chunks.begin(), _f.chunks.end(),
                [](const chunkInfo &a, const chunkInfo &b) -> bool
                { return a.nrgSelection > b.nrgSelection;});
      std::sort(_f.chunks.begin(), _f.chunks.begin()+nCeLikeChunks,
                [](const chunkInfo &a, const chunkInfo &b) -> bool
                { return a.avgTime < b.avgTime;});
      std::sort(_f.chunks.begin()+nCeLikeChunks, _f.chunks.end(),
                [](const chunkInfo &a, const chunkInfo &b) -> bool
                { return a.avgTime < b.avgTime;});
    }
    else {
      std::sort(_f.chunks.begin(), _f.chunks.end(),
                [](const chunkInfo &a, const chunkInfo &b) -> bool
                { return a.avgTime < b.avgTime;});
    }

  }


  //-----------------------------------------------------------------------------
  // logic for combining chunks of hits together
  //-----------------------------------------------------------------------------
  void TZClusterFinder::combineChunks() {

    double minValidDtFound    = 0.0;
    size_t validCombinesFound = 0;
    size_t chunkOneIdx        = 0;
    size_t chunkTwoIdx        = 0;

    for (size_t i=0; i<_f.chunks.size()-1; i++) {
      if (_f.chunks[i].nrgSelection != _f.chunks[i+1].nrgSelection) {continue;}
      _f.seedZpos = _f.chunks[i].avgZpos;
      _f.testZpos = _f.chunks[i+1].avgZpos;
      if (_f.chunks[i].nCombines == 0 && _f.chunks[i+1].nCombines == 0
          && std::abs(_f.seedZpos-_f.testZpos) > _maxCombineSep) {continue;}
      _f.seedTime = _f.chunks[i].avgTime;
      _f.testTime = _f.chunks[i+1].avgTime;
      double deltaTime = std::abs(_f.seedTime - _f.testTime);
      if (deltaTime > _combineWindow) {
        if (validCombinesFound != 0) {break;}
        else {continue;}
      }
      validCombinesFound++;
      if (validCombinesFound == 1 || deltaTime < minValidDtFound) {
        minValidDtFound = deltaTime;
        chunkOneIdx = i;
        chunkTwoIdx = i+1;
      }
    }

    if (validCombinesFound != 0) {
      for (size_t i=0; i<_f.chunks[chunkTwoIdx].hIndices.size(); i++) {
        _f.chunks[chunkOneIdx].hIndices.push_back(_f.chunks[chunkTwoIdx].hIndices[i]);
      }
      _f.chunks[chunkOneIdx].fitter.addSum(_f.chunks[chunkTwoIdx].fitter);
      _f.chunks[chunkOneIdx].avgTime = _f.chunks[chunkOneIdx].avgTime * _f.chunks[chunkOneIdx].nHits;
      _f.chunks[chunkOneIdx].avgTime += _f.chunks[chunkTwoIdx].avgTime * _f.chunks[chunkTwoIdx].nHits;
      _f.chunks[chunkOneIdx].avgTime /= _f.chunks[chunkOneIdx].nHits + _f.chunks[chunkTwoIdx].nHits;
      _f.chunks[chunkOneIdx].avgZpos = _f.chunks[chunkOneIdx].avgZpos * _f.chunks[chunkOneIdx].nHits;
      _f.chunks[chunkOneIdx].avgZpos += _f.chunks[chunkTwoIdx].avgZpos * _f.chunks[chunkTwoIdx].nHits;
      _f.chunks[chunkOneIdx].avgZpos /= _f.chunks[chunkOneIdx].nHits + _f.chunks[chunkTwoIdx].nHits;
      _f.chunks[chunkOneIdx].nHits = _f.chunks[chunkOneIdx].nHits + _f.chunks[chunkTwoIdx].nHits;
      _f.chunks[chunkOneIdx].nStrawHits = _f.chunks[chunkOneIdx].nStrawHits + _f.chunks[chunkTwoIdx].nStrawHits;
      _f.chunks[chunkOneIdx].nCombines++;
      _f.chunks.erase(_f.chunks.begin()+chunkTwoIdx);
    }

    else { _f.moreCombines = false; }

  }

  //-----------------------------------------------------------------------------
  // for recovering hits
  //-----------------------------------------------------------------------------
  void TZClusterFinder::recoverHits() {

    double minValidDtFound = 0.0;
    size_t validLinesFound = 0;
    size_t chunkIndex      = 0;

    for (int i=(int)_f.cHits.size()-1; i>=0; i--) {
      for (size_t j=0; j<_f.cHits[i].plnHits.size(); j++) {
        if (_f.cHits[i].plnHits[j].hIsUsed == 1) {continue;}
        _f.testTime = _f.cHits[i].plnHits[j].hTime;
        _f.testWeight = _f.cHits[i].plnHits[j].hWeight;
        _f.testZpos = _f.cHits[i].plnHits[j].hZpos;
        _f.testIndice = _f.cHits[i].plnHits[j].hIndex;
        const StrawHitFlag flag = _data._chColl->at(_f.testIndice).flag();
        if (flag.hasAnyProperty(StrawHitFlag::energysel)) { _f.testNRGselection = 1; }
        else { _f.testNRGselection = 0; }
        validLinesFound = 0;
        minValidDtFound = 0.0;
        for (size_t k=0; k<_f.chunks.size(); k++) {
          if ((int)_f.chunks[k].hIndices.size() < _chunkFitThresh) {continue;}
          if (_f.testNRGselection != _f.chunks[k].nrgSelection) {continue;}
          double deltaTtest = std::abs(_f.testTime - (_f.chunks[k].fitter.dydx()*_f.testZpos + _f.chunks[k].fitter.y0()));
          if (deltaTtest > _recoverWindow) {continue;}
          validLinesFound++;
          if (validLinesFound == 1 || deltaTtest < minValidDtFound) {
            minValidDtFound = deltaTtest;
            chunkIndex = k;
          }
        }
        if (validLinesFound != 0) {
          _f.chunks[chunkIndex].hIndices.push_back(_f.testIndice);
          _f.chunks[chunkIndex].nStrawHits = _f.chunks[chunkIndex].nStrawHits + _f.cHits[i].plnHits[j].nStrawHits;
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

  // logic to use calo clusters
  //-----------------------------------------------------------------------------
  void TZClusterFinder::checkCaloClusters() {

    const CaloCluster* cc;
    const ComboHit*    hit;

    double ccTime    = 0.0;
    int    ncc       = _data._ccColl->size();
    int    nchunks   = _f.chunks.size();
    int    addedToTC = 0;

    for (int i=0; i<ncc; i++) {
      cc = &_data._ccColl->at(i);
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
          for (size_t k=0; k<_data._chColl->size(); k++) {
            const StrawHitFlag flag = _data._chColl->at(k).flag();
            if (!flag.hasAnyProperty(StrawHitFlag::radsel) && _radSelect == 1) {continue;}
            if (bkgHit(flag)) {continue;}
            if (!flag.hasAnyProperty(StrawHitFlag::energysel)) {continue;}
            hit = &_data._chColl->at(k);
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
  // logic to flag bad TC candidates
  //-----------------------------------------------------------------------------
  void TZClusterFinder::refineChunks() {

    if(_diagLevel !=0) {
      _data._nTCRef.clear();
      _data._nChuncks = (int) _f.chunks.size();
    }

    for (size_t i=0; i<_f.chunks.size(); i++) {
      // first continue on chunks that are already not saved
      if (_f.chunks[i].nrgSelection == 0) {continue;}
      if ((int)_f.chunks[i].nStrawHits < _clusterThresh) {continue;}
      for (size_t j=0; j<_f.chunks[i].hIndices.size(); j++) {
        // ethan, put logic here
      }
      // ethan, set _f.chunks[i].goodCluster = false or true here based on result of logic you put above

      // pattern recognition to select the TC candidates containing CE
      fillClusterHits(i);
      if(_debugLevel > 0) clusterInfo(i);
      bool popCheck = checkPopulation(i);
      if(!popCheck) {
        _f.chunks[i].goodCluster = false;
        if(_diagLevel !=0) _data._nTCRef.push_back(0);
        continue;
      }
      bool pattern = checkPattern(i);
      if(!pattern) {
        _f.chunks[i].goodCluster = false;
        if(_diagLevel !=0) _data._nTCRef.push_back(0);
        continue;
      }
      _f.chunks[i].goodCluster = true;
      if(_diagLevel !=0) _data._nTCRef.push_back(1);
    }

  }

  //-----------------------------------------------------------------------------
  // check the population of the TC candidate
  //-----------------------------------------------------------------------------
  bool TZClusterFinder::checkPopulation(size_t tc) {

    // the time cluster is rejected if it has too few straw hits
    if(_f.chunks[tc].nStrawHits < _minSHsInCluster) {
      if(_debugLevel > 0) {
        std::cout << ">>> WARNING in TZClusterFinder::checkPopulation:" << std::endl;
        std::cout << ">>> The Time Cluster " << tc << " contains less than " << _minSHsInCluster <<
                                                            " Straw Hits and is rejected." << std::endl;
      }
      return false;
    }

    // the time cluster is rejected if it has too few combo hits
    if(_f.chunks[tc].nHits < _minCHsInCluster) {
      if(_debugLevel > 0) {
        std::cout << ">>> WARNING in TZClusterFinder::checkPopulation:" << std::endl;
        std::cout << ">>> The Time Cluster " << tc << " contains less than " << _minCHsInCluster <<
                                                            " Combo Hits and is rejected." << std::endl;
      }
      return false;
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  // recognize the pattern in the time cluster
  //-----------------------------------------------------------------------------
  bool TZClusterFinder::checkPattern(size_t tc) {

    // loop on the stations
    std::vector<size_t> pairs;
    for(size_t i=0; i<StrawId::_nstations-1; i++) {
      // take the pairs of neighboring stations with at least _minCHs combo hits per station
      if(_f.chStn[i].stnHits.size() < _minCHs || _f.chStn[i+1].stnHits.size() < _minCHs) continue;

      if(_debugLevel > 1) {
        std::cout << ">>> INFORMATION in TZClusterFinder::checkPattern:" << std::endl;
        std::cout << ">>> Candidate pair: Time Cluster " << tc << ", Stations " << i << "/" << i+1 <<
            " contain " <<_f.chStn[i].stnHits.size() << "/" << _f.chStn[i+1].stnHits.size() << " Combo Hits." << std::endl;
      }

      bool PairCheck[4] = {false};

      // check whether each station of the pair has at least one
      // rectangle of side _deltaPhi_1 that contains at least 2 CHs
      PairCheck[0] = checkStation(i);
      if(!PairCheck[0]) continue;
      PairCheck[1] = checkStation(i+1);
      if(!PairCheck[1]) continue;

      // check whether the pair of neighboring stations has at least
      // two rectangles within _deltaPhi_2 (one rectangle per station)
      PairCheck[2] = checkPair(i, i+1);
      if(!PairCheck[2]) continue;

      pairs.push_back(i);
    }

    // loop on the pairs of neighboring stations: check whether near a pair
    // there is a third station with at least a combo hit within _deltaPhi_2
    bool TSCheck = false;
    if(_thirdStation && !pairs.empty()) {
      for(size_t i=0; i<pairs.size(); i++) {
        TSCheck = checkThirdStation(pairs[i], pairs[i]+1);
        if(TSCheck) break;
      }
    }

    if(_debugLevel > 0) {
      std::cout << ">>> INFORMATION in TZClusterFinder::checkPattern:" << std::endl;
      std::cout << ">>> Found " << pairs.size() << " pairs ( ";
      for(size_t i=0; i<pairs.size(); i++) std::cout << pairs[i] << "-" << pairs[i]+1 << " ";
      std::cout << ") in Time Cluster " << tc << "." << std::endl;
    }

    if(pairs.size() < _nPairs) return false;
    if(_thirdStation && !TSCheck) return false;

    return true;
  }

  //-----------------------------------------------------------------------------
  // check the combo hits in the station
  //-----------------------------------------------------------------------------
  bool TZClusterFinder::checkStation(size_t stn) {

    bool rectangle = false;
    size_t nCHs = 0;

    // combinatorial computation: for each combo hit on the station count
    // the number of combo hits inside the cube of sides: _deltaR1 and _deltaPhi_1
    for(size_t i=0; i<_f.chStn[stn].stnHits.size(); i++) {
      for(size_t j=0; j<_f.chStn[stn].stnHits.size(); j++) {
        // same hit
        if(i==j) {
          nCHs++;
          continue;
        }
        double a[2] = {_f.chStn[stn].stnHits.at(i).x, _f.chStn[stn].stnHits.at(i).y};
        double b[2] = {_f.chStn[stn].stnHits.at(j).x, _f.chStn[stn].stnHits.at(j).y};
        double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
        double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
        double deltaR = std::abs(Rb - Ra);
        if(deltaR > _deltaR1) continue;
        // scalar product
        double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
        double deltaPhi = acos(c); // [rad]
        if(deltaPhi <= _deltaPhi_1) nCHs++;
      }
      if(_debugLevel > 1) {
        std::cout << ">>> INFORMATION in TZClusterFinder::checkStation:" << std::endl;
        std::cout << ">>> Station " << stn << ", Cubes: CentralCH(hIndex)/nCHs = " <<
                        i << "(" << _f.chStn[stn].stnHits.at(i).hIndex << ")/" << nCHs << std::endl;
      }
      if(nCHs >= _minCHs) {
        _f.chStn[stn].stnHits.at(i).center = true;
        rectangle = true;
      }
      nCHs = 0;
    }

    return rectangle;
  }

  //-----------------------------------------------------------------------------
  // check the pair of neighboring stations
  //-----------------------------------------------------------------------------
  bool TZClusterFinder::checkPair(size_t stn1, size_t stn2) {

    // combinatorial computation: check whether the pair of neighboring stations
    // has at least two cubes within _deltaR2 and _deltaPhi_2 (one cube per station)
    // note: the centers of the cubes are combo hits.
    for(size_t i=0; i<_f.chStn[stn1].stnHits.size(); i++) {
      if(_f.chStn[stn1].stnHits.at(i).center == false) continue;
      for(size_t j=0; j<_f.chStn[stn2].stnHits.size(); j++) {
        if(_f.chStn[stn2].stnHits.at(j).center == false) continue;
        double a[2] = {_f.chStn[stn1].stnHits.at(i).x, _f.chStn[stn1].stnHits.at(i).y};
        double b[2] = {_f.chStn[stn2].stnHits.at(j).x, _f.chStn[stn2].stnHits.at(j).y};
        double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
        double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
        double deltaR = std::abs(Rb - Ra);
        if(deltaR > _deltaR2) continue;
        // scalar product
        double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
        double deltaPhi = acos(c); // [rad]
        if(deltaPhi <= _deltaPhi_2) {
          if(_debugLevel > 1) {
            std::cout << ">>> INFORMATION in TZClusterFinder::checkPair:" << std::endl;
            std::cout << ">>> Found pair: Stations " << stn1 << "/" << stn2 <<
                            ", Cubes: CentralCH(hIndex)/CentralCH(hIndex) = " <<
                                  i << "(" << _f.chStn[stn1].stnHits.at(i).hIndex << ")/" <<
                                  j << "(" << _f.chStn[stn2].stnHits.at(j).hIndex << ")" << std::endl;
          }
          return true;
        }
      }
    }

    return false;
  }

  //-----------------------------------------------------------------------------
  // check the third station near the pair of neighboring stations
  //-----------------------------------------------------------------------------
  bool TZClusterFinder::checkThirdStation(size_t stn1, size_t stn2) {

    // combinatorial computation: check whether near the pair of neighboring stations
    // there is a third station with at least a combo hit within _deltaR2 and _deltaPhi_2
    // note: the centers of the cubes are combo hits.
    if(stn1 > 0) {
      for(size_t i=0; i<_f.chStn[stn1].stnHits.size(); i++) {
        if(_f.chStn[stn1].stnHits.at(i).center == false) continue;
        for(size_t j=0; j<_f.chStn[stn1-1].stnHits.size(); j++) {
          double a[2] = {_f.chStn[stn1].stnHits.at(i).x, _f.chStn[stn1].stnHits.at(i).y};
          double b[2] = {_f.chStn[stn1-1].stnHits.at(j).x, _f.chStn[stn1-1].stnHits.at(j).y};
          double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
          double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
          double deltaR = std::abs(Rb - Ra);
          if(deltaR > _deltaR2) continue;
          // scalar product
          double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
          double deltaPhi = acos(c); // [rad]
          if(deltaPhi <= _deltaPhi_2) {
            if(_debugLevel > 0) {
              std::cout << ">>> INFORMATION in TZClusterFinder::checkThirdStation:" << std::endl;
              std::cout << ">>> Found third station before the pair " << stn1 << "/" << stn2 <<
                  ", ComboHit/hIndex: " << j << "/" << _f.chStn[stn1-1].stnHits.at(j).hIndex << std::endl;
            }
            return true;
          }
        }
      }
    }

    if(stn2 < StrawId::_nstations-1) {
      for(size_t i=0; i<_f.chStn[stn2].stnHits.size(); i++) {
        if(_f.chStn[stn2].stnHits.at(i).center == false) continue;
        for(size_t j=0; j<_f.chStn[stn2+1].stnHits.size(); j++) {
          double a[2] = {_f.chStn[stn2].stnHits.at(i).x, _f.chStn[stn2].stnHits.at(i).y};
          double b[2] = {_f.chStn[stn2+1].stnHits.at(j).x, _f.chStn[stn2+1].stnHits.at(j).y};
          double Ra = std::sqrt(a[0]*a[0]+a[1]*a[1]);
          double Rb = std::sqrt(b[0]*b[0]+b[1]*b[1]);
          double deltaR = std::abs(Rb - Ra);
          if(deltaR > _deltaR2) continue;
          // scalar product
          double c = (a[0]*b[0]+a[1]*b[1])/(Ra*Rb);
          double deltaPhi = acos(c); // [rad]
          if(deltaPhi <= _deltaPhi_2) {
            if(_debugLevel > 0) {
              std::cout << ">>> INFORMATION in TZClusterFinder::checkThirdStation:" << std::endl;
              std::cout << ">>> Found third station after the pair " << stn1 << "/" << stn2 <<
                  ", ComboHit/hIndex: " << j << "/" << _f.chStn[stn2+1].stnHits.at(j).hIndex << std::endl;
            }
            return true;
          }
        }
      }
    }

    return false;
  }

  //-----------------------------------------------------------------------------
  // logic to find clusters
  //-----------------------------------------------------------------------------
  void TZClusterFinder::findClusters(TimeClusterCollection& TimeClusterColl) {

    // clear relevant data members that carry over from processing of previous event
    _f.clear_cHits();
    _f.clear_chunks();
    _f._clusterInfo._caloCluster = art::Ptr<mu2e::CaloCluster>();
    if (_diagLevel != 0 || _runDisplay != 0) { _data.clearDiagInfo(); }

    // fill cHits array to loop over
    cHitsFill();

    // loop over cHits looking for hits that can be chunked together
    chunkHits();

    // look for chunks that can be combined together
    if (_f.chunks.size() != 0) {
      _f.moreCombines = true;
      while(_f.moreCombines) {
        combineChunks();
      }
    }

    // recover hits that were missed
    recoverHits();

    // count number of protons
    countProtons(*_data._iiTC);

    // use calo clusters
    if (_useCaloClusters == 1) { checkCaloClusters(); }

    // flag bad clusters
    if (_doRefine == 1) { refineChunks(); }

    // fill TimeClusterColl and perform final fit to save
    for (size_t i=0; i<_f.chunks.size(); i++) {
      if (_f.chunks[i].nrgSelection == 0) {continue;} // don't save proton clusters
      if ((int)_f.chunks[i].nStrawHits < _clusterThresh) {continue;} // only save chunks with enough straw hits
      if (_doRefine == 1 && _f.chunks[i].goodCluster == false) {continue;} // filter out bad TCs if requested
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
      _f._clusterInfo._nsh = _f.chunks[i].nStrawHits;
      TimeClusterColl.push_back(_f._clusterInfo);
      if (_diagLevel != 0 || _runDisplay != 0) {
        _data.lineSlope.push_back(_f.chunks[i].fitter.dydx());
        _data.lineIntercept.push_back(_f.chunks[i].fitter.y0());
        _data.chi2DOF.push_back(_f.chunks[i].fitter.chi2Dof());
      }
    }

    if (_diagLevel != 0 || _runDisplay != 0) { _data._nTZClusters = (int)TimeClusterColl.size(); }
    if (_runDisplay != 0 ) { plotTZ(); }

  }

  //-----------------------------------------------------------------------------
  // print information on the TC candidates
  //-----------------------------------------------------------------------------
  void TZClusterFinder::clusterInfo(size_t tc) {

    std::cout << ">>> INFORMATION in TZClusterFinder::clusterInfo: " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << " Time Cluster " << tc                        << std::endl;
    std::cout << " StrawHits: " << _f.chunks[tc].nStrawHits    << std::endl;
    std::cout << " ComboHits: " << _f.chunks[tc].nHits         << std::endl;
    std::cout << "===========================================" << std::endl;

    if(_debugLevel < 2) return;

    for(size_t i=0; i<_f.chStn.size(); i++) {
      for(size_t j=0; j<_f.chStn[i].stnHits.size(); j++) {
        std::cout << "Station/ComboHit/hIndex/nStrawHits: " << i << "/" << j << "/" <<
                    _f.chStn[i].stnHits[j].hIndex << "/" << _f.chStn[i].stnHits[j].nStrawHits << std::endl;
        std::cout << "x/y/z/phi: " << _f.chStn[i].stnHits[j].x << "/" << _f.chStn[i].stnHits[j].y <<
                        "/" << _f.chStn[i].stnHits[j].z << "/" << _f.chStn[i].stnHits[j].phi << std::endl;
      }
    }
  }

}

using mu2e::TZClusterFinder;
DEFINE_ART_MODULE(TZClusterFinder)
