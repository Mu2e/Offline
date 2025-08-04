///////////////////////////////////////////////////////////////////////////////
// RecoNPOTMaker
// H. Applegate
///////////////////////////////////////////////////////////////////////////////


#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

//MC DataProducts
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

//Reco DataProducts
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"

#include <cmath>
#include <iostream>

namespace mu2e {

  class RecoNPOTMaker: public art::EDProducer {

  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>               debugLevel       {Name("debugLevel"       ), Comment("turn on/off debug"           ) };
      // fhicl::Atom<art::InputTag>     chCollTag        {Name("chCollTag"         ), Comment("ComboHitCollection"         ) };
      // fhicl::Atom<art::InputTag>     tcCollTag        {Name("tcCollTag"         ), Comment("TimeClusterCollection"      ) };
      // fhicl::Atom<art::InputTag>     hsCollTag        {Name("hsCollTag"         ), Comment("HelixSeedCollection"        ) };
      fhicl::Atom<art::InputTag>     caloIntInfoTag   {Name("caloIntInfoTag"    ), Comment("IntensityInfoCalo"          ) };
      // fhicl::Atom<art::InputTag>     tcIntInfoTZTag   {Name("tcIntInfoTZTag"    ), Comment("IntensityInfoTimeCluster"   ) };
      // fhicl::Atom<art::InputTag>     tcIntInfoFlagTag {Name("tcIntInfoFlagTag"  ), Comment("IntensityInfoTimeCluster"   ) };
      // fhicl::Atom<art::InputTag>     trkHitIntInfoTag {Name("trkHitIntInfoTag"  ), Comment("IntensityInfoTrackerHits"   ) };
      fhicl::Atom<art::InputTag>     evtWeightTag     {Name("evtWeightTag"      ), Comment("ProtonBunchIntensity"       ) };

    };

    //-----------------------------------------------------------------------------
    // Data Storage SetUp
    //-----------------------------------------------------------------------------
    struct Data {
      const art::Event* _event;
      unsigned long long actualNPOT = 0; //nPOT from ProtonBunchIntensity MC

      struct PredictionResult {
        std::string observable; //which observable gave this prediction
        double  observableValue = 0.0; //the value of the observable
        unsigned long long predictedNPOT = 0; //nPOT observable predicted
      };

      std::vector<PredictionResult> predictions;

    };





  private:

    //-----------------------------------------------------------------------------
    // Data Members
    //-----------------------------------------------------------------------------
    int              _debugLevel;

    //-----------------------------------------------------------------------------
    // Event Object Labels
    //-----------------------------------------------------------------------------
    // art::InputTag   _chCollTag ;
    // art::InputTag   _tcCollTag ;
    // art::InputTag   _hsCollTag ;
    art::InputTag   _caloIntInfoTag ;
    //  art::InputTag   _tcIntInfoTZTag ;
    // art::InputTag   _tcIntInfoFlagTag ;
    // art::InputTag   _trkHitIntInfoTag ;
    art::InputTag   _evtWeightTag ;

    const art::Event*                  _event;
    Data                               _Data;



    //-----------------------------------------------------------------------------
    // collections
    //-----------------------------------------------------------------------------
    const ProtonBunchIntensity*      _evtWeight;
    const IntensityInfoCalo*         _caloIntInfo;
    //const IntensityInfoTimeCluster*  _tcIntInfoTZ;
    //const IntensityInfoTimeCluster*  _tcIntInfoFlag;
    //const IntensityInfoTrackerHits*  _trkHitIntInfo;
    //const TimeClusterCollection*     _tcColl;
    //const HelixSeedCollection*       _hsColl;
    //const ComboHitCollection*        _chColl;


    ProtonBunchIntensity*            recoPBI; //Pointer to Collection Module Produces

    //-----------------------------------------------------------------------------
    // functions
    //-----------------------------------------------------------------------------

  public:

    explicit RecoNPOTMaker(const art::EDProducer::Table<Config>& config);
    virtual ~RecoNPOTMaker();

    virtual void beginJob ();
    virtual void beginRun (art::Run&);
    virtual void produce  (art::Event& e);
    virtual void endJob   ();

    //-----------------------------------------------------------------------------
    // helper functions
    //-----------------------------------------------------------------------------
    bool findData                                        (const art::Event& evt);


    //Add Prediction to Data Storage
    void addPrediction                                   (Data& data, const std::string& obsName, double& obsVal, unsigned long long& predicted);
    //Functions to Reconstruct nPOT from Observable
    unsigned long long predictnPOTfromCaloEnergy         (const art::Event& evt, double& caloEnergy);

  };


  //-----------------------------------------------------------------------------
  // module constructor
  //-----------------------------------------------------------------------------
  RecoNPOTMaker::RecoNPOTMaker(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    _debugLevel               (config().debugLevel()                                  ),
    // _chCollTag                (config().chCollTag()                                   ),
    // _tcCollTag                (config().tcCollTag()                                   ),
    // _hsCollTag                (config().hsCollTag()                                   ),
    _caloIntInfoTag           (config().caloIntInfoTag()                              ),
    // _tcIntInfoTZTag           (config().tcIntInfoTZTag()                              ),
    // _tcIntInfoFlagTag         (config().tcIntInfoFlagTag()                            ),
    // _trkHitIntInfoTag         (config().trkHitIntInfoTag()                            ),
    _evtWeightTag             (config().evtWeightTag()                                )
  {

    consumes<ProtonBunchIntensity>(_evtWeightTag);
    // consumes<ComboHitCollection>(_chCollTag);
    // consumes<TimeClusterCollection>(_tcCollTag);
    // consumes<HelixSeedCollection>(_hsCollTag);
    // consumes<IntensityInfoTrackerHits>(_trkHitIntInfoTag);
    consumes<IntensityInfoCalo>(_caloIntInfoTag);
    // consumes<IntensityInfoTimeCluster>(_tcIntInfoTZTag);
    // consumes<IntensityInfoTimeCluster>(_tcIntInfoFlagTag);

    produces<ProtonBunchIntensity>();

  }

  //-----------------------------------------------------------------------------
  // destructor
  //-----------------------------------------------------------------------------
  RecoNPOTMaker::~RecoNPOTMaker() {}

  //-----------------------------------------------------------------------------
  // beginJob
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::beginJob(){}

  //-----------------------------------------------------------------------------
  // beginRun
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::beginRun(art::Run& ) {}

  //-----------------------------------------------------------------------------
  // find input things
  //-----------------------------------------------------------------------------
  bool RecoNPOTMaker::findData(const art::Event& evt) {

    //ProtonBunchIntensity
    auto evtWeightH = evt.getValidHandle<ProtonBunchIntensity>(_evtWeightTag);
    if (evtWeightH.product() != 0){_evtWeight = evtWeightH.product(); }
    else {
      _evtWeight = 0;
      std::cout << ">>> ERROR in RecoNPOTMaker::findData: ProtonBunchIntensity not found." << std::endl;
      return false;
    }

    //IntensityInfoCalo
    auto caloIntInfoH = evt.getValidHandle<IntensityInfoCalo>(_caloIntInfoTag);
    if (caloIntInfoH.product() != 0){_caloIntInfo = caloIntInfoH.product(); }
    else {
      _caloIntInfo = 0;
      std::cout << ">>> ERROR in RecoNPOTMaker::findData: IntensityInfoCalo not found." << std::endl;
      return false;
    }

    /* //IntensityInfoTimeCluster, TTTZClusterFinder
       auto tcIntInfoTZH = evt.getValidHandle<IntensityInfoTimeCluster>(_tcIntInfoTZTag);
       if (tcIntInfoTZH.product() != 0){_tcIntInfoTZ = tcIntInfoTZH.product(); }
       else {
       _tcIntInfoTZ = 0;
       std::cout << ">>> ERROR in RecoNPOTMaker::findData: IntensityInfoTimeCluster with TTTZClusterFinder not found." << std::endl;
       return false;
       }

       //IntensityInfoTimeCluster, TTflagPH
       auto tcIntInfoFlagH = evt.getValidHandle<IntensityInfoTimeCluster>(_tcIntInfoFlagTag);
       if (tcIntInfoFlagH.product() != 0){_tcIntInfoFlag = tcIntInfoFlagH.product(); }
       else {
       _tcIntInfoFlag = 0;
       std::cout << ">>> ERROR in RecoNPOTMaker::findData: IntensityInfoTimeCluster with TTflagPH not found." << std::endl;
       return false;
       }

       //IntensityInfoTrackerHits
       auto trkHitIntInfoH = evt.getValidHandle<IntensityInfoTrackerHits>(_trkHitIntInfoTag);
       if (trkHitIntInfoH.product() != 0){_trkHitIntInfo = trkHitIntInfoH.product(); }
       else {
       _trkHitIntInfo = 0;
       std::cout << ">>> ERROR in RecoNPOTMaker::findData: IntensityInfoTrackerHits not found." << std::endl;
       return false;
       }

       //ComboHitCollection
       auto chCollH = evt.getValidHandle<ComboHitCollection>(_chCollTag);
       if (chCollH.product() != 0){ _chColl = chCollH.product(); }
       else{
       _chColl = 0;
       std::cout << ">>> ERROR in RecoNPOTMaker::findData: ComboHitCollection not found." << std::endl;
       return false;
       }

       //TimeClusterCollection
       auto tcCollH = evt.getValidHandle<TimeClusterCollection>(_tcCollTag);
       if (tcCollH.product() != 0){ _tcColl = tcCollH.product(); }
       else{
       _tcColl = 0;
       std::cout << ">>> ERROR in RecoNPOTMaker::findData: TimeClusterCollection not found." << std::endl;
       return false;
       }

       //HelixSeedCollection
       auto hsCollH = evt.getValidHandle<HelixSeedCollection>(_hsCollTag);
       if (hsCollH.product() != 0) {_hsColl = hsCollH.product(); }
    */
    return true;
  }
  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::produce(art::Event& event) {

    unsigned long long caloPredicted(0.0);

    bool dataexist = findData(event);
    if (dataexist) {
      //ProtonBunchIntensity
      if(_evtWeight != nullptr){
        unsigned long long npot = _evtWeight->intensity();
        _Data.actualNPOT = npot;
        //std::cout << "[RecoNPOTMaker::produce] >> nPOT = " << npot << std::endl;
      }
      else {std::cout << "[RecoNPOTMaker::produce] Did Not find Proton Bunch Intensity data" << std::endl;}

      //IntensityInfoCalo
      if(_caloIntInfo != nullptr){
        double caloEnergy = _caloIntInfo->caloEnergy();
        // int nCalo = _caloIntInfo->nCaloHits();
        // int nCaphri = _caloIntInfo->nCaphriHits();


        Data::PredictionResult caloPred;
        caloPred.observable = "caloEnergy";
        caloPred.observableValue = caloEnergy;
        caloPredicted = predictnPOTfromCaloEnergy(event, caloEnergy);
        addPrediction(_Data, "caloEnergy", caloEnergy, caloPredicted);
      }
      else {std:: cout <<"[RecoNPOTMaker::produce] Did Not find IntensityInfoCalo data" << std::endl;}

      /*
      //IntensityInfoTimeCluster, TTTZClusterFinder
      if(_tcIntInfoTZ != nullptr){std::cout << "[RecoNPOTMaker::produce] >> nProtonTCs = " <<  _tcIntInfoTZ->nProtonTCs() << std::endl; }
      else {std:: cout <<"[RecoNPOTMaker::produce] Did Not find IntensityInfoTimeCluster with TTTZClusterFinder data" << std::endl;}

      //IntensityInfoTimeCluster,TTflagPH
      if(_tcIntInfoFlag != nullptr){std::cout << "[RecoNPOTMaker::produce] >> nProtonTCs = " <<  _tcIntInfoFlag->nProtonTCs() << std::endl; }
      else {std:: cout <<"[RecoNPOTMaker::produce] Did Not find IntensityInfoTimeCluster with TTflagPH data" << std::endl;}


      //IntensityInfoTrackerHits
      if(_trkHitIntInfo != nullptr){std::cout << "[RecoNPOTMaker::produce] >> nTrackerHits = " <<  _trkHitIntInfo->nTrackerHits() << std::endl; }
      else {std:: cout <<"[RecoNPOTMaker::produce] Did Not find IntensityInfoTrackerHits data" << std::endl;}


      //ComboHitCollection
      if (_chColl != nullptr){ std::cout << "[RecoNPOTMaker::produce] >> total nStrawHits = " << _chColl->nStrawHits() << std::endl; }
      else {std::cout << "[RecoNPOTMaker::produce] Did Not find ComboHitCollection data" << std::endl;}


      //TimeClusterCollection
      if(_tcColl != nullptr) {
      size_t totalHits = 0;
      for (const auto& cluster : *_tcColl) { totalHits += cluster.nhits(); }
      std::cout << "[RecoNPOTMaker::produce] >> Total nhits = " << totalHits << std::endl;
      }
      else {std::cout << "[RecoNPOTMaker::produce] Did Not find TimeClusterCollection data" << std::endl;}



      //HelixSeedCollection
      if(_hsColl != nullptr){
      size_t i = 0;
      for (const auto& helixSeed : *_hsColl) {
      double radius = helixSeed.helix().radius();
      std::cout << "[RecoNPOTMaker::produce] >> Helix[" << i++ << "] Radius = " << radius << std::endl;
      }
      }
      else {std::cout << "[RecoNPOTMaker::produce] Did Not find HelixSeedCollection data" << std::endl;}
      */
    }

    else {std::cout << "[RecoNPOTMaker::produce] Data does not Exist" << std::endl; }



    std::unique_ptr<ProtonBunchIntensity> recoPBI (new ProtonBunchIntensity);
    recoPBI->set(caloPredicted);


    event.put(std::move(recoPBI));

  }

  //-----------------------------------------------------------------------------
  // endJob
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::endJob(){}


  //------------------------ My Helper Functions---------------------------------

  //-----------------------------------------------------------------------------
  // predict nPOT from CaloEnergy Observable
  //-----------------------------------------------------------------------------

  unsigned long long RecoNPOTMaker:: predictnPOTfromCaloEnergy(const art::Event& evt, double& caloEnergy){
    unsigned long long predictednPOT;

    if (caloEnergy >= 0.000000 && caloEnergy <= 3750.000000) {
      predictednPOT = 700948 + 12627.9*pow(caloEnergy,1) + 0.293111*pow(caloEnergy,2);
    }
    else if (caloEnergy >= 3750.000000 && caloEnergy <= 6300.000000) {
      predictednPOT = 3.25191e+06 + 11524*pow(caloEnergy,1) + 0.403671*pow(caloEnergy,2);
    }
    else {
      predictednPOT = -5.36566e+06 + 15314.2*pow(caloEnergy,1);
    }

    return predictednPOT;
  }
  //-----------------------------------------------------------------------------
  //add predicted nPOT to my vector of predictions of type PredictionResult in my Data Struct
  //-----------------------------------------------------------------------------

  void RecoNPOTMaker::addPrediction(Data& data, const std::string& obsName,  double& obsVal, unsigned long long& predicted) {
    Data::PredictionResult pred;
    pred.observable = obsName;
    pred.observableValue = obsVal;
    pred.predictedNPOT = predicted;
    data.predictions.push_back(pred);
  }

}
using mu2e::RecoNPOTMaker;
DEFINE_ART_MODULE(RecoNPOTMaker)
