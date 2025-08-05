///////////////////////////////////////////////////////////////////////////////
// RecoNPOTMaker
// H. Applegate, 2025
///////////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"

//Reco DataProducts
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
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
      fhicl::Atom<art::InputTag>     caloIntInfoTag   {Name("caloIntInfoTag"    ), Comment("IntensityInfoCalo"          ) };
      // fhicl::Atom<art::InputTag>     tcIntInfoTZTag   {Name("tcIntInfoTZTag"    ), Comment("IntensityInfoTimeCluster"   ) };
      // fhicl::Atom<art::InputTag>     tcIntInfoFlagTag {Name("tcIntInfoFlagTag"  ), Comment("IntensityInfoTimeCluster"   ) };
      // fhicl::Atom<art::InputTag>     trkHitIntInfoTag {Name("trkHitIntInfoTag"  ), Comment("IntensityInfoTrackerHits"   ) };
      fhicl::Atom<int>               debugLevel       {Name("debugLevel"       ), Comment("turn on/off debug"           ), 0};

    };

    //-----------------------------------------------------------------------------
    // Data Storage SetUp
    //-----------------------------------------------------------------------------
    struct Data {
      const art::Event* _event;
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

    //-----------------------------------------------------------------------------
    // Event Object Labels
    //-----------------------------------------------------------------------------
    art::InputTag   _caloIntInfoTag ;
    // art::InputTag   _trkHitIntInfoTag ;
    // art::InputTag   _tcIntInfoTZTag ;
    // art::InputTag   _tcIntInfoFlagTag ;

    int              _debugLevel;
    Data             _Data;



    //-----------------------------------------------------------------------------
    // collections
    //-----------------------------------------------------------------------------
    const IntensityInfoCalo*         _caloIntInfo;
    //const IntensityInfoTrackerHits*  _trkHitIntInfo;
    //const IntensityInfoTimeCluster*  _tcIntInfoTZ;
    //const IntensityInfoTimeCluster*  _tcIntInfoFlag;


    ProtonBunchIntensity*            recoPBI; //Pointer to collection module produces

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
    void findData                                        (const art::Event& evt);
    //Add prediction to data storage
    void addPrediction                                   (Data& data, const std::string& obsName, const double obsVal, const unsigned long long predicted);
    //Functions to reconstruct nPOT from observable
    unsigned long long predictnPOTfromCaloEnergy         (const art::Event& evt, const double caloEnergy);

  };

  //-----------------------------------------------------------------------------
  // module constructor
  //-----------------------------------------------------------------------------
  RecoNPOTMaker::RecoNPOTMaker(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}
    , _caloIntInfoTag           (config().caloIntInfoTag()                              )
    // , _trkHitIntInfoTag         (config().trkHitIntInfoTag()                            )
    // , _tcIntInfoTZTag           (config().tcIntInfoTZTag()                              )
    // , _tcIntInfoFlagTag         (config().tcIntInfoFlagTag()                            )
    , _debugLevel               (config().debugLevel()                                  )
  {

    consumes<IntensityInfoCalo>(_caloIntInfoTag);
    // consumes<IntensityInfoTrackerHits>(_trkHitIntInfoTag);
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
  void RecoNPOTMaker::findData(const art::Event& evt) {

    //IntensityInfoCalo
    auto caloIntInfoH = evt.getValidHandle<IntensityInfoCalo>(_caloIntInfoTag);
    _caloIntInfo = caloIntInfoH.product();

    /*
       //IntensityInfoTrackerHits
       auto trkHitIntInfoH = evt.getValidHandle<IntensityInfoTrackerHits>(_trkHitIntInfoTag);
       _trkHitIntInfo = trkHitIntInfoH.product();

       //IntensityInfoTimeCluster, TTTZClusterFinder
       auto tcIntInfoTZH = evt.getValidHandle<IntensityInfoTimeCluster>(_tcIntInfoTZTag);
       _tcIntInfoTZ = tcIntInfoTZH.product();

       //IntensityInfoTimeCluster, TTflagPH
       auto tcIntInfoFlagH = evt.getValidHandle<IntensityInfoTimeCluster>(_tcIntInfoFlagTag);
       _tcIntInfoFlag = tcIntInfoFlagH.product();
    */
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::produce(art::Event& event) {

    unsigned long long caloPredicted(0);

    findData(event);

    //IntensityInfoCalo
    if(_caloIntInfo) {
      const double caloEnergy = _caloIntInfo->caloEnergy();
      // int nCalo = _caloIntInfo->nCaloHits();
      // int nCaphri = _caloIntInfo->nCaphriHits();

      Data::PredictionResult caloPred;
      caloPred.observable = "caloEnergy";
      caloPred.observableValue = caloEnergy;
      caloPredicted = predictnPOTfromCaloEnergy(event, caloEnergy);
      addPrediction(_Data, "caloEnergy", caloEnergy, caloPredicted);
    }
    else {std:: cout <<"[RecoNPOTMaker::produce] Did Not find IntensityInfoCalo data" << std::endl;}

    std::unique_ptr<ProtonBunchIntensity> recoPBI (new ProtonBunchIntensity);
    recoPBI->set(caloPredicted);


    event.put(std::move(recoPBI));

  }

  //-----------------------------------------------------------------------------
  // endJob
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::endJob(){}


  //-----------------------------------------------------------------------------
  // predict nPOT from CaloEnergy Observable
  //-----------------------------------------------------------------------------

  unsigned long long RecoNPOTMaker:: predictnPOTfromCaloEnergy(const art::Event& evt, const double caloEnergy){
    unsigned long long predictednPOT(0);

    // FIXME: Using a piece-wise fit result for now
    if (caloEnergy < 3750.) {
      predictednPOT = 700948. + 12627.9*std::pow(caloEnergy,1) + 0.293111*std::pow(caloEnergy,2);
    } else if (caloEnergy < 6300.) {
      predictednPOT = 3.25191e+06 + 11524*std::pow(caloEnergy,1) + 0.403671*std::pow(caloEnergy,2);
    } else { // linear extrapolation to high POT
      predictednPOT = -5.36566e+06 + 15314.2*std::pow(caloEnergy,1);
    }

    return predictednPOT;
  }

  //-----------------------------------------------------------------------------
  //add predicted nPOT to the vector of predictions of type PredictionResult in Data Struct
  //-----------------------------------------------------------------------------

  void RecoNPOTMaker::addPrediction(Data& data, const std::string& obsName,  const double obsVal, const unsigned long long predicted) {
    Data::PredictionResult pred;
    pred.observable = obsName;
    pred.observableValue = obsVal;
    pred.predictedNPOT = predicted;
    data.predictions.push_back(pred);
  }
}
using mu2e::RecoNPOTMaker;
DEFINE_ART_MODULE(RecoNPOTMaker)
