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
      fhicl::Atom<art::InputTag>     caloIntInfoTag   {Name("caloIntInfoTag"    ), Comment("IntensityInfoCalo tag"                          ) };
      fhicl::Atom<art::InputTag>     tcIntInfoTZTag   {Name("tcIntInfoTZTag"    ), Comment("IntensityInfoTimeCluster TZClusterFinder tag"   ) };
      fhicl::Atom<art::InputTag>     tcIntInfoDFTag   {Name("tcIntInfoDFTag"    ), Comment("IntensityInfoTimeCluster DeltaFinder tag"       ) };
      fhicl::Atom<art::InputTag>     trkHitIntInfoTag {Name("trkHitIntInfoTag"  ), Comment("IntensityInfoTrackerHits tag"                   ) };
      fhicl::Atom<int>               debugLevel       {Name("debugLevel"        ), Comment("turn on/off debug"                              ), 0};
    };

    //-----------------------------------------------------------------------------
    // Data Storage SetUp
    //-----------------------------------------------------------------------------
    struct Data {
      const art::Event* _event;
      struct PredictionResult {
        std::string observable; //which observable gave this prediction
        double  observableValue = 0.0; //the value of the observable
        unsigned long long NPOT = 0; //nPOT observable predicted
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
    art::InputTag   _trkHitIntInfoTag ;
    art::InputTag   _tcIntInfoTZTag ;
    art::InputTag   _tcIntInfoDFTag ;

    int              _debugLevel;
    Data             _Data;



    //-----------------------------------------------------------------------------
    // collections
    //-----------------------------------------------------------------------------
    const IntensityInfoCalo*         _caloIntInfo;
    const IntensityInfoTrackerHits*  _trkHitIntInfo;
    const IntensityInfoTimeCluster*  _tcIntInfoTZ;
    const IntensityInfoTimeCluster*  _tcIntInfoDF;


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
    void addPrediction                                   (Data& data, const std::string& obsName, const double obsVal, const unsigned long long POT);
    //Functions to reconstruct nPOT from observable
    unsigned long long POTfromCaloEnergy         (const art::Event& evt, const double caloEnergy);
    unsigned long long POTfromCaloHits           (const art::Event& evt, const int caloHits);
    unsigned long long POTfromTrackerHits        (const art::Event& evt, const int trackerHits);
    unsigned long long POTfromTZTCs              (const art::Event& evt, const int TCs);
    unsigned long long POTfromDFTCs              (const art::Event& evt, const int TCs);

  };

  //-----------------------------------------------------------------------------
  // module constructor
  //-----------------------------------------------------------------------------
  RecoNPOTMaker::RecoNPOTMaker(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}
    , _caloIntInfoTag           (config().caloIntInfoTag()                              )
    , _trkHitIntInfoTag         (config().trkHitIntInfoTag()                            )
    , _tcIntInfoTZTag           (config().tcIntInfoTZTag()                              )
    , _tcIntInfoDFTag           (config().tcIntInfoDFTag()                              )
    , _debugLevel               (config().debugLevel()                                  )
  {

    consumes<IntensityInfoCalo>(_caloIntInfoTag);
    consumes<IntensityInfoTrackerHits>(_trkHitIntInfoTag);
    consumes<IntensityInfoTimeCluster>(_tcIntInfoTZTag);
    consumes<IntensityInfoTimeCluster>(_tcIntInfoDFTag);

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
    art::Handle<IntensityInfoCalo> caloIntInfoH;
    if(evt.getByLabel(_caloIntInfoTag, caloIntInfoH)) { _caloIntInfo = caloIntInfoH.product(); }
    else                                              { _caloIntInfo = nullptr;                }

    //IntensityInfoTrackerHits
    art::Handle<IntensityInfoTrackerHits> trkHitIntInfoH;
    if(evt.getByLabel(_trkHitIntInfoTag, trkHitIntInfoH)) { _trkHitIntInfo = trkHitIntInfoH.product(); }
    else                                                  { _trkHitIntInfo = nullptr  ;                }

    //IntensityInfoTimeCluster, TZClusterFinder
    art::Handle<IntensityInfoTimeCluster> tcIntInfoTZH;
    if(evt.getByLabel(_tcIntInfoTZTag, tcIntInfoTZH)) { _tcIntInfoTZ = tcIntInfoTZH.product(); }
    else                                              { _tcIntInfoTZ = nullptr;                }

    //IntensityInfoTimeCluster, DeltaFinder
    art::Handle<IntensityInfoTimeCluster> tcIntInfoDFH;
    if(evt.getByLabel(_tcIntInfoDFTag, tcIntInfoDFH)) { _tcIntInfoDF = tcIntInfoDFH.product(); }
    else                                              { _tcIntInfoDF = nullptr;                }
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::produce(art::Event& event) {

    unsigned long long POT_caloEnergy(0), POT_caloHits(0), POT_trkHits(0), POT_TZTCs(0), POT_DFTCs(0);

    findData(event);

    //calo info
    if(_caloIntInfo) {
      const double caloEnergy = _caloIntInfo->caloEnergy();
      Data::PredictionResult caloEnergyPred;
      caloEnergyPred.observable = "caloEnergy";
      caloEnergyPred.observableValue = caloEnergy;
      POT_caloEnergy = POTfromCaloEnergy(event, caloEnergy);
      addPrediction(_Data, "caloEnergy", caloEnergy, POT_caloEnergy);

      const int nCaloHits = _caloIntInfo->nCaloHits();
      Data::PredictionResult caloHitsPred;
      caloHitsPred.observable = "caloHits";
      caloHitsPred.observableValue = nCaloHits;
      POT_caloHits = POTfromCaloHits(event, nCaloHits);
      addPrediction(_Data, "caloHits", nCaloHits, POT_caloHits);
    }
    else if(_debugLevel > 0) {std:: cout <<"[RecoNPOTMaker::produce] Did Not find IntensityInfoCalo data" << std::endl;}

    // tracker hits
    if(_trkHitIntInfo) {
      const int nHits = _trkHitIntInfo->nTrackerHits();
      Data::PredictionResult pred;
      pred.observable = "trackerHits";
      pred.observableValue = nHits;
      POT_trkHits = POTfromTrackerHits(event, nHits);
      addPrediction(_Data, "trackerHits", nHits, POT_trkHits);
    }
    else if(_debugLevel > 0) {std:: cout <<"[RecoNPOTMaker::produce] Did Not find IntensityInfoTrackerHits data" << std::endl;}

    // TZClusterFinder proton time clusters
    if(_tcIntInfoTZ) {
      const int nTCs = _tcIntInfoTZ->nProtonTCs();
      Data::PredictionResult pred;
      pred.observable = "TZTimeClusters";
      pred.observableValue = nTCs;
      POT_TZTCs = POTfromTZTCs(event, nTCs);
      addPrediction(_Data, "TZTimeClusters", nTCs, POT_TZTCs);
    }
    else if(_debugLevel > 0) {std:: cout <<"[RecoNPOTMaker::produce] Did Not find TZClusterFinder IntensityInfoTimeCluster data" << std::endl;}

    // DeltaFinder proton time clusters
    if(_tcIntInfoDF) {
      const int nTCs = _tcIntInfoDF->nProtonTCs();
      Data::PredictionResult pred;
      pred.observable = "DFTimeClusters";
      pred.observableValue = nTCs;
      POT_DFTCs = POTfromDFTCs(event, nTCs);
      addPrediction(_Data, "DFTimeClusters", nTCs, POT_DFTCs);
    }
    else if(_debugLevel > 0) {std:: cout <<"[RecoNPOTMaker::produce] Did Not find TZClusterFinder IntensityInfoTimeCluster data" << std::endl;}

    std::unique_ptr<ProtonBunchIntensity> recoPBI (new ProtonBunchIntensity);
    recoPBI->set(POT_caloEnergy); // only using calo energy to make the estimate for now


    event.put(std::move(recoPBI));

  }

  //-----------------------------------------------------------------------------
  // endJob
  //-----------------------------------------------------------------------------
  void RecoNPOTMaker::endJob(){}


  //-----------------------------------------------------------------------------
  // predict nPOT from CaloEnergy Observable
  //-----------------------------------------------------------------------------

  unsigned long long RecoNPOTMaker:: POTfromCaloEnergy(const art::Event& evt, const double caloEnergy){
    unsigned long long POT(0);

    // FIXME: Using a piece-wise fit result for now
    if (caloEnergy < 3750.) {
      POT = 700948. + 12627.9*std::pow(caloEnergy,1) + 0.293111*std::pow(caloEnergy,2);
    } else if (caloEnergy < 6300.) {
      POT = 3.25191e+06 + 11524*std::pow(caloEnergy,1) + 0.403671*std::pow(caloEnergy,2);
    } else { // linear extrapolation to high POT
      POT = -5.36566e+06 + 15314.2*std::pow(caloEnergy,1);
    }
    return POT;
  }

  unsigned long long RecoNPOTMaker:: POTfromTrackerHits(const art::Event& evt, const int /*nHits*/){
    unsigned long long POT(0);
    return POT;
  }

  unsigned long long RecoNPOTMaker:: POTfromTZTCs(const art::Event& evt, const int /*nTCs*/){
    unsigned long long POT(0);
    return POT;
  }

  unsigned long long RecoNPOTMaker:: POTfromDFTCs(const art::Event& evt, const int /*nTCs*/){
    unsigned long long POT(0);
    return POT;
  }

  unsigned long long RecoNPOTMaker:: POTfromCaloHits(const art::Event& evt, const int /*nHits*/){
    unsigned long long POT(0);
    return POT;
  }

  //-----------------------------------------------------------------------------
  //add predicted nPOT to the vector of predictions of type PredictionResult in Data Struct
  //-----------------------------------------------------------------------------

  void RecoNPOTMaker::addPrediction(Data& data, const std::string& obsName,  const double obsVal, const unsigned long long predicted) {
    Data::PredictionResult pred;
    pred.observable = obsName;
    pred.observableValue = obsVal;
    pred.NPOT = predicted;
    data.predictions.push_back(pred);
  }
}
using mu2e::RecoNPOTMaker;
DEFINE_ART_MODULE(RecoNPOTMaker)
