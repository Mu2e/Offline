//
// An EDAnalyzer module that reads back the hits created by the Calorimeter Digitization chain
//
//
// Original author

// ROOT includes
#include "TH1F.h"
#include "TF1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TFolder.h"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// Mu2e includes.
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace mu2e {
  struct DAQ_t {
    TH1F*       _hNCaloDigi;
    TH1F*       _hCDT0;
    TH1F*       _hCDROId;
    TH1F*       _hCDPeakPos;
    TH1F*       _hCDAmp;
    TH1F*       _hNSamplesPerDigi ;
    TH1F*       _hNSamplesPerEvent;

    // TH2F*       _hNSamplesVsRadius;
    // TH2F*       _hNHitsVsRadius   ;
   };

  class CaloDigiAna : public art::EDAnalyzer {

  public:

    explicit CaloDigiAna(fhicl::ParameterSet const& pset);
    virtual ~CaloDigiAna() { }

    virtual void beginRun(art::Run const& ) override;

    virtual void beginJob() override;
    virtual void endJob() override ;

    // This is called for each event.
    virtual void analyze(const art::Event& e) override;

    int waveformMaximumIndex(std::vector<int>const * waveform);

  private:

    std::string                _caloDigisModuleLabel;

    //some histograms for DAQ analyses purposes
    DAQ_t                      _histDisk[2];
    const Calorimeter*         _calorimeter; // cached pointer to the calorimeter geometry

  };



  CaloDigiAna::CaloDigiAna(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    _caloDigisModuleLabel          (pset.get<std::string>("digisTag")){}

  void CaloDigiAna::beginRun(art::Run const& ){
    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
  }
  void CaloDigiAna::beginJob(){

    art::ServiceHandle<art::TFileService> tfs;

    for (int i=0; i<2; ++i){
      art::TFileDirectory tfdir = tfs->mkdir(Form("disk%i",i));

      _histDisk[i]._hNCaloDigi         = tfdir.make<TH1F>(Form("hNHitsDisk%i", i)      ,
                                                          "Disk0: N caloDigis ", 1e4, 0, 1e4);
      _histDisk[i]._hCDT0              = tfdir.make<TH1F>(Form("hCDT0Disk%i", i)      ,
                                                          "Disk0: calo digi t0; caloDigi_t0 [ns]", 100, 0, 2e3);
      _histDisk[i]._hCDROId            = tfdir.make<TH1F>(Form("hCDROIdDisk%i", i)      ,
                                                          "Disk0: calo digi roId; caloDigi_roId", 4e3, 0, 4e3);
      _histDisk[i]._hCDPeakPos         = tfdir.make<TH1F>(Form("hCDPeakPosDisk%i", i)      ,
                                                          "Disk0: calo digi t0; caloDigi_peakPos", 100, 0, 100);
      _histDisk[i]._hCDAmp             = tfdir.make<TH1F>(Form("hCDAmpDisk%i", i)      ,
                                                          "Disk0: calo digi amplitude; caloDigi_amplitude [adc]", 250, 0, 2500);
      _histDisk[i]._hNSamplesPerDigi   = tfdir.make<TH1F>(Form("hNSampleHitsDisk%i", i),
                                                          "Disk0: N of words per caloDigi ",
                                                          200, 0., 200.);
      _histDisk[i]._hNSamplesPerEvent  = tfdir.make<TH1F>(Form("hNSampleDisk%i", i),
                                                          "Disk0: N of words per event ",
                                                          5e2, 0, 5e4);
    }
  }

  void CaloDigiAna::endJob(){}

  void CaloDigiAna::analyze(const art::Event& event) {

    art::Handle<CaloDigiCollection>   caloDigiHandle;
    event.getByLabel(_caloDigisModuleLabel, caloDigiHandle);

    const     CaloDigiCollection*       caloDigiCol(0);
    int       nCalodigi(0);
    if (caloDigiHandle.isValid()){
       caloDigiCol = caloDigiHandle.product();
       nCalodigi   = caloDigiCol->size();
    }

    //--------------------------------------------------------------------------------
    const CaloDigi   *caloDigi;

    //fill DAQ histograms
    double     amplitude(0);
    int        roId(0), crystalID(0), diskId(0), peakPos(-1);

    const std::vector<int>   *pulse;
    int                nDigi[2] = {0};
    int                nDigiWords[2] = {0};

    for (int i=0; i< nCalodigi; ++i){
      caloDigi   = &caloDigiCol->at(i);
      roId       = caloDigi->SiPMID();
      crystalID  = CaloSiPMId(roId).crystal().id();
      diskId     = _calorimeter->crystal(crystalID).diskID();
      ++nDigi[diskId];

      _histDisk[diskId]._hCDT0           ->Fill(caloDigi->t0());
      _histDisk[diskId]._hCDROId         ->Fill(roId);

      pulse      = &caloDigi->waveform();
      nDigiWords[diskId] += pulse->size();
      _histDisk[diskId]._hNSamplesPerDigi->Fill(pulse->size());

      //FIXME!
      peakPos    = waveformMaximumIndex(pulse);//caloDigi->peakPos();
      _histDisk[diskId]._hCDPeakPos      ->Fill(peakPos);
      amplitude  = pulse->at(peakPos);
      _histDisk[diskId]._hCDAmp          ->Fill(amplitude);
    }

    int nDisks(2);
    for(int i=0; i<nDisks; ++i){
      _histDisk[i]._hNCaloDigi       ->Fill(nDigi[i]);
      _histDisk[i]._hNSamplesPerEvent->Fill(nDigiWords[i]);
    }
  }

  //--------------------------------------------------------------------------------
  // temporary function used to find the location of the waveform peak in the
  // calorimeter digitized waveform
  //--------------------------------------------------------------------------------
  int CaloDigiAna::waveformMaximumIndex(std::vector<int>const * waveform) {
    int  indexMax(-1), content(0);
    for (size_t i = 0; i < waveform->size(); ++i) {
      if (waveform->at(i) > content) {
        content = waveform->at(i);
        indexMax = i;
      }
    }

    return indexMax;
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloDigiAna)
