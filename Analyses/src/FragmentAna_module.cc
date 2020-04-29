// ======================================================================
//
// StrawAndCaloDigisFromFragments_plugin:  Add tracker/cal data products to the event
//
// ======================================================================


// ROOT includes
#include "TH1F.h"
//#include "TFolder.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Principal/Handle.h"
#include "mu2e-artdaq-core/Overlays/FragmentType.hh"
#include "mu2e-artdaq-core/Overlays/ArtFragmentReader.hh"

#include <artdaq-core/Data/Fragment.hh>

// Mu2e includes
#include "DataProducts/inc/TrkTypes.hh"
#include "DataProducts/inc/StrawId.hh"

#include <iostream>
#include <string>
#include <memory>

// ======================================================================
namespace mu2e {
  
  class FragmentAna : public art::EDAnalyzer
  {

  public:
    struct Config 
    {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>            diagLevel             { Name("diagLevel"),         Comment("diagnostic level")};
      fhicl::Atom<int>            parseCAL              { Name("parseCAL"),          Comment("parseCAL")};
      fhicl::Atom<int>            parseTRK              { Name("parseTRK"),          Comment("parseTRK")};
      fhicl::Atom<art::InputTag>  caloTag               { Name("caloTag"),           Comment("caloTag") };
      fhicl::Atom<art::InputTag>  trkTag                { Name("trkTag"),            Comment("trkTag") };
    };
    using EventNumber_t = art::EventNumber_t;
    using adc_t = mu2e::ArtFragmentReader::adc_t;
  
    explicit  FragmentAna(const art::EDAnalyzer::Table<Config>& config);
    virtual  ~FragmentAna()  { }

    virtual void beginJob() override;
    virtual void endJob() override ;

    virtual void analyze(const art::Event& e) override;

  private:
    int   diagLevel_;

    int   parseCAL_;
    int   parseTRK_;

    art::InputTag trkFragmentsTag_;
    art::InputTag caloFragmentsTag_;

    TH1F*  _hTrkNFragment;
    TH1F*  _hTrkStrawId; 
    TH1F*  _hTrkTDC[4];	 
    TH1F*  _hTrkTOT;	 
    TH1F*  _hTrkMeanADC;
    TH1F*  _hTrkMaxADC;
    TH1F*  _hTrkWfSize;	 
	                 
    TH1F*  _hCalNFragment;
    TH1F*  _hCalROId;	 
    TH1F*  _hCalT0;	 
    TH1F*  _hCalPeakPos; 
    TH1F*  _hCalWfSize;  
  
  };

  // ======================================================================

  FragmentAna::FragmentAna(const art::EDAnalyzer::Table<Config>& config):
    art::EDAnalyzer{ config },
    diagLevel_       (config().diagLevel()),
    parseCAL_        (config().parseCAL()),
    parseTRK_        (config().parseTRK()),
    trkFragmentsTag_ (config().trkTag()),
    caloFragmentsTag_(config().caloTag()){}

  //--------------------------------------------------------------------------------
  // create the histograms
  //--------------------------------------------------------------------------------
  void FragmentAna::
  beginJob(){
    art::ServiceHandle<art::TFileService> tfs; 

    art::TFileDirectory calDir = tfs->mkdir("calorimeter");
    art::TFileDirectory trkDir = tfs->mkdir("tracker");

    _hTrkNFragment = trkDir.make<TH1F>("hTrkNFragment", "n fragments from the tracker; nTrkFragments", 1000, 0., 10000.);
    _hTrkStrawId   = trkDir.make<TH1F>("hTrkStrawId"  , "trk fragment strawId; strawId", 20000, 0., 20000.); 
    _hTrkTDC[0]    = trkDir.make<TH1F>("hTrkTDC0"     , "trk fragment TDC0; TDC[0]", 264,   0.,  264000);
    _hTrkTDC[1]    = trkDir.make<TH1F>("hTrkTDC1"     , "trk fragment TDC1; TDC[1]", 264,   0.,  264000); 
    _hTrkTDC[2]    = trkDir.make<TH1F>("hTrkTDCMean"  , "trk fragment average TDC; (TDC[0]+TDC[1])/2", 2000, 0., 20000.);	     
    _hTrkTDC[3]    = trkDir.make<TH1F>("hTrkTDCDelta" , "trk fragment delta TDC; TDC[1]-TDC[0]", 220, -100., 10000.);
    _hTrkTOT       = trkDir.make<TH1F>("hTrkTOT"      , "trk fragment average TOT; (TOT[0]+TOT[1])/2", 100, 0., 200.);	 
    _hTrkMeanADC   = trkDir.make<TH1F>("hTrkMeanADC"  , "trk fragment Mean ADC; <ADC>", 250,  0., 2500. );
    _hTrkMaxADC    = trkDir.make<TH1F>("hTrkMaxADC"   , "trk fragment Max ADC; Max_ADC",  250,  0., 2500. );
    _hTrkWfSize    = trkDir.make<TH1F>("hTrkWfSize"   , "trk fragment waveform size; trkFragment_wf_size", 20, 0., 20.);  

    _hCalNFragment = calDir.make<TH1F>("hCalNFragment", "n fragments from the calorimeter; nCalFragments", 400, 0., 400.);
    _hCalROId 	   = calDir.make<TH1F>("hCalROId", "calo fragment roId; calFragment_roId", 4000, 0., 4000.);
    _hCalT0	   = calDir.make<TH1F>("hCalT0", "calo fragment t0; calFragment_t0 [ns]", 200, 0., 2000.);
    _hCalPeakPos   = calDir.make<TH1F>("hCalPeakPos", "calo fragment peakPos; calFragment_peakPos", 100, 0., 100.);
    _hCalWfSize    = calDir.make<TH1F>("hCalWfSize", "calo fragment waveform size; calFragment_wf_size", 100, 0., 100.);
  }

  void FragmentAna::
  endJob(){}

  //--------------------------------------------------------------------------------
  void
  FragmentAna::
  analyze( const art::Event & event )
  {
    art::EventNumber_t eventNumber = event.event();

    art::Handle<artdaq::Fragments> trkFragments, calFragments;
    size_t numTrkFrags(0), numCalFrags(0);
    if (parseTRK_){
      event.getByLabel(trkFragmentsTag_ , trkFragments);
      if (!trkFragments.isValid()){            
	return;
      }
      numTrkFrags = trkFragments->size();
      _hTrkNFragment->Fill(numTrkFrags);
    }
    if (parseCAL_){
      event.getByLabel(caloFragmentsTag_, calFragments);
      if (!calFragments.isValid()){
	return;
      }
      numCalFrags = calFragments->size();
      _hCalNFragment->Fill(numCalFrags);
    }

    if( diagLevel_ > 1 ) {
      std::cout << std::dec << "Producer: Run " << event.run() << ", subrun " << event.subRun()
		<< ", event " << eventNumber << " has " << std::endl;
      std::cout << numTrkFrags << " TRK fragments, and ";
      std::cout << numCalFrags << " CAL fragments." << std::endl;

      size_t totalSize = 0;
      for(size_t idx = 0; idx < numTrkFrags; ++idx) {
	auto size = ((*trkFragments)[idx]).size();// * sizeof(artdaq::RawDataType);
	totalSize += size;
	//      std::cout << "\tTRK Fragment " << idx << " has size " << size << std::endl;
      }
      for(size_t idx = 0; idx < numCalFrags; ++idx) {
	auto size = ((*calFragments)[idx]).size();// * sizeof(artdaq::RawDataType);
	totalSize += size;
	//      std::cout << "\tCAL Fragment " << idx << " has size " << size << std::endl;
      }
      totalSize *= sizeof(artdaq::RawDataType);

      std::cout << "\tTotal Size: " << (int)totalSize << " bytes." << std::endl;  
    }
    std::string curMode = "TRK";

    // Loop over the TRK and CAL fragments
    for (size_t idx = 0; idx < numTrkFrags+numCalFrags; ++idx) {

      auto curHandle = trkFragments;
      size_t curIdx = idx;
      if(idx>=numTrkFrags) {
	curIdx = idx-numTrkFrags;
	curHandle = calFragments;
      }
      const auto& fragment((*curHandle)[curIdx]);
    
      mu2e::ArtFragmentReader cc(fragment);
    
      if( diagLevel_ > 1 ) {
	std::cout << std::endl;
	std::cout << "ArtFragmentReader: ";
	std::cout << "\tBlock Count: " << std::dec << cc.block_count() << std::endl;
	std::cout << "\tByte Count: " << cc.byte_count() << std::endl;
	std::cout << std::endl;
	std::cout << "\t" << "====== Example Block Sizes ======" << std::endl;
	for(size_t i=0; i<10; i++) {
	  if(i <cc.block_count()) {
	    std::cout << "\t" << i << "\t" << cc.blockIndexBytes(i) << "\t" << cc.blockSizeBytes(i) << std::endl;
	  }
	}
	std::cout << "\t" << "=========================" << std::endl;
      }
    
      mu2e::FragmentType mode_;

      for(size_t curBlockIdx=0; curBlockIdx<cc.block_count(); curBlockIdx++) {

	size_t blockStartBytes = cc.blockIndexBytes(curBlockIdx);
	size_t blockEndBytes = cc.blockEndBytes(curBlockIdx);

	if( diagLevel_ > 1 ) {
	  std::cout << "BLOCKSTARTEND: " << blockStartBytes << " " << blockEndBytes << " " << cc.blockSizeBytes(curBlockIdx)<< std::endl;
	  std::cout << "IndexComparison: " << cc.blockIndexBytes(0)+16*(0+3*curBlockIdx) << "\t";
	  std::cout                        << cc.blockIndexBytes(curBlockIdx)+16*(0+3*0) << std::endl;
	}

	auto hdr = cc.GetHeader(curBlockIdx);
	if(hdr == nullptr) {
	  mf::LogError("FragmentAna") << "Unable to retrieve header from block " << curBlockIdx << "!" << std::endl;
	  continue;
	}

	eventNumber = hdr->GetTimestamp();
      
	if(idx < numTrkFrags){
	  mode_ = mu2e::FragmentType::TRK;//"TRK";
	}else {
	  mode_ = mu2e::FragmentType::CAL;//"CAL";
	}

	// Parse phyiscs information from TRK packets
	if(mode_ == mu2e::FragmentType::TRK && hdr->PacketCount>0 && parseTRK_>0) {

	  // Create the StrawDigi data products
	  auto trkData = cc.GetTrackerData(curBlockIdx);
	  if(trkData == nullptr) {
	    mf::LogError("FragmentAna") << "Error retrieving Tracker data from DataBlock " << curBlockIdx << "! Aborting processing of this block!";
	    continue;
	  }

	  mu2e::StrawId sid(trkData->StrawIndex);
	  mu2e::TrkTypes::TDCValues   tdc = {trkData->TDC0 , trkData->TDC1};
	  mu2e::TrkTypes::TOTValues   tot = {trkData->TOT0 , trkData->TOT1};
	  mu2e::TrkTypes::ADCWaveform adcs  = trkData->Waveform();	
	  int sum{0};
	  unsigned short maxadc{0};
	  for ( auto adc : adcs ){
	    sum += adc;
	    maxadc = std::max( maxadc, adc);
	  }
	  // Fill the StrawDigiCollection
	  _hTrkStrawId->Fill(sid.asUint16());
	  _hTrkTDC[0] ->Fill(tdc[0]);
	  _hTrkTDC[1] ->Fill(tdc[1]);
	  _hTrkTDC[2] ->Fill((tdc[0] + tdc[1])/2.);
	  _hTrkTDC[3] ->Fill(tdc[1] - tdc[0]);
	  _hTrkTOT    ->Fill((tot[0] + tot[1])/2.);
	  int mean = ( adcs.size() != 0 ) ? sum/adcs.size() : -1.;
	  _hTrkMeanADC->Fill(mean);
	  _hTrkMaxADC ->Fill(maxadc);
	  _hTrkWfSize ->Fill(adcs.size());


	} else if(mode_ == mu2e::FragmentType::CAL && hdr->PacketCount>0 && parseCAL_>0) {	// Parse phyiscs information from CAL packets
	
	  auto calData = cc.GetCalorimeterData(curBlockIdx);
	  if(calData == nullptr) {
	    mf::LogError("FragmentAna") << "Error retrieving Calorimeter data from block " << curBlockIdx << "! Aborting processing of this block!";
	    continue;
	  }

	  if( diagLevel_ > 0 ) {
	    std::cout <<"[StrawAndCaloDigiFromFragments] NEW CALDATA: NumberOfHits "<< calData->NumberOfHits << std::endl;
	  }

	  bool err = false;
	  for(size_t hitIdx = 0; hitIdx<calData->NumberOfHits; hitIdx++) {

	    // Fill the CaloDigiCollection
	    const mu2e::ArtFragmentReader::CalorimeterHitReadoutPacket* hitPkt(0);
	    hitPkt = cc.GetCalorimeterReadoutPacket(curBlockIdx, hitIdx);
	    if(hitPkt == nullptr) {
	      mf::LogError("FragmentAna") << "Error retrieving Calorimeter data from block " << curBlockIdx << " for hit " << hitIdx << "! Aborting processing of this block!";
	      err = true;
	      break;
	    }
	  
	    auto first = cc.GetCalorimeterReadoutSample(curBlockIdx,hitIdx,0);
	    auto last  = cc.GetCalorimeterReadoutSample(curBlockIdx, hitIdx, hitPkt->NumberOfSamples - 1);
	    if(first == nullptr || last == nullptr) {
	      mf::LogError("FragmentAna") << "Error retrieving Calorimeter samples from block " << curBlockIdx << " for hit " << hitIdx << "! Aborting processing of this block!";
	      err = true;
	      break;
	    }

	    //the second argument is not included in the vector, so we need to add "+1"
	    // because we want the "last" item included
	    std::vector<int> cwf(first,last+1);

	    // IMPORTANT NOTE: we don't have a final
	    // mapping yet so for the moment, the BoardID field (described in docdb 4914) is just a
	    // placeholder. Because we still need to know which crystal a hit belongs to, we are
	    // temporarily storing the 4-bit apdID and 12-bit crystalID in the Reserved DIRAC A slot.
	    // Also, note that until we have an actual map, channel index does not actually correspond
	    // to the physical readout channel on a ROC.
	    adc_t crystalID  = hitPkt->DIRACB & 0x0FFF;
	    adc_t roId      = hitPkt->DIRACB >> 12;
	    _hCalROId   ->Fill(crystalID*2 + roId);
	    _hCalT0     ->Fill(hitPkt->Time);
	    _hCalPeakPos->Fill(hitPkt->IndexOfMaxDigitizerSample);
	    _hCalWfSize ->Fill(cwf.size());

	  } // End loop over readout channels in DataBlock
	  if(err) continue;
	
	} // End Cal Mode
      
      } // End loop over DataBlocks within fragment 
      
    } // Close loop over fragments

    //  }  // Close loop over the TRK and CAL collections

    if( diagLevel_ > 0 ) {
      std::cout << "mu2e::FragmentAna::produce exiting eventNumber=" << (int)(event.event()) << " / timestamp=" << (int)eventNumber <<std::endl;

    }
  }  // analyze()

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::FragmentAna);
