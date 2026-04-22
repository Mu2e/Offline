// =====================================================================
//
// STMDigisFromFragments: create all types of STMDigis from STMFragments
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <list>
#include <numeric>
#include <random>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TLine.h>
#include <TGraph.h>
#include <stdbool.h>

namespace art
{
  class STMDigisFromFragments;
}

using art::STMDigisFromFragments;

class art::STMDigisFromFragments : public EDProducer
{
public:
  struct Config {
    fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("Input module")};
    fhicl::Atom<bool> saveRawWithHeaderFile {fhicl::Name("saveRawWithHeaderFile"), false};
    fhicl::OptionalAtom<int> verbosityLevel{fhicl::Name("verbosityLevel"), fhicl::Comment("Verbosity level")};
  };
  
  explicit STMDigisFromFragments(const art::EDProducer::Table<Config>& config); // constructor created, config via fcl
  virtual void produce(Event &) override;
  void endJob() override;//Final prinout summary

private:
  art::InputTag _stmFragmentsTag;
  
  //Metrics-> Will change this soon 
  size_t _totalEvents{0};
  size_t _totalFragments{0};
  size_t _totalContainers{0};
  size_t _totalInner{0};
  size_t _totalRaw{0};
  size_t _totalZS{0};
  size_t _totalPH{0};
  size_t _totalZeroRaw{0};
  size_t _totalZeroZS{0};
  size_t _totalZeroPH{0};
  size_t _totalEmptyRaw{0};
  size_t _totalEmptyZS{0};
  size_t _totalEmptyPH{0};

  size_t _unreadInnerFrags{0};

  //Additional metrics
  //HPGe
  size_t _totalRawHPGe{0};
  size_t _totalZSHPGe{0};
  size_t _totalPHHPGe{0};
  size_t _totalZeroRawHPGe{0};
  size_t _totalZeroZSHPGe{0};
  size_t _totalZeroPHHPGe{0};
  size_t _totalEmptyRawHPGe{0};
  size_t _totalEmptyZSHPGe{0};
  size_t _totalEmptyPHHPGe{0};

  //LaBr
  size_t _totalRawLaBr{0};
  size_t _totalZSLaBr{0};
  size_t _totalPHLaBr{0};
  size_t _totalZeroRawLaBr{0};
  size_t _totalZeroZSLaBr{0};
  size_t _totalZeroPHLaBr{0};
  size_t _totalEmptyRawLaBr{0};
  size_t _totalEmptyZSLaBr{0};
  size_t _totalEmptyPHLaBr{0};

  //fhicl varibales
  int _verbosityLevel = 0;
  bool _saveRawWithHeaderFile{false};
}; // STMDigisFromFragments

// ======================================================================


STMDigisFromFragments::STMDigisFromFragments(const art::EDProducer::Table<Config>& config)
  : art::EDProducer{config}
  ,_stmFragmentsTag(config().stmTag())
  ,_saveRawWithHeaderFile(config().saveRawWithHeaderFile())
  ,_verbosityLevel(config().verbosityLevel() ? *(config().verbosityLevel()) : 0)

{
  // Set the size of vectors for HPGe
  produces<mu2e::STMWaveformDigiCollection>("rawHPGe");//Waveforms
  produces<mu2e::STMWaveformDigiCollection>("zsHPGe");
  produces<mu2e::STMPHDigiCollection>("phHPGe"); // digi series
  if(_saveRawWithHeaderFile){ produces<mu2e::STMWaveformDigiCollection>("rawWithHeaderHPGe");}

  //Set the size of vectors for LaBr
  produces<mu2e::STMWaveformDigiCollection>("rawLaBr");
  produces<mu2e::STMWaveformDigiCollection>("zsLaBr");
  produces<mu2e::STMPHDigiCollection>("phLaBr"); // digi series                                                                                                          
  if(_saveRawWithHeaderFile){ produces<mu2e::STMWaveformDigiCollection>("rawWithHeaderLaBr");}

}

// ----------------------------------------------------------------------

void STMDigisFromFragments::produce(Event& event)
{
  
  ++_totalEvents; //Increment Total Event Counter

  //HPGe
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_HPGe_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> zs_HPGe_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMPHDigiCollection> ph_HPGe_digis(new mu2e::STMPHDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_HPGe_header_waveform_digis(new mu2e::STMWaveformDigiCollection);
  //std::unique_ptr<mu2e::STMWaveformDigiCollection> ph_waveform_digis(new mu2e::STMWaveformDigiCollection);//Original

  //LaBr
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_LaBr_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> zs_LaBr_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMPHDigiCollection> ph_LaBr_digis(new mu2e::STMPHDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_LaBr_header_waveform_digis(new mu2e::STMWaveformDigiCollection);
  //std::unique_ptr<mu2e::STMWaveformDigiCollection> ph_waveform_digis(new mu2e::STMWaveformDigiCollection);//Original   
  
  art::Handle<artdaq::Fragments> STMFragmentsH;
  event.getByLabel(_stmFragmentsTag, STMFragmentsH);
  const auto STMFragments = STMFragmentsH.product();

  //Event Metrics
  
  //HPGe
  size_t localRawHPGe_frags{0};
  size_t localZSHPGe_frags{0};
  size_t localPHHPGe_frags{0};
  size_t zeroRawHPGe_frags{0}; 
  size_t zeroZSHPGe_frags{0};
  size_t zeroPHHPGe_frags{0};
  size_t emptyRawHPGe_frags{0};
  size_t emptyZSHPGe_frags{0};
  size_t emptyPHHPGe_frags{0};
  
  //LaBr
  size_t localRawLaBr_frags{0};
  size_t localZSLaBr_frags{0};
  size_t localPHLaBr_frags{0};
  size_t zeroRawLaBr_frags{0};
  size_t zeroZSLaBr_frags{0};
  size_t zeroPHLaBr_frags{0};
  size_t emptyRawLaBr_frags{0};
  size_t emptyZSLaBr_frags{0};
  size_t emptyPHLaBr_frags{0};

  //Additional
  size_t unread_InnerFrags{0};
  uint16_t ContainerFragID{0};
  
  //loop over frags
  for (const auto& frag : *STMFragments) {
    ++_totalFragments; //Increment Total Frag counter
    bool isLaBr{false};
    bool isHPGe{false};

    uint16_t ZSfromRaw{0};
    uint16_t expectedZSRegions{0};
    bool readRawZSinfo{false};

    ContainerFragID = frag.fragmentID();
    if (_verbosityLevel >= 3){std::cout << "\nFrag_id : " << ContainerFragID << "\n";}

    if(ContainerFragID == 103) {isHPGe = true;} else if (ContainerFragID == 203){isLaBr = true ;}
    
    //Check if this is a container fragment
    if (frag.type() == artdaq::Fragment::ContainerFragmentType){

      artdaq::ContainerFragment cont_frag(frag);
      ++_totalContainers;
      size_t blocks = cont_frag.block_count();
      _totalInner += blocks;
      //loop over container
      // i index corresponds to inner frag 
      for (size_t i = 0; i < cont_frag.block_count(); ++i){

        auto inner_frag = cont_frag.at(i);
        mu2e::STMFragment stm_frag(*inner_frag);
        mu2e::STMWaveformDigi stm_waveform;

        //auto ptr = stm_frag.payloadBegin();//Outer variable
        //auto words = stm_frag.payloadWords();

        if (stm_frag.isRaw()) {
	  
	  //Job Counter
          ++_totalRaw;	
	  //Conitional Job and Event Counter
	  if(isHPGe){ ++_totalRawHPGe ; ++localRawHPGe_frags;} else if (isLaBr){++_totalRawLaBr ; ++localRawLaBr_frags;}
	  
          auto payloadPtr = stm_frag.payloadBegin();
          auto payloadWords = stm_frag.payloadWords();
          bool allZeros = true;//Assume raw frag is zero filled
	  
          //Checks for empty frag
          if (payloadWords == 0) {
	    if(_verbosityLevel >=3){std::cout<< "\nFound an empty frag, i = " << i <<" @Raw\n";}
	    ++_totalEmptyRaw;//Job counter
	    //Eventcounter
	    if(isHPGe){ ++_totalEmptyRawHPGe; ++emptyRawHPGe_frags;}else if (isLaBr){ ++_totalEmptyRawLaBr; ++emptyRawLaBr_frags;}
	    continue;//stops this loop check and goes to next frag
            }

          //Check if any data points are not zero
          for (size_t k =0; k < payloadWords; ++k){
            if (payloadPtr[k] != 0){
              allZeros = false;
              break;
            }
          }

          //Check if zero filled
          if (allZeros){
	    if (_verbosityLevel >=3){std::cout << "\nFound a zero filled frag, i = " << i << " @Raw\n";}
	    ++_totalZeroRaw;//Increment counters
	    //Eventcounter                                                                                                                                                
            if(isHPGe){++_totalZeroRawHPGe; ++zeroRawHPGe_frags;} else if (isLaBr){ ++_totalZeroRawLaBr; ++zeroRawLaBr_frags;}
            continue;
          }

          //Print first 20 values
	  if (_verbosityLevel >=3){std::cout << "\nFound a good frag, i = " << i <<" @Raw\n";}

	  if (_verbosityLevel >= 4){
	    std::cout << "\nFirst 20 adcs: ";
	    for (size_t kk = 0; kk < std::min<size_t>(payloadWords,20); ++kk){
	      std::cout << payloadPtr[kk] << " ,";
	    }
	    if(_verbosityLevel >=5){
	      std::cout << "\nRaw header : Raw Length = " << stm_frag.rawLength()
			<<" , ZS Length = " << stm_frag.zsLength()
			<< " , ZS Regions  = " << stm_frag.zsRegions()
			<< " , i = " << i << " @Raw\n";
	    }
	  }
	  expectedZSRegions = stm_frag.zsRegions();
	  ZSfromRaw = stm_frag.zsLength();//Stores ZS length on outside variable
	  readRawZSinfo = true;//Shows we were able to read the Raw header information
	  
          //Ideally only good frags get up to here
	  if (isHPGe){
	    if (_saveRawWithHeaderFile){
	      auto dataPtr = stm_frag.dataBegin();
	      auto dataWords = stm_frag.dataWords();
	      stm_waveform.set_data(dataWords, dataPtr);
	      raw_HPGe_header_waveform_digis->emplace_back(stm_waveform);
	    }

	    {
	      stm_waveform.set_data(payloadWords, payloadPtr);
	      raw_HPGe_waveform_digis->emplace_back(stm_waveform);
	    }

	  } else if (isLaBr){
            if (_saveRawWithHeaderFile){
              auto dataPtr = stm_frag.dataBegin();
              auto dataWords = stm_frag.dataWords();
              stm_waveform.set_data(dataWords, dataPtr);
              raw_LaBr_header_waveform_digis->emplace_back(stm_waveform);
            }

            {
              stm_waveform.set_data(payloadWords, payloadPtr);
              raw_LaBr_waveform_digis->emplace_back(stm_waveform);
            }
	    
	  }//End of if CFID=103,203
	  
        }//End of isRaw
	
        else if (stm_frag.isZS()){
          ++_totalZS; //Incremenet ZS counter
          if (isHPGe){ ++_totalZSHPGe; ++localZSHPGe_frags; } else if (isLaBr){ ++_totalZSLaBr; ++localZSLaBr_frags;}
	  	  
	  auto payloadPtr = stm_frag.payloadBegin();
	  auto payloadWords = stm_frag.payloadWords();
	  bool allZeros = true; //assumes all adcs are zero

	  //Check if payload is empty
	  if (payloadWords == 0) {
	    if (_verbosityLevel >=3){std::cout << "\nFound an empty frag, i = " << i << " @ZS\n";}
	    ++_totalEmptyZS;
	    if (isHPGe){ ++_totalEmptyZSHPGe; ++emptyZSHPGe_frags; } else if (isLaBr){ ++_totalEmptyZSLaBr; ++emptyZSLaBr_frags;}
	    continue;
	  }
	  //Check if any adc are non-zero
	  for (size_t k = 0 ; k < payloadWords ; ++k){
	    if(payloadPtr[k] != 0 ){
	      allZeros = false;
	      break;
	    }
	  }
	  //Check if zero filled
	  if (allZeros) {
	    if (_verbosityLevel >=3){std::cout << "\nFound a zero filled frag, i = "<< i << " @ZS\n";}
	    ++_totalZeroZS;
	    if (isHPGe){ ++_totalZeroZSHPGe; ++zeroZSHPGe_frags; } else if (isLaBr){ ++_totalZeroZSLaBr; ++zeroZSLaBr_frags;}
	    continue;
	  }
	  //Print first 20 payload adcs
	  if (_verbosityLevel >=3){std::cout << "\nFound a good frag, i = " << i << " @ZS\n";}

	  if (_verbosityLevel >=4){
	    std::cout << "\nFirst 20 adcs: ";
	    for (size_t kk = 0 ; kk < std::min<size_t>(payloadWords,20) ; ++kk){
	      std::cout << payloadPtr[kk] << " , ";
	    }
	    std::cout << "i = " << i << " @ZS\n";
	  }
	  
	  //Defintions for payload references
	  auto dataPtr = stm_frag.dataBegin();
	  auto dataWords = stm_frag.dataWords();
	  auto dataEnd = dataPtr + dataWords;
	  size_t seg = 0;
	  size_t totalLen = 0;
	  uint16_t lastZSindex = 0; //keeps track of last recorded index from header -> with respect to what?
	  uint16_t lastLen = 0; //keeps track of last recoded length from header

	  if(_verbosityLevel >= 6){std::cout << "dataWords  : " << dataWords
					     << " dataWords%4 : " << dataWords%4
					     <<" @ZS" << "\n"; }
	  
	  while (dataPtr + 2 <= dataEnd){
	    if ( readRawZSinfo && seg >= expectedZSRegions ) break;
	    uint16_t current_zs_location = static_cast<uint16_t>(dataPtr[0]);
	    uint16_t current_zs_size = static_cast<uint16_t>(dataPtr[1]);
	    auto adc = dataPtr + 2;
	    if (adc + current_zs_size > dataEnd)
	      break;

	    uint32_t trigTimeOffset = current_zs_location;
	    std::vector<int16_t> segADCS(adc, adc +  current_zs_size); //1D array, contains adcs to this_zs_Size - 1
	    mu2e::STMWaveformDigi stm_waveform(trigTimeOffset, segADCS); //New constructore use

	    //emplacing 
	    if (isHPGe){zs_HPGe_waveform_digis->emplace_back(stm_waveform);
	    } else if (isLaBr){zs_LaBr_waveform_digis->emplace_back(stm_waveform);
	    }

	    if (_verbosityLevel >=6){
	      //A print check per segment
	      std::cout << "Region = " << seg << " , zs_index = " << current_zs_location << " , zs_size = " << current_zs_size
			<< " , trigTimeOffset = " << trigTimeOffset << "\n" ;
	    }
	    
	    //Update variables
	    lastZSindex = current_zs_location;
	    lastLen = current_zs_size;
	    totalLen += lastLen;
	    ++seg;
	    dataPtr = adc + current_zs_size;
	  }

	  if (_verbosityLevel >=5){
	    //Summary
	    std::cout << "ZS Regions = " << seg
		      << " , lastZSindex = " << lastZSindex
		      << " , lastZSLen = " << lastLen
		      << " , ZS total length = " << totalLen
		      << " , i =  " << i <<  " @ZS\n";
	  }

	  //Throw out if ZSLengthfromRaw != totalLen
	  if (readRawZSinfo && ZSfromRaw != totalLen){
	    throw cet::exception("STM_UNPACKING")
	      << "\n=== ZS Length mismatch ===\n"
	      << "ZS length from Raw header : " << ZSfromRaw << "\n"
	      << "ZS length calculated from file : " << totalLen << "\n"
	      << "Found at inner frag i : " << i <<"\n"
	      << "Encountered at event : " << _totalEvents <<"\n"      ;
	  }
	  
        }//End of isZS
	
        else if (stm_frag.isPH()){
	  
	  ++_totalPH;
	  if (isHPGe){ ++_totalPHHPGe; ++localPHHPGe_frags; } else if (isLaBr){ ++_totalPHLaBr; ++localPHLaBr_frags; }

	  //Check if zero filled
	  auto payloadPtr = stm_frag.payloadBegin();
	  auto payloadWords = stm_frag.payloadWords();
	  bool allZeros = true;

	  if (payloadWords ==0){
            if (_verbosityLevel >= 3){std::cout << "\nFound an empty frag, i = " << i<< " @PH\n";}
	    ++_totalEmptyPH;
	    if (isHPGe){ ++_totalEmptyPHHPGe; ++emptyPHHPGe_frags; } else if (isLaBr){ ++_totalEmptyPHLaBr; ++emptyPHLaBr_frags; }
	    continue;
	  }

	  //Check if zero filled
	  for (size_t k = 0; k< payloadWords; ++k) {
	    if (payloadPtr[k] !=0){
	      allZeros = false;
	      break;
	    }
	  }
	  //count if zero filled
	  if (allZeros){
	    if (_verbosityLevel >=3){ std::cout<< "\nFound a zero filled frag, i = " << i<< " @PH\n";}
	    ++_totalZeroPH;
	    if (isHPGe){ ++_totalZeroPHHPGe; ++zeroPHHPGe_frags; } else if (isLaBr){ ++_totalZeroPHLaBr; ++zeroPHLaBr_frags; }
	    continue;
	  }

	  if (_verbosityLevel >= 3){std::cout << "\nFound a good frag, i = " << i <<" @PH\n";}
	  
	  size_t digiWords = stm_frag.payloadWords();	  
          auto const* digiPtr = stm_frag.payloadBegin();

	  for (size_t i_PH = 0; i_PH < digiWords ; ++i_PH){
	    int16_t PH = digiPtr[i_PH];
	    mu2e::STMPHDigi PH_digi(0, PH);

	    if(isHPGe){ ph_HPGe_digis->emplace_back(PH_digi);
	    } else if (isLaBr){ ph_LaBr_digis->emplace_back(PH_digi);
	    }
	    
	  }

	}//End of isPH and is checks
	else{
	  if(_verbosityLevel >=3){
	    //std::cout << "\n Inner frag i = " << i
	    //	      << "\nDid not read Raw/ZS/PH\n"
	    //	      << "Inner Frag_id : " << inner_frag->fragmentID() << "\n" ;

	    ++ _unreadInnerFrags; //For Job summary
	    ++ unread_InnerFrags; //For event summary


	  }
	}
      }
    } else {
      //fallback (non-container case)
      mu2e::STMFragment stm_frag(frag);
      auto ptr = stm_frag.payloadBegin();
      auto words = stm_frag.payloadWords();

      //if (stm_frag.isRaw()){}	//        writePayload(_rawOut, ptr, words);
    }
  } //End of frags loop---George suggestions

  if (_verbosityLevel >= 2 ){ 
    //Event Summary -> tells us what happens per event
    std::cout << "\n========== STM EVENT SUMMARY - (Unpacking Module) ==========\n";
    std::cout << "Extracted Raw waveforms     : " << raw_HPGe_waveform_digis->size() <<"\n";
    std::cout << "Extracted ZS waveforms      : " << zs_HPGe_waveform_digis->size() <<"\n";
    std::cout << "Extracted PH digis          : " << ph_HPGe_digis->size() <<"\n";

    std::cout << "Extracted Raw waveforms     : " << raw_LaBr_waveform_digis->size() <<"\n";
    std::cout << "Extracted ZS waveforms      : " << zs_LaBr_waveform_digis->size() <<"\n";
    std::cout << "Extracted PH digis          : " << ph_LaBr_digis->size() <<"\n";

    
    std::cout << "\n--- Frags Read ---\n";
    std::cout << "Raw HPGE frags         : " << localRawHPGe_frags << "\n";
    std::cout << "ZS  HPGe frags         : " << localZSHPGe_frags << "\n";
    std::cout << "PH  HPGe frags         : " << localPHHPGe_frags << "\n";
    
    std::cout << "Raw LaBr frags         : " << localRawLaBr_frags << "\n";
    std::cout << "ZS  LaBr frags         : " << localZSLaBr_frags << "\n";
    std::cout << "PH  LaBr frags         : " << localPHLaBr_frags << "\n";

  
    std::cout << "\n--- Filter results ---\n";
    std::cout << "Zero Raw  HPGe frags   : " << zeroRawHPGe_frags << "\n";
    std::cout << "Zero ZS   HPGe frags   : " << zeroZSHPGe_frags << "\n";
    std::cout << "Zero PH   HPGe frags   : " << zeroPHHPGe_frags << "\n";
    
    std::cout << "Zero Raw  LaBr frags   : " << zeroRawLaBr_frags << "\n";
    std::cout << "Zero ZS   LaBr frags   : " << zeroZSLaBr_frags << "\n";
    std::cout << "Zero PH   LaBr frags   : " << zeroPHLaBr_frags << "\n";


    std::cout << "Empty Raw HPGe frags   : " << emptyRawHPGe_frags <<"\n";
    std::cout << "Empty ZS  HPGe frags   : " << emptyZSHPGe_frags << "\n";
    std::cout << "Empty PH  HPGe frags   : " << emptyPHHPGe_frags << "\n";

    std::cout << "Empty Raw LaBr frags   : " << emptyRawLaBr_frags <<"\n";
    std::cout << "Empty ZS  LaBr frags   : " << emptyZSLaBr_frags << "\n";
    std::cout << "Empty PH  LaBr frags   : " << emptyPHLaBr_frags << "\n";

    std::cout << "Unread frags           : " << unread_InnerFrags << "\n";
  
    std::cout << "=================================\n";
  }
  
  //Final move
  //HPGe
  event.put(std::move(raw_HPGe_waveform_digis), "rawHPGe");
  event.put(std::move(zs_HPGe_waveform_digis), "zsHPGe");
  event.put(std::move(ph_HPGe_digis), "phHPGe");
  if (_saveRawWithHeaderFile){ event.put(std::move(raw_HPGe_header_waveform_digis), "rawWithHeaderHPGe"); }


  //LaBr
  event.put(std::move(raw_LaBr_waveform_digis), "rawLaBr");
  event.put(std::move(zs_LaBr_waveform_digis), "zsLaBr");
  event.put(std::move(ph_LaBr_digis), "phLaBr");
  if (_saveRawWithHeaderFile){ event.put(std::move(raw_LaBr_header_waveform_digis), "rawWithHeaderLaBr"); }
    
} // produce()

// ======================================================================


void STMDigisFromFragments::endJob() {
  if (_verbosityLevel >= 1){
    //Tells us what happened at the very end
    std::cout << "\n========== STM JOB SUMMARY - (Unpacking Module) ==========\n";

    std::cout << "Total events             : " << _totalEvents << "\n";
    std::cout << "Total frags              : " << _totalFragments << "\n";
    std::cout << "Container frags          : " << _totalContainers << "\n";
    std::cout << "Inner frags              : " << _totalInner << "\n";

    std::cout << "\n--- Data types read ---\n";
    std::cout << "RAW Total                : " << _totalRaw << "\n";
    std::cout << "ZS  Total                : " << _totalZS << "\n";
    std::cout << "PH  Total                : " << _totalPH << "\n";

    std::cout << "RAW LaBr                 : " << _totalRawLaBr << "\n";
    std::cout << "ZS  LaBr                 : " << _totalZSLaBr << "\n";
    std::cout << "PH  LaBr                 : " << _totalPHLaBr << "\n";
    
    std::cout << "RAW HPGe                 : " << _totalRawHPGe << "\n";
    std::cout << "ZS  HPGe                 : " << _totalZSHPGe << "\n";
    std::cout << "PH  HPGe                 : " << _totalPHHPGe << "\n";

    std::cout << "\n--- Data types filtered ---\n";
    std::cout << "Zero RAW   HPGe frags    : " << _totalZeroRawHPGe << "\n";
    std::cout << "Zero ZS    HPGe frags    : " << _totalZeroZSHPGe << "\n";
    std::cout << "Zero PH    HPGe frags    : " << _totalZeroPHHPGe << "\n";
    std::cout << "Empty Raw  HPGe frags    : " << _totalEmptyRawHPGe << "\n";
    std::cout << "Empty ZS   HPGe frags    : " << _totalEmptyZSHPGe << "\n";
    std::cout << "Empty PH   HPGe frags    : " << _totalEmptyPHHPGe << "\n";
  
    std::cout << "Zero RAW   LaBr frags    : " << _totalZeroRawLaBr << "\n";
    std::cout << "Zero ZS    LaBr frags    : " << _totalZeroZSLaBr << "\n";
    std::cout << "Zero PH    LaBr frags    : " << _totalZeroPHLaBr << "\n";
    std::cout << "Empty Raw  LaBr frags    : " << _totalEmptyRawLaBr << "\n";
    std::cout << "Empty ZS   LaBr frags    : " << _totalEmptyZSLaBr << "\n";
    std::cout << "Empty PH   LaBr frags    : " << _totalEmptyPHLaBr << "\n";
    
    std::cout << "Unread frags             : " << _unreadInnerFrags << "\n";


    
    std::cout << "=================================\n";

  }
}

DEFINE_ART_MODULE(STMDigisFromFragments)
