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
#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"
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
#include <stdbool.h> //temporary

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
    fhicl::Atom<std::string> rawFile {fhicl::Name("rawFile"), "raw.bin"};
    fhicl::Atom<std::string> zsFile {fhicl::Name("zsFile"), "zs.bin"};
    fhicl::Atom<std::string> phFile {fhicl::Name("phFile"), "ph.bin"};
    fhicl::Atom<std::string> rawHeaderFile {fhicl::Name("rawHeaderFile"), "rawWithHeader.bin"}; 
    fhicl::Atom<std::string> eventFile {fhicl::Name("eventFile"), "event.bin"};
  };
  
  explicit STMDigisFromFragments(const art::EDProducer::Table<Config>& config); // constructor created, config via fcl
  virtual ~STMDigisFromFragments(); //declares destructor

  virtual void produce(Event &) override;
  void endJob() override;//Final prinout summary

private:
  art::InputTag _stmFragmentsTag;
  std::ofstream _rawOut; //Files which can be acessesed throughout
  std::ofstream _zsOut;
  std::ofstream _phOut;
  std::ofstream _rawHeaderOut;
  std::ofstream _eventOut;

  //Metrics
  size_t _totalEvents{0}; //Another way to initialize to zero
  size_t _totalFragments{0};
  size_t _totalContainers{0};
  size_t _totalInner{0};
  size_t _totalRaw{0};
  size_t _totalZS{0};
  size_t _totalPH{0};

  //Additional metrics
  size_t _totalZeroRaw{0};
  size_t _totalZeroZS{0};
  size_t _totalZeroPH{0};
  size_t _totalEmptyRaw{0};
  size_t _totalEmptyZS{0};
  size_t _totalEmptyPH{0};
}; // STMDigisFromFragments

// ======================================================================


STMDigisFromFragments::STMDigisFromFragments(const art::EDProducer::Table<Config>& config)
  : art::EDProducer{config}
  ,_stmFragmentsTag(config().stmTag())

{
  // Set the size of the vector
  produces<mu2e::STMWaveformDigiCollection>("raw");//Waveforms
  produces<mu2e::STMWaveformDigiCollection>("zs");
  produces<mu2e::STMWaveformDigiCollection>("rawWithHeader");
  produces<mu2e::STMMWDDigiCollection>("ph"); // digi series

  //turns file into binary
  _rawOut.open(config().rawFile(), std::ios::binary);
  _zsOut.open(config().zsFile(), std::ios::binary);
  _phOut.open(config().phFile(), std::ios::binary);
  _rawHeaderOut.open(config().rawHeaderFile(), std::ios::binary);
  _eventOut.open(config().eventFile(), std::ios::binary);

  //Check to make sure we are reading the file
  if(!_rawOut || !_zsOut || !_phOut || !_rawHeaderOut || !_eventOut){
    throw cet::exception("FILEOPEN")<< "Failed to open one or more output files\n";
  }
  
}

STMDigisFromFragments::~STMDigisFromFragments(){
  if ( _rawOut.is_open() )
    _rawOut.close();
  if ( _zsOut.is_open() )
    _zsOut.close();
  if ( _phOut.is_open() )
    _phOut.close();
  if( _rawHeaderOut.is_open() )
    _rawHeaderOut.close();
  if( _eventOut.is_open() )
    _eventOut.close();    
}//Closing files

// ----------------------------------------------------------------------

void STMDigisFromFragments::produce(Event& event)
{
  
  ++_totalEvents; //Increment Total Event Counter
  
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> zs_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMMWDDigiCollection> ph_digis(new mu2e::STMMWDDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_header_waveform_digis(new mu2e::STMWaveformDigiCollection);
  //std::unique_ptr<mu2e::STMWaveformDigiCollection> ph_waveform_digis(new mu2e::STMWaveformDigiCollection);//Original

  art::Handle<artdaq::Fragments> STMFragmentsH;
  event.getByLabel(_stmFragmentsTag, STMFragmentsH);
  const auto STMFragments = STMFragmentsH.product();

  auto writePayload = [](std::ofstream& out, // [] is a capture list that only works with these parameters
                         const int16_t* data,
                         size_t words){
    out.write(reinterpret_cast<const char*>(data),//Writing to file stream in binary
              words * sizeof(int16_t));
  };

  //Event Metrics
  size_t localRaw_frags{0};
  size_t localZS_frags{0};
  size_t localPH_frags{0};
  size_t zeroRaw_frags{0}; 
  size_t zeroZS_frags{0};
  size_t zeroPH_frags{0};
  size_t emptyRaw_frags{0};
  size_t emptyZS_frags{0};
  size_t emptyPH_frags{0};

  //loop over frags
  for (const auto& frag : *STMFragments) {
    ++_totalFragments; //Increment Total Frag counter
    
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
          ++_totalRaw; //Increment job counter
	  ++localRaw_frags;//Increment event counter
	    
          auto payloadPtr = stm_frag.payloadBegin();
          auto payloadWords = stm_frag.payloadWords();
          bool allZeros = true;//Assume raw frag is zero filled
	  
          //Checks for empty frag
          if (payloadWords == 0) {
            std::cout<< "\nFound an empty frag, i = " << i <<" @Raw\n";
	    ++_totalEmptyRaw;//Job counter
	    ++emptyRaw_frags;//Event counter
            continue;//stops this loop check and goes to next frag
            }

          //Check if any data points are not zero
          for (size_t k =0; k < payloadWords; ++k){
            if (payloadPtr[k] != 0){
              allZeros = false;
              break;
            }
          }

          //Chekcs if zero filled
          if (allZeros){
            std::cout << "\nFound a zero filled frag, i = " << i << " @Raw\n";
	    ++_totalZeroRaw;//Increment counters
	    ++zeroRaw_frags;//Increment Zero Raw Counter
            continue;
          }

          //Print first 20 values 
          std::cout << "\nFirst 20 adcs: ";
          for (size_t kk = 0; kk < std::min<size_t>(payloadWords,20); ++kk){
            std::cout << payloadPtr[kk] << " ,";
          }
          std::cout << "\nRaw header : ZS Length = " << stm_frag.zsLength()
		    << " , ZS Regions  = " << stm_frag.zsRegions()
		    << " , i = " << i << " @Raw\n";

          //Ideally only good frags get up to here
          //------full data (header + payload)
          { 
            //Waveform with data creation
            auto dataPtr = stm_frag.dataBegin();//Inner variable
            auto dataWords = stm_frag.dataWords();
            stm_waveform.set_data(dataWords, dataPtr);
            raw_header_waveform_digis->emplace_back(stm_waveform);
            writePayload( _rawHeaderOut , dataPtr, dataWords); 
          }

          //----- payload-only
          {
	    //auto payloadPtr = stm_frag.payloadBegin();
            //auto payloadWords = stm_frag.payloadWords();
            stm_waveform.set_data(payloadWords, payloadPtr);
            raw_waveform_digis->emplace_back(stm_waveform);
            writePayload(_rawOut, payloadPtr, payloadWords);
          }

        }//End of isRaw
	
        else if (stm_frag.isZS()){
          ++_totalZS; //Incremenet ZS counter
	  ++localZS_frags;
	  
	  auto payloadPtr = stm_frag.payloadBegin();
	  auto payloadWords = stm_frag.payloadWords();
	  bool allZeros = true; //assumes all adcs are zero

	  //Check if payload is empty
	  if (payloadWords == 0) {
	    std::cout << "\nFound an empty frag, i = " << i<< " @ZS\n";
	    ++_totalEmptyZS;
	    ++emptyZS_frags; //Increment empty ZS counter
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
	    std::cout << "\nFound a zero filled frag, i = "<< i<< " @ZS\n";
	    ++_totalZeroZS;
	    ++zeroZS_frags;
	    continue;
	  }
	  //Print first 20 payload adcs
	  std::cout << "\nFirst 20 adcs: ";
	  for (size_t kk = 0 ; kk < std::min<size_t>(payloadWords,20) ; ++kk){
	    std::cout << payloadPtr[kk] << " , ";
	  }
	  std::cout << "i = " << i << " @ZS\n";

	  //Defintions for payload references
	  auto dataPtr = stm_frag.dataBegin();
	  auto dataWords = stm_frag.dataWords();
	  auto dataEnd = dataPtr + dataWords;
	  size_t seg = 0;
	  size_t totalLen = 0;
	  uint16_t lastZSindex = 0; //keeps track of last recorded index from header -> with respect to what?
	  uint16_t lastLen = 0; //keeps track of last recoded length from header

	  writePayload(_zsOut, payloadPtr, payloadWords);
	  
	  while (dataPtr + 2 <= dataEnd){
	    uint16_t current_zs_location = static_cast<uint16_t>(dataPtr[0]);
	    uint16_t current_zs_size = static_cast<uint16_t>(dataPtr[1]);
	    auto adc = dataPtr + 2;
	    if (adc + current_zs_size > dataEnd)
	      break;

	    uint32_t trigTimeOffset = current_zs_location;
	    std::vector<int16_t> segADCS(adc, adc +  current_zs_size); //1D array, contains adcs to this_zs_Size - 1
	    mu2e::STMWaveformDigi stm_waveform(trigTimeOffset, segADCS); //New constructore use
	    zs_waveform_digis->emplace_back(stm_waveform); //emplace
	    
	    //A print check per segment
	    std::cout << "Region = " << seg << " , zs_location = " << current_zs_location << " , current_ZS_size = " << current_zs_size
		      << " , trigTimeOffset = " << trigTimeOffset << "\n" ;
	    
	    lastZSindex = current_zs_location;
	    lastLen = current_zs_size;
	    totalLen += lastLen;
	    ++seg;

	    dataPtr = adc + current_zs_size;
	  }
     	  //Summary
	  std::cout << "ZS Regions = " << seg
		    << " , lastZSindex = " << lastZSindex
		    << " , lastZSLen = " << lastLen
		    << " , ZS total length = " << totalLen
		    << " , i =  " << i <<  " @ZS\n";
	  
        }//End of isZS
	
        else if (stm_frag.isMWD()){
	  ++_totalPH;
	  ++localPH_frags;
	  //Check if zero filled
	  auto payloadPtr = stm_frag.payloadBegin();
	  auto payloadWords = stm_frag.payloadWords();
	  bool allZeros = true;

	  if (payloadWords ==0){
	    ++_totalEmptyPH;
	    ++emptyPH_frags;
	    std::cout << "\nFound an empty frag, i = " << i<< " @PH\n";
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
	    ++_totalZeroPH;
	    ++zeroPH_frags;
	    std::cout<< "\nFound a zero filled frag, i = " << i<< " @PH\n";
	    continue;
	  }

	  std::cout << "\nFound a good frag, i = " << i <<" @PH\n";
	  size_t digiWords = stm_frag.payloadWords();	  
          auto const* digiPtr = stm_frag.payloadBegin();

	  writePayload(_phOut,digiPtr, digiWords);
	  
	  for (size_t i_PH = 0; i_PH < digiWords ; ++i_PH){
	    int16_t PH = digiPtr[i_PH];
	    mu2e::STMMWDDigi PH_digi(0, PH);
	    ph_digis->emplace_back(PH_digi);
	  }

      }//End of isMWD and is checks

        //---Combined stream write w. order preserved ----
        {
          const int16_t* cptr = nullptr;
          size_t cwords = 0;

          if (stm_frag.isRaw()){
            cptr = stm_frag.dataBegin();
            cwords = stm_frag.dataWords();
          }
          else if (stm_frag.isZS()){
            cptr = stm_frag.payloadBegin();
            cwords = stm_frag.payloadWords();
          }
          else if (stm_frag.isMWD()){
            cptr = stm_frag.payloadBegin();
            cwords = stm_frag.payloadWords();
          }

          if (cptr != nullptr && cwords >0) {
            _eventOut.write(reinterpret_cast<const char*>(cptr),
                            cwords * sizeof(int16_t) );
          }
        }

      }//End Container loop
    } else {
      //fallback (non-container case)
      mu2e::STMFragment stm_frag(frag);
      auto ptr = stm_frag.payloadBegin();
      auto words = stm_frag.payloadWords();

      if (stm_frag.isRaw()) {
        writePayload(_rawOut, ptr, words);
      }
    }
  } //End of frags loop---George suggestions

  //Event Summary -> tells us what happens per event
  std::cout << "\n========== STM EVENT SUMMARY ==========\n";
  std::cout << "Extracted Raw waveforms     : "<< raw_waveform_digis->size() <<"\n";
  std::cout << "Extracted ZS waveforms      : "<< zs_waveform_digis->size() <<"\n";
  std::cout << "Extracted PH digis          : " << ph_digis->size() <<"\n";

  std::cout << "\n--- Frags Read ---\n";
  std::cout << "Raw frags   : " << localRaw_frags << "\n";
  std::cout << "ZS frags    : " << localZS_frags << "\n";
  std::cout << "PH frags    : " << localPH_frags << "\n";
  
  std::cout << "\n--- Filter results ---\n";
  std::cout << "Zero Raw frags   : " << zeroRaw_frags << "\n";
  std::cout << "Zero ZS frags    : " << zeroZS_frags << "\n";
  std::cout << "Zero PH frags    : " << zeroPH_frags << "\n";

  std::cout << "Empty Raw frags  : " << emptyRaw_frags <<"\n";
  std::cout << "Empty ZS frags   : " << emptyZS_frags << "\n";
  std::cout << "Empty PH frags   : " << emptyPH_frags << "\n";
  
  std::cout << "=================================\n";

  //Final move
  event.put(std::move(raw_waveform_digis), "raw");
  event.put(std::move(zs_waveform_digis), "zs");
  event.put(std::move(ph_digis), "ph");
  event.put(std::move(raw_header_waveform_digis), "rawWithHeader");
    
} // produce()

// ======================================================================


void STMDigisFromFragments::endJob() {

  //Tells us what happened at the very end
  std::cout << "\n========== STM JOB SUMMARY ==========\n";

  std::cout << "Total events       : " << _totalEvents << "\n";
  std::cout << "Total fragments    : " << _totalFragments << "\n";
  std::cout << "Container frags    : " << _totalContainers << "\n";
  std::cout << "Inner fragments    : " << _totalInner << "\n";

  std::cout << "\n--- Data types read ---\n";
  std::cout << "RAW               : " << _totalRaw << "\n";
  std::cout << "ZS                : " << _totalZS << "\n";
  std::cout << "PH                : " << _totalPH << "\n";

  std::cout << "\n-- -Data types filtered ---\n";
  std::cout << "Zero RAW frags    : " << _totalZeroRaw << "\n";
  std::cout << "Zero ZS frags     : " << _totalZeroZS << "\n";
  std::cout << "Zero PH frags     : " << _totalZeroPH << "\n";
  std::cout << "Empty Raw frags   : " << _totalEmptyRaw << "\n";
  std::cout << "Empty ZS frags    : " << _totalEmptyZS << "\n";
  std::cout << "Empty PH frags    : " << _totalEmptyPH << "\n";


  std::cout << "=================================\n";
}


DEFINE_ART_MODULE(STMDigisFromFragments)
