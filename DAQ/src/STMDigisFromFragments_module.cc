// =====================================================================
//
// STMDigisFromFragments: create all types of STMDigis from STMFragments
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/RecoDataProducts/inc/STMFragmentSummary.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMPHDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/ContainerFragment.hh>
#include <artdaq-core/Data/Fragment.hh>
#include "canvas/Persistency/Common/Ptr.h"
//#include "canvas/Persistency/Common/Assns.h"

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
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
    fhicl::OptionalAtom<int> verbosityLevel{fhicl::Name("verbosityLevel"), fhicl::Comment("Verbosity level")};
    fhicl::Atom<bool> saveSTMFragSummary{fhicl::Name("saveSTMFragSummary"),false};

    fhicl::Atom<bool> saveRawWithHeaderWaveform_HPGe {fhicl::Name("saveRawWithHeaderWaveform_HPGe"), false};
    fhicl::Atom<bool> saveRawWaveform_HPGe {fhicl::Name("saveRawWaveform_HPGe"), false};
    fhicl::Atom<bool> saveZSWaveform_HPGe{fhicl::Name("saveZSWaveform_HPGe"), false};

    fhicl::Atom<bool> saveRawWithHeaderWaveform_LaBr{fhicl::Name("saveRawWithHeaderWaveform_LaBr"), false};
    fhicl::Atom<bool> saveRawWaveform_LaBr{fhicl::Name("saveRawWaveform_LaBr"), false};
    fhicl::Atom<bool> saveZSWaveform_LaBr{fhicl::Name("saveZSWaveform_LaBr"), false};

  };

  explicit STMDigisFromFragments(const art::EDProducer::Table<Config>& config); // constructor created, config via fcl
  virtual void produce(Event &) override;
  void endJob() override;

private:
  art::InputTag _stmFragmentsTag;

  //Metrics
  size_t _totalEvents{0};
  size_t _totalFragments{0};
  size_t _totalContainers{0};
  size_t _totalContainersHPGe{0};
  size_t _totalContainersLaBr{0};
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
  size_t _totalEventsWithHPGeLaBr{0};
  size_t _totalEventsWithOnlyHPGe{0};
  size_t _totalEventsWithOnlyLaBr{0};
  size_t _totalEventsWithNone{0};
  size_t _totalNonContainerFrags{0};

  //Additional metrics
  //HPGe
  size_t _totalRawHPGe{0};
  size_t _totalZSHPGe{0};
  size_t _totalPHHPGe{0};
  size_t _totalGoodRawHPGe{0};
  size_t _totalGoodZSHPGe{0};
  size_t _totalGoodPHHPGe{0};
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
  size_t _totalGoodRawLaBr{0};
  size_t _totalGoodZSLaBr{0};
  size_t _totalGoodPHLaBr{0};
  size_t _totalZeroRawLaBr{0};
  size_t _totalZeroZSLaBr{0};
  size_t _totalZeroPHLaBr{0};
  size_t _totalEmptyRawLaBr{0};
  size_t _totalEmptyZSLaBr{0};
  size_t _totalEmptyPHLaBr{0};

  //fhicl varibales
  bool _saveRawWithHeaderWaveform_HPGe{false};
  bool _saveRawWaveform_HPGe{false};
  bool _saveZSWaveform_HPGe{false};
  bool _saveRawWithHeaderWaveform_LaBr{false};
  bool _saveRawWaveform_LaBr{false};
  bool _saveZSWaveform_LaBr{false};
  bool _saveSTMFragSummary{false};
  int _verbosityLevel{0};

}; // STMDigisFromFragments

// ======================================================================


STMDigisFromFragments::STMDigisFromFragments(const art::EDProducer::Table<Config>& config)
  : art::EDProducer{config}
  ,_stmFragmentsTag(config().stmTag())
  ,_saveRawWithHeaderWaveform_HPGe(config().saveRawWithHeaderWaveform_HPGe())
  ,_saveRawWaveform_HPGe(config().saveRawWaveform_HPGe())
  ,_saveZSWaveform_HPGe(config().saveZSWaveform_HPGe())
  ,_saveRawWithHeaderWaveform_LaBr(config().saveRawWithHeaderWaveform_LaBr())
  ,_saveRawWaveform_LaBr(config().saveRawWaveform_LaBr())
  ,_saveZSWaveform_LaBr(config().saveZSWaveform_LaBr())
  ,_saveSTMFragSummary(config().saveSTMFragSummary())
  ,_verbosityLevel(config().verbosityLevel() ? *(config().verbosityLevel()) : 0)

{
  if (_saveSTMFragSummary) {
    produces<mu2e::STMFragmentSummaryCollection>("stmFragSummaryHPGe");
    produces<mu2e::STMFragmentSummaryCollection>("stmFragSummaryLaBr");
  }
  // Set the size of vectors for HPGe
  if (_saveRawWithHeaderWaveform_HPGe){produces<mu2e::STMWaveformDigiCollection>("rawWithHeaderHPGe");}
  if (_saveRawWaveform_HPGe){produces<mu2e::STMWaveformDigiCollection>("rawHPGe");}//Waveforms
  if (_saveZSWaveform_HPGe){produces<mu2e::STMWaveformDigiCollection>("zsHPGe");}
  produces<mu2e::STMPHDigiCollection>("phHPGe"); // digi series


  //Set the size of vectors for LaBr
  if (_saveRawWithHeaderWaveform_LaBr){produces<mu2e::STMWaveformDigiCollection>("rawWithHeaderLaBr");}
  if (_saveRawWaveform_LaBr){produces<mu2e::STMWaveformDigiCollection>("rawLaBr");}
  if (_saveZSWaveform_LaBr){produces<mu2e::STMWaveformDigiCollection>("zsLaBr");}
  produces<mu2e::STMPHDigiCollection>("phLaBr"); // digi series

}

// ----------------------------------------------------------------------

void STMDigisFromFragments::produce(Event& event)
{

  art::Ptr<mu2e::STMWaveformDigi> lastRawHPGePtr;
  art::Ptr<mu2e::STMWaveformDigi> lastRawLaBrPtr;

  bool haveLastRawHPGePtr{false};
  bool haveLastRawLaBrPtr{false};

  std::vector<int> zsToRawIndexHPGe; // creates vector
  std::vector<int> zsToRawIndexLaBr;
  ++_totalEvents; //Increment Total Event Counter

  //Frag Summaries
  std::unique_ptr<mu2e::STMFragmentSummaryCollection> stmFragSummaryHPGe(new mu2e::STMFragmentSummaryCollection);
  std::unique_ptr<mu2e::STMFragmentSummaryCollection> stmFragSummaryLaBr(new mu2e::STMFragmentSummaryCollection);
  //HPGe
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_HPGe_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> zs_HPGe_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMPHDigiCollection> ph_HPGe_digis(new mu2e::STMPHDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_HPGe_header_waveform_digis(new mu2e::STMWaveformDigiCollection);

  //LaBr
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_LaBr_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> zs_LaBr_waveform_digis(new mu2e::STMWaveformDigiCollection);
  std::unique_ptr<mu2e::STMPHDigiCollection> ph_LaBr_digis(new mu2e::STMPHDigiCollection);
  std::unique_ptr<mu2e::STMWaveformDigiCollection> raw_LaBr_header_waveform_digis(new mu2e::STMWaveformDigiCollection);

  art::Handle<artdaq::Fragments> STMFragmentsH;
  event.getByLabel(_stmFragmentsTag, STMFragmentsH);
  const auto& STMFragments = STMFragmentsH.product();

  auto rawHPGeProductID = event.getProductID<mu2e::STMWaveformDigiCollection>("rawHPGe");
  auto rawLaBrProductID = event.getProductID<mu2e::STMWaveformDigiCollection>("rawLaBr");

  //Event Metrics

  //HPGe
  size_t totalRawHPGeFrags{0};
  size_t totalZSHPGeFrags{0};
  size_t totalPHHPGeFrags{0};
  size_t goodRawHPGeFrags{0};
  size_t goodZSHPGeFrags{0};
  size_t goodPHHPGeFrags{0};
  size_t zeroRawHPGeFrags{0};
  size_t zeroZSHPGeFrags{0};
  size_t zeroPHHPGeFrags{0};
  size_t emptyRawHPGeFrags{0};
  size_t emptyZSHPGeFrags{0};
  size_t emptyPHHPGeFrags{0};

  //LaBr
  size_t totalRawLaBrFrags{0};
  size_t totalZSLaBrFrags{0};
  size_t totalPHLaBrFrags{0};
  size_t goodRawLaBrFrags{0};
  size_t goodZSLaBrFrags{0};
  size_t goodPHLaBrFrags{0};
  size_t zeroRawLaBrFrags{0};
  size_t zeroZSLaBrFrags{0};
  size_t zeroPHLaBrFrags{0};
  size_t emptyRawLaBrFrags{0};
  size_t emptyZSLaBrFrags{0};
  size_t emptyPHLaBrFrags{0};

  //Additional
  size_t unread_InnerFrags{0};
  uint16_t outerFragID{0};
  bool eventHasHPGe{false};
  bool eventHasLaBr{false};

  //HPGe
  uint16_t expectedZSLengthHPGe{0};
  uint16_t expectedZSRegionsHPGe{0};
  bool readZSinfoFromRawHeaderHPGe{false};

  //LaBr
  uint16_t expectedZSLengthLaBr{0};
  uint16_t expectedZSRegionsLaBr{0};
  bool readZSinfoFromRawHeaderLaBr{false};

  //Frag counters
  size_t nContainerFragsThisEvent{0};
  size_t nInnerFragsThisEvent{0};

  //loop over outer frags
  for (const auto& frag : *STMFragments) {
    ++_totalFragments; //Increment Total Frag counter

    outerFragID = frag.fragmentID();
    if (_verbosityLevel >= 3){std::cout << "\nFrag_id : " << outerFragID << "\n";}


    //Check if this is a container fragment
    if (frag.type() == artdaq::Fragment::ContainerFragmentType){

      mu2e::STMFragment container_frag(frag);
      if ( container_frag.isHPGeContainer()){ ++_totalContainersHPGe; eventHasHPGe = true; }
      else if ( container_frag.isLaBrContainer() ) { ++_totalContainersLaBr; eventHasLaBr = true; }
      else {
        if (_verbosityLevel >= 1){
          std::cout << "Encounter an unknown STM Container frag ID\n"
                    << "Frag ID : " << frag.fragmentID() <<"\n"
                    << "Event   : " << _totalEvents <<"\n";
        }
        continue;
      }
      artdaq::ContainerFragment cont_frag(frag);
      ++_totalContainers;
      ++nContainerFragsThisEvent;
      size_t blocks = cont_frag.block_count();
      _totalInner += blocks;
      nInnerFragsThisEvent += blocks;

      //loop over container where i corresponds to inner frag
      for (size_t i = 0; i < cont_frag.block_count(); ++i){

        auto inner_frag = cont_frag.at(i);
        mu2e::STMFragment stm_frag(*inner_frag);
        mu2e::STMWaveformDigi stm_waveform;

        if ( stm_frag.isRaw() ) {

          if (stm_frag.isHPGe()){
            readZSinfoFromRawHeaderHPGe = false;
            expectedZSLengthHPGe = 0;
            expectedZSRegionsHPGe = 0;
          } else if (stm_frag.isLaBr()){
            readZSinfoFromRawHeaderLaBr = false;
            expectedZSLengthLaBr = 0;
            expectedZSRegionsLaBr = 0;
          }
          //Job Counter
          ++_totalRaw;

          //Conditional Job and Event Counter
          if( stm_frag.isHPGe() ){
            ++_totalRawHPGe ;
            ++totalRawHPGeFrags;

          } else if ( stm_frag.isLaBr() ){
            ++_totalRawLaBr ;
            ++totalRawLaBrFrags;
          }

          auto payloadPtr = stm_frag.payloadBegin();
          auto payloadWords = stm_frag.payloadWords();
          bool allZeros = true;//Assume raw frag is zero filled

          //Checks for empty frag
          if (payloadWords == 0) {
            if(_verbosityLevel >=3){std::cout<< "\nFound an empty frag, i = " << i <<" @Raw\n";}
            ++_totalEmptyRaw;//Job counter
            //Eventcounter
            if( stm_frag.isHPGe() ){
              ++_totalEmptyRawHPGe; ++emptyRawHPGeFrags; }
            else if (stm_frag.isLaBr()){
              ++_totalEmptyRawLaBr; ++emptyRawLaBrFrags; }
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
            if( stm_frag.isHPGe() ){
              ++_totalZeroRawHPGe; ++zeroRawHPGeFrags; }
            else if ( stm_frag.isLaBr() ){
              ++_totalZeroRawLaBr; ++zeroRawLaBrFrags; }
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

          //Ideally only good frags get up to here
          if ( stm_frag.isHPGe() ){
            expectedZSRegionsHPGe = stm_frag.zsRegions();
            expectedZSLengthHPGe = stm_frag.zsLength();
            readZSinfoFromRawHeaderHPGe = true;

            ++ _totalGoodRawHPGe;
            ++ goodRawHPGeFrags;

            if (_saveRawWithHeaderWaveform_HPGe){
              auto dataPtr = stm_frag.dataBegin();//memory check for emplace_back
              auto dataWords = stm_frag.dataWords();
              stm_waveform.set_data(dataWords, dataPtr);
              raw_HPGe_header_waveform_digis->emplace_back(stm_waveform);
            }
            if(_saveRawWaveform_HPGe) {
              stm_waveform.set_data(payloadWords, payloadPtr);
              raw_HPGe_waveform_digis->emplace_back(stm_waveform);
              size_t parentRawWfmIdx = raw_HPGe_waveform_digis->size()-1;

              lastRawHPGePtr = art::Ptr<mu2e::STMWaveformDigi>(rawHPGeProductID,parentRawWfmIdx,
                                                         event.productGetter(rawHPGeProductID));
              haveLastRawHPGePtr = true;
            }

          } else if ( stm_frag.isLaBr() ){
            expectedZSRegionsLaBr = stm_frag.zsRegions();
            expectedZSLengthLaBr = stm_frag.zsLength();
            readZSinfoFromRawHeaderLaBr = true;

            ++ _totalGoodRawLaBr;
            ++ goodRawLaBrFrags;

            if (_saveRawWithHeaderWaveform_LaBr){
              auto dataPtr = stm_frag.dataBegin();
              auto dataWords = stm_frag.dataWords();
              stm_waveform.set_data(dataWords, dataPtr);
              raw_LaBr_header_waveform_digis->emplace_back(stm_waveform);
            }
            if(_saveRawWaveform_LaBr){
              stm_waveform.set_data(payloadWords, payloadPtr);
              raw_LaBr_waveform_digis->emplace_back(stm_waveform);
              size_t parentRawWfmIdx = raw_LaBr_waveform_digis->size()-1;
              lastRawLaBrPtr = art::Ptr<mu2e::STMWaveformDigi>(rawLaBrProductID,parentRawWfmIdx,
                                                         event.productGetter(rawLaBrProductID));
              haveLastRawLaBrPtr = true;

            }

          }

        } // End of isRaw Block


        else if (stm_frag.isZS()){
          ++_totalZS; //Incremenet ZS counter

          if ( stm_frag.isHPGe() ){
            ++_totalZSHPGe; ++totalZSHPGeFrags; }
          else if ( stm_frag.isLaBr() ){ ++_totalZSLaBr; ++totalZSLaBrFrags;
          }

          auto payloadPtr = stm_frag.payloadBegin();
          auto payloadWords = stm_frag.payloadWords();
          bool allZeros = true; //assumes all adcs are zero

          // Extract variables from Raw Header
          bool readZSinfoFromRawHeader = stm_frag.isHPGe() ? readZSinfoFromRawHeaderHPGe : readZSinfoFromRawHeaderLaBr;
          uint16_t expectedZSRegions = stm_frag.isHPGe() ? expectedZSRegionsHPGe : expectedZSRegionsLaBr;
          uint16_t expectedZSLength = stm_frag.isHPGe() ? expectedZSLengthHPGe : expectedZSLengthLaBr;

          //After copying varibales, reset detector specific variables
          if (stm_frag.isHPGe()){
            readZSinfoFromRawHeaderHPGe = false;
            expectedZSLengthHPGe = 0;
            expectedZSRegionsHPGe = 0;
          } else if (stm_frag.isLaBr()){
            readZSinfoFromRawHeaderLaBr = false;
            expectedZSLengthLaBr = 0;
            expectedZSRegionsLaBr = 0;
          }

          //Check if payload is empty
          if ( payloadWords == 0) {
            if (_verbosityLevel >=3){std::cout << "\nFound an empty frag, i = " << i << " @ZS\n";}
            ++_totalEmptyZS;
            if ( stm_frag.isHPGe() ){
              ++_totalEmptyZSHPGe; ++emptyZSHPGeFrags; }
            else if ( stm_frag.isLaBr() ){
              ++_totalEmptyZSLaBr; ++emptyZSLaBrFrags; }
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
            if ( stm_frag.isHPGe() ){
              ++_totalZeroZSHPGe; ++zeroZSHPGeFrags; }
            else if ( stm_frag.isLaBr() ){
              ++_totalZeroZSLaBr; ++zeroZSLaBrFrags; }
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
            if ( readZSinfoFromRawHeader && seg >= expectedZSRegions ) break;
            uint16_t current_zs_location = static_cast<uint16_t>(dataPtr[0]);
            uint16_t current_zs_size = static_cast<uint16_t>(dataPtr[1]);
            auto adc = dataPtr + 2;
            if (adc + current_zs_size > dataEnd)
              break;

            uint32_t trigTimeOffset = current_zs_location;

            //emplacing
            if ( stm_frag.isHPGe() && _saveZSWaveform_HPGe){
              std::vector<int16_t> segADCS(adc, adc + current_zs_size); //1D array, contains adcs to this_zs_size - 1
              mu2e::STMWaveformDigi zsDigi(trigTimeOffset, segADCS);

              if (readZSinfoFromRawHeader && haveLastRawHPGePtr) {
                zsDigi.setParent(lastRawHPGePtr);
              }

              zs_HPGe_waveform_digis->emplace_back(zsDigi);
              if(_verbosityLevel >= 3 && zsDigi.hasParent()){
                std::cout << "ZS HPGe segment has parent Raw index "
                          << zsDigi.parent().key()
                          << "with offset "
                          << zsDigi.trigTimeOffset()
                          << "\n";
              }
            }
            else if ( stm_frag.isLaBr() && _saveZSWaveform_LaBr){
              std::vector<int16_t> segADCS(adc, adc + current_zs_size);
              mu2e::STMWaveformDigi zsDigi(trigTimeOffset, segADCS);

              if (readZSinfoFromRawHeader && haveLastRawLaBrPtr){
                zsDigi.setParent(lastRawLaBrPtr);
              }

              zs_LaBr_waveform_digis->emplace_back(zsDigi);

              }

            if (_verbosityLevel >=6){

              //A print check per segment
              std::cout << "Region = " << seg << " , zs_index = "
                        << current_zs_location << " , zs_size = " << current_zs_size
                        << " , trigTimeOffset = " << trigTimeOffset << "\n" ;
            }

            //Update variables
            lastZSindex = current_zs_location;
            lastLen = current_zs_size;
            totalLen += lastLen;
            ++seg;
            dataPtr = adc + current_zs_size;
          } // end of while loop

          if (_verbosityLevel >=5){
            //Summary
            std::cout << "ZS Regions = " << seg
                      << " , lastZSindex = " << lastZSindex
                      << " , lastZSLen = " << lastLen
                      << " , ZS total length = " << totalLen
                      << " , i =  " << i <<  " @ZS\n";
          }

          //Throw out if ZSLengthfromRaw != totalLen, throw out if length mismatch
          if (readZSinfoFromRawHeader){
            if (expectedZSLength != totalLen){
              throw cet::exception("STM_UNPACKING")
                << "\n=== ZS Length count mismatch ===\n"
                << "ZS length from Raw header : " << expectedZSLength << "\n"
                << "ZS length calculated from file : " << totalLen << "\n"
                << "Found at inner frag i : " << i << "\n"
                << "Encountered at event : " << _totalEvents << "\n"      ;
            }//Also throw out if regions for ZS pulse do not match
            if ( expectedZSRegions != seg){
              throw cet::exception("STM_UNPACKING")
                << "\n=== ZS Region count mismatch ===\n"
                << "ZS Region count from Raw header : " << expectedZSRegions << "\n"
                << "ZS Regionc count calculated from file : " << seg << "\n"
                << "Found at inner frag i : " << i << "\n"
                << "Encountered at event : " << _totalEvents << "\n"      ;
            }
          }

          //Increment counter after checking if there is a match
          if ( stm_frag.isHPGe()){ ++_totalGoodZSHPGe; ++goodZSHPGeFrags; }
          else if(stm_frag.isLaBr()){++ _totalGoodZSLaBr; ++goodZSLaBrFrags; }

        }//End of isZS

        else if ( stm_frag.isPH() ){

          ++_totalPH;
          if ( stm_frag.isHPGe() ){ ++_totalPHHPGe; ++totalPHHPGeFrags;}
          else if ( stm_frag.isLaBr() ){ ++_totalPHLaBr; ++totalPHLaBrFrags; }

          //Check if zero filled
          auto payloadPtr = stm_frag.payloadBegin();
          auto payloadWords = stm_frag.payloadWords();
          bool allZeros = true;


          if (payloadWords ==0){
            if (_verbosityLevel >= 3){std::cout << "\nFound an empty frag, i = " << i<< " @PH\n";}
            ++_totalEmptyPH;
            if ( stm_frag.isHPGe() ){ ++_totalEmptyPHHPGe; ++emptyPHHPGeFrags; }
            else if ( stm_frag.isLaBr() ){ ++_totalEmptyPHLaBr; ++emptyPHLaBrFrags;}
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

            if ( stm_frag.isHPGe() ){ ++_totalZeroPHHPGe; ++zeroPHHPGeFrags; }
            else if ( stm_frag.isLaBr() ){ ++_totalZeroPHLaBr; ++zeroPHLaBrFrags; }
            continue;
          }

          if (_verbosityLevel >= 3){std::cout << "\nFound a good frag, i = " << i <<" @PH\n";}

          if ( stm_frag.isHPGe() ){ ++ _totalGoodPHHPGe; ++goodPHHPGeFrags; }
          else if ( stm_frag.isLaBr() ) { ++ _totalGoodPHLaBr; ++goodPHLaBrFrags; }

          size_t digiWords = stm_frag.payloadWords();
          auto const* digiPtr = stm_frag.payloadBegin();

          for (size_t i_PH = 0; i_PH < digiWords ; ++i_PH){
            int16_t PH = digiPtr[i_PH];
            mu2e::STMPHDigi PH_digi(0, PH);

            if( stm_frag.isHPGe() ){
              ph_HPGe_digis->emplace_back(PH_digi); }
            else if ( stm_frag.isLaBr() ){
              ph_LaBr_digis->emplace_back(PH_digi); }

          }

        }//End of isPH and is checks
        else{
          ++ _unreadInnerFrags; //For Job summary
          ++ unread_InnerFrags; //For event summary
          if(_verbosityLevel >=3){
            std::cout << "Unread Inner fragment" << "\n"
                      << "Inner frag i : " << i << "\n"
                      << "Frag ID      : " << inner_frag->fragmentID() << "\n"
                      << "Event        : " <<_totalEvents << "\n";
          }
        }
      }

    } else {
      //fallback (non-container case)
      if (_verbosityLevel >=1){
        std::cout << "Found non-container STM Fragment " <<"\n"
                  << "Fragment ID : " << frag.fragmentID() << "\n"
                  << "Event : " << _totalEvents << "\n";
      }
      ++ _totalNonContainerFrags;
      continue;
    }
  } //End of frags loop

  if (eventHasHPGe && eventHasLaBr){ ++_totalEventsWithHPGeLaBr; }
  else if (eventHasHPGe){ ++_totalEventsWithOnlyHPGe; }
  else if (eventHasLaBr){ ++_totalEventsWithOnlyLaBr; }
  else { ++_totalEventsWithNone; }

  if(_saveSTMFragSummary){
    stmFragSummaryHPGe->emplace_back(nContainerFragsThisEvent,nInnerFragsThisEvent,
                                     goodRawHPGeFrags, goodZSHPGeFrags, goodPHHPGeFrags,
                                     zeroRawHPGeFrags, zeroZSHPGeFrags, zeroPHHPGeFrags,
                                     emptyRawHPGeFrags, emptyZSHPGeFrags, emptyPHHPGeFrags);

    stmFragSummaryLaBr->emplace_back(nContainerFragsThisEvent,nInnerFragsThisEvent,
                                     goodRawLaBrFrags, goodZSLaBrFrags, goodPHLaBrFrags,
                                     zeroRawLaBrFrags, zeroZSLaBrFrags, zeroPHLaBrFrags,
                                     emptyRawLaBrFrags, emptyZSLaBrFrags, emptyPHLaBrFrags);

  }
  if (_verbosityLevel >= 2 ){
    //Event Summary -> tells us what happens per art event
    std::cout << "\n========== STM EVENT SUMMARY - (Unpacking Module) ==========\n";
    std::cout << "\n--- Products extracted ---\n";
    std::cout << "Extracted Raw HPGe waveforms     : " << raw_HPGe_waveform_digis->size() <<"\n";
    std::cout << "Extracted ZS  HPGe waveforms     : " << zs_HPGe_waveform_digis->size() <<"\n";
    std::cout << "Extracted PH  HPGe digis         : " << ph_HPGe_digis->size() <<"\n";
    std::cout << "\n";
    std::cout << "Extracted Raw LaBr waveforms     : " << raw_LaBr_waveform_digis->size() <<"\n";
    std::cout << "Extracted ZS  LaBr waveforms     : " << zs_LaBr_waveform_digis->size() <<"\n";
    std::cout << "Extracted PH  LaBr digis         : " << ph_LaBr_digis->size() <<"\n";


    std::cout << "\n--- Frags Seen ---\n"; // Counts all HPGe and LaBr frags, no fiklter
    std::cout << "Raw HPGe frags seen    : " << totalRawHPGeFrags << "\n";
    std::cout << "ZS  HPGe frags seen    : " << totalZSHPGeFrags << "\n";
    std::cout << "PH  HPGe frags seen    : " << totalPHHPGeFrags << "\n";
    std::cout << "\n";
    std::cout << "Raw LaBr frags seen    : " << totalRawLaBrFrags << "\n";
    std::cout << "ZS  LaBr frags seen    : " << totalZSLaBrFrags << "\n";
    std::cout << "PH  LaBr frags seen    : " << totalPHLaBrFrags << "\n";
    std::cout << "\nUnread frags           : " << unread_InnerFrags << "\n";


    std::cout << "\n--- Filtered results ---\n";
    std::cout << "Good Raw  HPGe frags   : " << goodRawHPGeFrags << "\n";
    std::cout << "Good ZS   HPGe frags   : " << goodZSHPGeFrags << "\n";
    std::cout << "Good PH   HPGe frags   : " << goodPHHPGeFrags << "\n";
    std::cout << "\n";
    std::cout << "Good Raw  LaBr frags   : " << goodRawLaBrFrags << "\n";
    std::cout << "Good ZS   LaBr frags   : " << goodZSLaBrFrags << "\n";
    std::cout << "Good PH   LaBr frags   : " << goodPHLaBrFrags << "\n";
    std::cout << "\n";
    std::cout << "Zero Raw  HPGe frags   : " << zeroRawHPGeFrags << "\n";
    std::cout << "Zero ZS   HPGe frags   : " << zeroZSHPGeFrags << "\n";
    std::cout << "Zero PH   HPGe frags   : " << zeroPHHPGeFrags << "\n";
    std::cout << "\n";
    std::cout << "Zero Raw  LaBr frags   : " << zeroRawLaBrFrags << "\n";
    std::cout << "Zero ZS   LaBr frags   : " << zeroZSLaBrFrags << "\n";
    std::cout << "Zero PH   LaBr frags   : " << zeroPHLaBrFrags << "\n";
    std::cout << "\n";
    std::cout << "Empty Raw HPGe frags   : " << emptyRawHPGeFrags <<"\n";
    std::cout << "Empty ZS  HPGe frags   : " << emptyZSHPGeFrags << "\n";
    std::cout << "Empty PH  HPGe frags   : " << emptyPHHPGeFrags << "\n";
    std::cout << "\n";
    std::cout << "Empty Raw LaBr frags   : " << emptyRawLaBrFrags <<"\n";
    std::cout << "Empty ZS  LaBr frags   : " << emptyZSLaBrFrags << "\n";
    std::cout << "Empty PH  LaBr frags   : " << emptyPHLaBrFrags << "\n";

    std::cout << "=================================\n";
  }

  //Final move
  if (_saveSTMFragSummary) {
    event.put(std::move(stmFragSummaryHPGe),"stmFragSummaryHPGe");
    event.put(std::move(stmFragSummaryLaBr),"stmFragSummaryLaBr");
  }

  //HPGe
  if (_saveRawWithHeaderWaveform_HPGe){ event.put(std::move(raw_HPGe_header_waveform_digis), "rawWithHeaderHPGe"); }
  if (_saveRawWaveform_HPGe){event.put(std::move(raw_HPGe_waveform_digis), "rawHPGe");}
  if (_saveZSWaveform_HPGe){event.put(std::move(zs_HPGe_waveform_digis), "zsHPGe");}
  event.put(std::move(ph_HPGe_digis), "phHPGe");

  //LaBr
  if (_saveRawWithHeaderWaveform_LaBr){ event.put(std::move(raw_LaBr_header_waveform_digis), "rawWithHeaderLaBr"); }
  if (_saveRawWaveform_LaBr){event.put(std::move(raw_LaBr_waveform_digis), "rawLaBr");}
  if (_saveZSWaveform_LaBr){event.put(std::move(zs_LaBr_waveform_digis), "zsLaBr");}
  event.put(std::move(ph_LaBr_digis), "phLaBr");

} // produce()

// ======================================================================


void STMDigisFromFragments::endJob() {
  if (_verbosityLevel >= 1){
    //Art job summary for the unpacking
    std::cout << "\n========== STM JOB SUMMARY - (Unpacking Module) ==========\n";

    std::cout << "Total Art events         : " << _totalEvents << "\n";
    std::cout << "Total Art events w/HPGe&LaBr : " << _totalEventsWithHPGeLaBr << "\n";
    std::cout << "Total Art events w/only HPGe : " << _totalEventsWithOnlyHPGe << "\n";
    std::cout << "Total Art events w/only LaBr : " << _totalEventsWithOnlyLaBr << "\n";
    std::cout << "Total Art events w/Neither   : " << _totalEventsWithNone << "\n";
    std::cout << "Total frags              : " << _totalFragments << "\n";
    std::cout << "Total Container frags    : " << _totalContainers << "\n";
    std::cout << "Total HPGe Containers    : " << _totalContainersHPGe << "\n";
    std::cout << "Total LaBr Containers    : " << _totalContainersLaBr << "\n";
    std::cout << "Total Inner frags        : " << _totalInner << "\n";
    std::cout << "Unread frags             : " << _unreadInnerFrags << "\n";
    std::cout << "Non Container frags      : " << _totalNonContainerFrags << "\n";

    std::cout << "\n--- Data types read pre filtering ---\n";
    std::cout << "Total RAW frags seen              : " << _totalRaw << "\n";
    std::cout << "Total ZS  frags seen              : " << _totalZS << "\n";
    std::cout << "total PH  frags seen              : " << _totalPH << "\n";
    std::cout << "\n";
    std::cout << "Total RAW HPGe frags seen         : " << _totalRawHPGe << "\n";
    std::cout << "Total ZS  HPGe frags seen         : " << _totalZSHPGe << "\n";
    std::cout << "Total PH  HPGe frags seen         : " << _totalPHHPGe << "\n";
    std::cout << "\n";
    std::cout << "Total RAW LaBr frags seen         : " << _totalRawLaBr << "\n";
    std::cout << "Total ZS  LaBr frags seen         : " << _totalZSLaBr << "\n";
    std::cout << "Total PH  LaBr frags seen         : " << _totalPHLaBr << "\n";

    std::cout << "\n--- Data types filtered HPGe ---\n";
    std::cout << "Good RAW   HPGe frags    : " << _totalGoodRawHPGe << "\n";
    std::cout << "Good ZS    HPGe frags    : " << _totalGoodZSHPGe << "\n";
    std::cout << "Good PH    HPGe frags    : " << _totalGoodPHHPGe << "\n";
    std::cout << "\n";
    std::cout << "Zero RAW   HPGe frags    : " << _totalZeroRawHPGe << "\n";
    std::cout << "Zero ZS    HPGe frags    : " << _totalZeroZSHPGe << "\n";
    std::cout << "Zero PH    HPGe frags    : " << _totalZeroPHHPGe << "\n";
    std::cout << "\n";
    std::cout << "Empty Raw  HPGe frags    : " << _totalEmptyRawHPGe << "\n";
    std::cout << "Empty ZS   HPGe frags    : " << _totalEmptyZSHPGe << "\n";
    std::cout << "Empty PH   HPGe frags    : " << _totalEmptyPHHPGe << "\n";

    std::cout << "\n--- Data types filtered LaBr ---\n";
    std::cout << "Good RAW   LaBr frags    : " << _totalGoodRawLaBr << "\n";
    std::cout << "Good ZS    LaBr frags    : " << _totalGoodZSLaBr << "\n";
    std::cout << "Good PH    LaBr frags    : " << _totalGoodPHLaBr << "\n";
    std::cout << "\n";
    std::cout << "Zero RAW   LaBr frags    : " << _totalZeroRawLaBr << "\n";
    std::cout << "Zero ZS    LaBr frags    : " << _totalZeroZSLaBr << "\n";
    std::cout << "Zero PH    LaBr frags    : " << _totalZeroPHLaBr << "\n";
    std::cout << "\n";
    std::cout << "Empty Raw  LaBr frags    : " << _totalEmptyRawLaBr << "\n";
    std::cout << "Empty ZS   LaBr frags    : " << _totalEmptyZSLaBr << "\n";
    std::cout << "Empty PH   LaBr frags    : " << _totalEmptyPHLaBr << "\n";

    std::cout << "=================================\n";

  }
}

DEFINE_ART_MODULE(STMDigisFromFragments)
