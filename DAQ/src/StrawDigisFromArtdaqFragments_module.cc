// ======================================================================
// PM
// StrawDigisFromArtdaqFragments:  add tracker data products to the event
// debugMode_ > 0 : enables diagnostic printouts
// debugBit_[0]: ???
// debugBit_[1]: digis
// debugBit_[2]: waveforms
// ======================================================================
#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_RocDataHeaderPacket.h"

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"

#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>

#include <iostream>

#include <string>

#include <map>
#include <memory>

// #define TRACEMF_USE_VERBATIM 1

// #include "TRACE/tracemf.h"
// #define TRACE_NAME "StrawDigisFromArtdaqFragments"

#include "Offline/DAQ/inc/TrkPanelMap_t.hh"

namespace mu2e {
  class StrawDigisFromArtdaqFragments;
}
// ======================================================================

class mu2e::StrawDigisFromArtdaqFragments : public art::EDProducer {

public:

  struct Config {
    fhicl::Atom<int>              diagLevel    {fhicl::Name("diagLevel"    ), fhicl::Comment("diagnostic severity level, default = 0" ), 0};
    fhicl::Atom<int>              debugMode    {fhicl::Name("debugMode"    ), fhicl::Comment("debug mode, default = 0"                ), 0};
    fhicl::Sequence<std::string>  debugBits    {fhicl::Name("debugBits"    ), fhicl::Comment("debug bits"                             )   };
    fhicl::Atom<bool>             saveWaveforms{fhicl::Name("saveWaveforms"), fhicl::Comment("save StrawDigiADCWaveforms, default = true"), true};

    // individual tuple specifying a minnesota label, e.g. MN123,
    // with geographic plane/panel numbers, i.e. from DocDB-#888
    struct GeographicTuple{
      fhicl::Atom<std::string> minnesota{
        fhicl::Name("minnesota"),
        fhicl::Comment("Minnesota # label")
      };
      fhicl::Atom<uint16_t> plane{
        fhicl::Name("plane"),
        fhicl::Comment("Geographic plane [DocDB-888]")
      };
      fhicl::Atom<uint16_t> panel{
        fhicl::Name("panel"),
        fhicl::Comment("Geographic panel [DocDB-888]")
      };
    };
    // mandatory listing of geographic entries, to enable translation
    // of panel-labelings in the data to Offline StrawIds
    fhicl::Sequence< fhicl::Table<GeographicTuple> > geography{
      fhicl::Name("geography"),
      fhicl::Comment("Mapping of Minnesota numbers to geographic planes and panels")
    };

    // individual tuple specifying a minnesota label, e.g. MN123,
    // with logical channeling, i.e. DTC ID and associated Link #
    struct LogicalTuple{
      fhicl::Atom<uint16_t> dtc{
        fhicl::Name("dtc"),
        fhicl::Comment("DTC ID")
      };
      fhicl::Atom<uint16_t> link{
        fhicl::Name("link"),
        fhicl::Comment("Link #")
      };
      fhicl::Atom<std::string> minnesota{
        fhicl::Name("minnesota"),
        fhicl::Comment("Minnesota # label")
      };
    };
    // optional listing of logical entries, to provide a backup
    // translation in case of invalid labeling in the data
    fhicl::OptionalSequence< fhicl::Table<LogicalTuple> > channeling{
      fhicl::Name("channeling"),
      fhicl::Comment("Logical channeling of panels (optional)")
    };
  };

  // --- C'tor/d'tor:
  explicit StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config);
  virtual ~StrawDigisFromArtdaqFragments() {}

  void         print_(const std::string&  Message, int DiagLevel = -1,
                      const std::source_location& location = std::source_location::current());

  void         print_fragment(const artdaq::Fragment* Frag);

    // --- overloaded functions of the art producer
  virtual void produce (art::Event& ArtEvent) override;
  virtual void beginRun(art::Run&   ArtRun  ) override;

  enum {kNDebugBits = 100};

private:
                                        // talk-to parameters
  int       diagLevel_    ;
  int       debugMode_   ;
  
  std::vector<std::string> debugBits_;
  int                      debugBit_[kNDebugBits];  
  bool      saveWaveforms_;
                                        // the rest
  int       nADCPackets_{-1};           // N(ADC packets per hit)
  int       nSamples_   {-1};           // N(ADC samples per hit)
  int       np_per_hit_ {-1};           // N(data packets per hits)

  std::map<uint16_t, uint16_t> minnesota_map_; // mapping from minnesota number to upper bits of StrawId
  uint16_t channel_map_[36][6] ; // mapping from DTC link number to minnesota number

  const art::Event*        event_;
                                                // for now, IDTC=2*nodename+PCIE_ADDR
  const TrkPanelMap_t*     panel_map_[50][6];   // panel_map_[idtc][ilink] = pointer to future DB record
                                                // make sure DAQ22 fits in

  // less than 300 panels physically exist and are enumeratively labeled
  // hence, the max allowed word can act be used as a sentinel
  const static uint16_t invalid_minnesota_ = static_cast<uint16_t>(-1);
  uint16_t parse_minnesota_label(std::string label);
};

// ======================================================================
mu2e::StrawDigisFromArtdaqFragments::StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_    (config().diagLevel    ()),
    debugMode_    (config().debugMode    ()),
    debugBits_    (config().debugBits    ()),
    saveWaveforms_(config().saveWaveforms()),
    event_(nullptr)
{
  produces<mu2e::StrawDigiCollection>();
  if (saveWaveforms_) produces<mu2e::StrawDigiADCWaveformCollection>();

  produces<mu2e::IntensityInfoTrackerHits>();

  // initialize geographic mapping of minnesota-labled panels
  for (const auto& entry: config().geography()){
    const auto& minnesota = entry.minnesota();
    const auto mid = parse_minnesota_label(minnesota);
    auto plane = entry.plane();
    auto panel = entry.panel();
    mu2e::StrawId pid(plane, panel, 0);
    if (0 < minnesota_map_.count(mid)){
      std::string msg = "duplicate mapping of panel " + minnesota;
      throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
    }
    minnesota_map_[mid] = pid.getPanelId().asUint16();
  }

  // initialize fallback mapping of dtc links to minnesota-labeling
  for (size_t i = 0 ; i < mu2e::StrawId::_nplanes ; i++){
    for (size_t j = 0 ; j < mu2e::StrawId::_npanels ; j++){
      channel_map_[i][j] = StrawDigisFromArtdaqFragments::invalid_minnesota_;
    }
  }
  const auto channeling = config().channeling();
  if (channeling.has_value()){
    for (const auto& entry: channeling.value()){
      uint16_t dtc = entry.dtc();
      if (!(dtc < mu2e::StrawId::_nplanes)){
        std::string msg = "invalid DTC ID: " + std::to_string(dtc);
        throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
      }
      uint16_t link = entry.link();
      if (!(link < mu2e::StrawId::_npanels)){
        std::string msg = "invalid DTC Link number: " + std::to_string(link);
        throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
      }
      std::string minnesota = entry.minnesota();
      uint16_t mid = parse_minnesota_label(minnesota);
      if (minnesota_map_.count(mid) < 1){
        std::string msg = "dtc link mapping defined for unmapped panel " + minnesota;
        throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
      }
      channel_map_[dtc][link] = mid;
    }
  }
//-----------------------------------------------------------------------------
// initialize debug bits  : debugBits: [ "bit0:1" , "bit14:1" ]
//-----------------------------------------------------------------------------
  for (int i=0; i<kNDebugBits; ++i) debugBit_[i] = 0;

  const char* key;
  int nbits = debugBits_.size(); // from FCL

  for (int i=0; i<nbits; i++) {
    int index(0), value(0);
    key               = debugBits_[i].data();
    sscanf(key,"bit%i:%i",&index,&value);
    debugBit_[index]  = value;
        
    print_(std::format("StrawDigisFromArtdaqFragments: bit={:4d} is set to {}\n",index,debugBit_[index]));
  }
}


//-----------------------------------------------------------------------------
void mu2e::StrawDigisFromArtdaqFragments::print_(const std::string& Message, int DiagLevel,
                                                const std::source_location& location) {
  if (DiagLevel > diagLevel_) return;
  if (event_) {
    std::cout << std::format(" event:{}:{}:{} ",event_->run(),event_->subRun(),event_->event());
  }
  std::cout << location.file_name() << ":" << location.line()
    //            << location.function_name()
            << ": " << Message;
}

//-----------------------------------------------------------------------------
// HEX print of a fragment, the data has to be in 2-byte words
//-----------------------------------------------------------------------------
void mu2e::StrawDigisFromArtdaqFragments::print_fragment(const artdaq::Fragment* Frag) {
  ushort* buf = (ushort*) (Frag->dataBegin());
  int nw      = buf[0]/2;
  int loc     = 0;

  for (int i=0; i<nw; i++) {
    if (loc == 0) printf(" 0x%08x: ",i*2);

    ushort  word = buf[i];
    printf("0x%04x ",word);

    loc += 1;
    if (loc == 8) {
      printf("\n");
      loc = 0;
    }
  }

  if (loc != 0) printf("\n");
}

//-----------------------------------------------------------------------------
void mu2e::StrawDigisFromArtdaqFragments::beginRun(art::Run&  ArtRun) {
  // fill panel_map_ for a given run - should come from the database
  
  for (const TrkPanelMap_t* tpm = TrkPanelMap_data.begin(); tpm != TrkPanelMap_data.end(); ++tpm) {
    int dtc               = tpm->dtc;
    int link              = tpm->link;
    panel_map_[dtc][link] = tpm;
  }
}
// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void mu2e::StrawDigisFromArtdaqFragments::produce(art::Event& event) {
  int const packet_size(16); // in bytes

  event_ = &event;                      // cache for printouts
  
  if (debugMode_ > 0) print_("-- START\n",1);

   // Collection of StrawDigis for the event
  std::unique_ptr<mu2e::StrawDigiCollection> straw_digis(new mu2e::StrawDigiCollection);
  std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> straw_digi_adcs(new mu2e::StrawDigiADCWaveformCollection);

  // IntensityInfoTrackerHits
  std::unique_ptr<mu2e::IntensityInfoTrackerHits> intInfo(new mu2e::IntensityInfoTrackerHits);

//-----------------------------------------------------------------------------
// defined by the first hit
//-----------------------------------------------------------------------------
  artdaq::Fragments    fragments;
  artdaq::FragmentPtrs containerFragments;

  auto fragmentHandles = event.getMany<std::vector<artdaq::Fragment>>();

  for (auto handle : fragmentHandles) {
    if (!handle.isValid() || handle->empty())     continue;

    if (handle->front().type() == artdaq::Fragment::ContainerFragmentType) {
      for (const auto& cont : *handle) {
        artdaq::ContainerFragment contf(cont);
        for (size_t ii = 0; ii < contf.block_count(); ++ii) {
          containerFragments.push_back(contf[ii]);
          fragments.push_back(*containerFragments.back());
        }
      }
    }
    else {
//-----------------------------------------------------------------------------
// the 'handle' handles a list of artdaq fragments
// each artdaq fragment corresponds to a single DTC, or a plane
// loop over them
//-----------------------------------------------------------------------------
      int nfrag = handle->size();
      for (int ifrag=0; ifrag<nfrag; ifrag++) {
        const artdaq::Fragment* frag = &handle->at(ifrag);
        uint8_t* fdata   = (uint8_t*) (frag->dataBegin());
        if (debugMode_ and (debugBit_[0] > 0)) {
          print_(std::format("-- fragment number:{}\n",ifrag));
          print_fragment(frag);
        }
//-----------------------------------------------------------------------------
// skip non-tracker fragments
// after a recent format change, a DTC fragment may contain ROC data from different
// subdetectors, make sure that at least one of them is the tracker ROC
//-----------------------------------------------------------------------------
        DTCLib::DTC_SubEventHeader* seh = (DTCLib::DTC_SubEventHeader*) fdata;
        if ((seh->link0_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link1_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link2_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link3_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link4_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) and
            (seh->link5_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker)     )
                                                            continue;
        int dtc_id = seh->source_dtc_id;
//-----------------------------------------------------------------------------
// this is a tracker DTC fragment, loop over the ROCs
//-----------------------------------------------------------------------------
        ushort*  buf          = (ushort*) fdata;
        int      nbytes       = buf[0];             // frag.dataSizeBytes() includes extra 0x20
        uint8_t* roc_data     = fdata+sizeof(*seh);
        uint8_t* last_address = fdata+nbytes;

        while (roc_data < last_address) {
          RocDataHeaderPacket_t* rdh = (RocDataHeaderPacket_t*) roc_data;
          int nhits = 0;
          int header_printed = 0;
//------------------------------------------------------------------------------
// skip empty ROC blocks
//------------------------------------------------------------------------------
          if (rdh->packetCount > 1) {
            if (nADCPackets_ < 0) {
//-----------------------------------------------------------------------------
// take the number of ADC packets per hit from the hit data,
// trust that but watch if it changes
// so far, any corruptions we saw were contained withing the ROC payload, and nhits
// was a reliable number
//-----------------------------------------------------------------------------
              mu2e::TrackerDataDecoder::TrackerDataPacket* h0;
              h0           = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+packet_size);
              nADCPackets_ = h0->NumADCPackets;
              nSamples_    = 3+12*nADCPackets_;
              np_per_hit_  = nADCPackets_+1;
            }
            int link_id = rdh->linkID;
            nhits       = rdh->packetCount/(nADCPackets_+1);

            if (debugMode_) {
              print_(std::format("--- DTC:{} ROC:{} nhits:{}\n",dtc_id,link_id,nhits));
            }
            for (int ihit=0; ihit<nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
              mu2e::TrackerDataDecoder::TrackerDataPacket* hit_data ;
              
              int offset = (ihit*np_per_hit_+1)*packet_size;   // in bytes
              hit_data   = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+offset);
//-----------------------------------------------------------------------------
// at this point, check consistency between the channel_id, dtc_id and link_id for a given run
// panel ID is a derivative of the DTC ID and the link iD
// mn_id - 'MinnesotaID' of the panel
//-----------------------------------------------------------------------------
              mu2e::StrawDigiFlag digi_flag;
              uint16_t channel = static_cast<uint16_t>(hit_data->StrawIndex);
              uint16_t chid   = mu2e::StrawId(channel).straw(); // channel ID within the panel
              uint16_t mnid    = channel >> mu2e::StrawId::_panelsft;
              
              const TrkPanelMap_t* pm = panel_map_[dtc_id][link_id];  // DB here
              if (pm->mnid != mnid) {
                print_(std::format("ERROR: hit chid:{:04x} inconsistent with the dtc_id:{} and link_id:{}\n",
                                   hit_data->StrawIndex, dtc_id, link_id));
//-----------------------------------------------------------------------------
// in case of a single channel ID error no need to skip the rest of the ROC data -
// force geographical address and mark the produced digi
//-----------------------------------------------------------------------------
                digi_flag = mu2e::StrawDigiFlag::corrupted;
              }
              else if (hit_data->NumADCPackets != nADCPackets_) {
                int np = hit_data->NumADCPackets;
                print_(std::format("ERROR: wrong NADCpackets:{} , expected:{}, STOP PROCESSING HITS\n",
                                   np,nADCPackets_));
                break;
              }

//               uint16_t panel_id;
//               if (0 < minnesota_map_.count(mnid)){
//                 panel_id = minnesota_map_[mnid];
//               }
// //-----------------------------------------------------------------------------
// // in case of a single channel ID error no need to skip the rest of the ROC data -
// // force geographical address and mark the produced digi
// //-----------------------------------------------------------------------------
//               else{
//                 print_(std::format("ERROR: hit chid:{:04x} inconsistent with the dtc_id:{} and link_id:{}\n", hit_data->StrawIndex, dtc_id, link_id));
//                 mn_id = channel_map_[dtc_id][link_id];
//                 if (mn_id == StrawDigisFromArtdaqFragments::invalid_minnesota_){
//                   std::string msg = "encountered invalid PanelID";
//                   throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
//                 }
//                 if (minnesota_map_.count(mn_id) < 1){
//                   std::string msg = "undefined minnesota number in fallback mapping:" + std::to_string(mn_id);
//                   throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
//                 }
//                 panel_id = minnesota_map_[mn_id];
//                 digi_flag = mu2e::StrawDigiFlag::corrupted;
//               }

              if (hit_data->NumADCPackets != nADCPackets_) {
                int np = hit_data->NumADCPackets;
                print_(std::format("ERROR: wrong NADCpackets:{} , expected:{}, STOP PROCESSING HITS\n",
                                   np,nADCPackets_));
                break;
              }
//-----------------------------------------------------------------------------
// convert channel_id into a strawID
//-----------------------------------------------------------------------------
              mu2e::StrawId sid(pm->plane, pm->panel, chid);
              
              mu2e::TrkTypes::TDCValues tdc = {hit_data->TDC0(), hit_data->TDC1()};
              mu2e::TrkTypes::TOTValues tot = {hit_data->TOT0  , hit_data->TOT1  };
              mu2e::TrkTypes::ADCValue  pmp = hit_data->PMP;
              if (debugMode_ and debugBit_[1]) {
                if (header_printed == 0) {
                                        // print header
                  std::cout << "offset sid_data  mnID  plane panel    straw      TDC0       TDC1  TOT0  TOT1   PMP\n";
                  header_printed = 1;
                }
                std::cout << std::format("0x{:04x}   0x{:04x} MN{:03d}   {:3} {:3}      0x:{:04x}  {:9} {:9}   {:2}   {:2}  {:5}\n",
                                         offset,hit_data->StrawIndex,mnid,pm->plane,pm->panel,sid.straw(),hit_data->TDC0(),
                                         hit_data->TDC1(),tot[0],tot[1],pmp);
              }

              straw_digis->emplace_back(sid, tdc, tot, pmp);
//-----------------------------------------------------------------------------
// an if could be more disruptive
//-----------------------------------------------------------------------------
              auto digi = straw_digis->back();
              digi.digiFlag() = digi_flag;
//------------------------------------------------------------------------------
// the corresponding waveform
//-----------------------------------------------------------------------------
              if (saveWaveforms_) {
                std::vector<uint16_t> wf(nSamples_);

                wf[0] = hit_data->ADC00;
                wf[1] = hit_data->ADC01();
                wf[2] = hit_data->ADC02;

                auto npackets   = 0;
                auto idx        = 2;
                auto adc_packet = (mu2e::TrackerDataDecoder::TrackerADCPacket*)((char*) hit_data + 16); // the packet size is 16 bytes
                while (npackets < nADCPackets_) {
                  wf[++idx] = adc_packet->ADC0;
                  wf[++idx] = adc_packet->ADC1();
                  wf[++idx] = adc_packet->ADC2;
                  wf[++idx] = adc_packet->ADC3;
                  wf[++idx] = adc_packet->ADC4();
                  wf[++idx] = adc_packet->ADC5;
                  wf[++idx] = adc_packet->ADC6;
                  wf[++idx] = adc_packet->ADC7();
                  wf[++idx] = adc_packet->ADC8;
                  wf[++idx] = adc_packet->ADC9;
                  wf[++idx] = adc_packet->ADC10();
                  wf[++idx] = adc_packet->ADC11;
                  npackets++;
                  adc_packet++;
                }
                straw_digi_adcs->emplace_back(wf);
              }
            }
          }
//-----------------------------------------------------------------------------
// end fo ROC data processing, on to the next one
//-----------------------------------------------------------------------------
          roc_data += (nhits*np_per_hit_+1)*packet_size;
        }
      }
    }
  }

  intInfo->setNTrackerHits(straw_digis->size());
  event.put(std::move(intInfo));
//-----------------------------------------------------------------------------
// Store the straw digis in the event
//-----------------------------------------------------------------------------
  event.put(std::move(straw_digis));
  if (saveWaveforms_) {
//-----------------------------------------------------------------------------
// formatting the waveform printout takes time
// so use an aditional switch (diagLevel_)
//-----------------------------------------------------------------------------
    if (debugMode_ and (debugBit_[2] > 0)) {
//-----------------------------------------------------------------------------
// print waveforms - before moving, that invalidates the pointer...
// make sure that the case of 2 packets prints in one line, the rest is less important
//-----------------------------------------------------------------------------
      print_(std::format("--- waveforms: n:{}\n",straw_digi_adcs->size()));
      int iwf = 0;
      for (auto wf : *straw_digi_adcs) {
        std::string line = std::format("{:5d}",iwf);
        int loc = 0;
        for (int i=0; i<nSamples_; i++) {
          line += std::format(" {:5}",wf.samples()[i]);
          loc++;
          if (loc >= 27) {
            printf("%s\n",line.data());
            line = "     ";
            loc  = 0;
          }
        }
        if (loc > 0) printf("%s\n",line.data());
        iwf++;
      }
    }
    event.put(std::move(straw_digi_adcs));
  }

  if (debugMode_) print_("-- END\n",1);
}


//-----------------------------------------------------------------------------
uint16_t mu2e::StrawDigisFromArtdaqFragments::parse_minnesota_label(std::string label){
    if ((label.size() != 5) || (label[0] != 'M') || (label[1] != 'N')){
        std::string msg = "invalid minnesota label: " + label;
        throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
    }
    std::string substr = label.substr(2, 3);
    unsigned int parsed;
    int scanned = sscanf(substr.c_str(), "%u", &parsed);
    if (scanned != 1){
      std::string msg = "failed to parse minnesota label: " + label;
      throw cet::exception("StrawDigisFromArtdaqFragments") << msg << std::endl;
    }
    uint16_t rv = static_cast<uint16_t>(parsed);
    return rv;
}

// ======================================================================

DEFINE_ART_MODULE(mu2e::StrawDigisFromArtdaqFragments)

// ======================================================================
