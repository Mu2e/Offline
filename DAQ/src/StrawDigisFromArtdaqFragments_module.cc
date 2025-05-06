// ======================================================================
// PM
// StrawDigisFromArtdaqFragments:  add tracker data products to the event
// debugLevel_ > 0 : enables diagnostic printouts
// diagLevel_ : bitmask, determines the format of the printout
// ======================================================================
#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_RocDataHeaderPacket.h"

#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
// #include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>

#include <iostream>

#include <string>

#include <memory>

// #define TRACEMF_USE_VERBATIM 1

// #include "TRACE/tracemf.h"
// #define TRACE_NAME "StrawDigisFromArtdaqFragments"

#include "Offline/DAQ/inc/TrkPanelMap_t.hh"

namespace mu2e {

// ======================================================================

  class StrawDigisFromArtdaqFragments : public art::EDProducer {

public:

  struct Config {
    fhicl::Atom<int>             diagLevel    {fhicl::Name("diagLevel"    ), fhicl::Comment("diagnostic severity level, default = 0" )};
    fhicl::Atom<int>             debugLevel   {fhicl::Name("debugLevel"   ), fhicl::Comment("debug level, default = 0"               )};
    fhicl::Sequence<std::string> debugBits    {fhicl::Name("debugBits"    ), fhicl::Comment("debug bits"                             )};
    fhicl::Atom<int>             saveWaveforms{fhicl::Name("saveWaveforms"), fhicl::Comment("save StrawDigiADCWaveforms, default = 1")};
  };

  // --- C'tor/d'tor:
  explicit StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config);
  virtual ~StrawDigisFromArtdaqFragments() {}

  int          panelID(uint16_t ChannelID, int DtcID, int LinkID);

  void         print_(const std::string&          Message,
                      const std::source_location& location = std::source_location::current());

  void         print_fragment(const artdaq::Fragment* Frag);

    // --- overloaded functions of the art producer
  virtual void produce (art::Event& ArtEvent) override;
  virtual void beginRun(art::Run&   ArtRun  ) override;

  int offlineDtcID(int DtcID);
      
  private:
                                        // talk-to parameters
  int                      diagLevel_    ;
  int                      debugLevel_   ;
  std::vector<std::string> debugBits_;
  int                      debugBit_[100];  
  int                      saveWaveforms_;
                                        // the rest
  int                      nADCPackets_{-1};           // N(ADC packets per hit)
  int                      nSamples_   {-1};           // N(ADC samples per hit)
  int                      np_per_hit_ {-1};           // N(data packets per hits)

  const TrkPanelMap_t*     panel_map_[36][6];          // panel_map_[idtc][ilink] = panel MN number

  const art::Event* event_;

};

// ======================================================================
  StrawDigisFromArtdaqFragments::StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_    (config().diagLevel    ()),
    debugLevel_   (config().debugLevel   ()),
    debugBits_    (config().debugBits    ()),
    saveWaveforms_(config().saveWaveforms()),
    event_        (nullptr)                   {
      produces<mu2e::StrawDigiCollection>();
      if (saveWaveforms_) produces<mu2e::StrawDigiADCWaveformCollection>();

      produces<mu2e::IntensityInfoTrackerHits>();
  // FIXME!
  //  produces<mu2e::ProtonBunchTime>();
//-----------------------------------------------------------------------------
// parse debug bits
//-----------------------------------------------------------------------------
      const char* key;
                                        // a flag is an integer!
      int nbits = debugBits_.size();
      for (int i=0; i<nbits; i++) {
        int index(0), value(0);
        key               = debugBits_[i].data();
        sscanf(key,"bit%i:%i",&index,&value);
        debugBit_[index]  = value;
        
        print_(std::format("... StrawDigisFromArtdaqFragments: bit={:4d} is set to {}\n",index,debugBit_[index]));
      }
  }

//-----------------------------------------------------------------------------
void StrawDigisFromArtdaqFragments::print_(const std::string& Message, const std::source_location& location) {
  if (event_) std::cout << std::format(" event:{}:{}:{}",event_->run(),event_->subRun(),event_->event());
  std::cout << " " << location.file_name() << ":" << location.line()
    //            << location.function_name()
            << ": " << Message << std::endl;
}

//-----------------------------------------------------------------------------
// HEX print of a fragment, the data has to be in 2-byte words
//-----------------------------------------------------------------------------
void StrawDigisFromArtdaqFragments::print_fragment(const artdaq::Fragment* Frag) {
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
// make sure that can handle the test stand in EC and the station in Lab3
// to be written - the DB access
// returns truncated DTC ID - 0 or 1, don't need more for now - don't read out more
// thank one station
//-----------------------------------------------------------------------------
int StrawDigisFromArtdaqFragments::offlineDtcID(int DtcID) {
  int idtc(-1);

  if      (DtcID == 18) idtc = 0;   // daq09
  else if (DtcID == 19) idtc = 1;   // daq09
  else if (DtcID == 42) idtc = 1;   // daq09 old
  else if (DtcID == 73) idtc = 0;   // daq09 old
  else if (DtcID == 44) idtc = 0;   // daq22
  else if (DtcID == 45) idtc = 1;   // daq22
  else {
                                    // unknown DTC, return an error
    print_(std::format("ERROR: unknown DtcID:{}",DtcID));
    return -1;
  }

//   if (MnID == panel_map_[idtc][LinkID]) return MnID;
//   else {
// //-----------------------------------------------------------------------------
// //  handle test stands and the tower
// //  panelID 400-405 : the tower (bottom to top)
// //  panelID 410     : TS0
// //  panelID 411     : TS1
// //  panelID 412     : TS2
// //  it would be interesting to see MnID's resulting in a failure
// //-----------------------------------------------------------------------------
//     if ((MnID >= 400) and (MnID < 412)) {
//       panel_id = MnID-400;
//     }
//     else {
//       print_(std::format("ERROR: Minnesota ID:{:04x} inconsistent with the dtc_id:{} and link_id:{}",
//                          MnID, DtcID, LinkID));
//     }
//   }
  return idtc;
}

//-----------------------------------------------------------------------------
void StrawDigisFromArtdaqFragments::beginRun(art::Run&  ArtRun) {
  // fill panel_map_ for a given run - should come from the database
  
    for (const TrkPanelMap_t* tpm = TrkPanelMap_data.begin(); tpm != TrkPanelMap_data.end(); ++tpm) {
      int dtc  = tpm->dtc;
      int link = tpm->link;
      panel_map_[dtc][link] = tpm;
    }

}
 
// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
 void StrawDigisFromArtdaqFragments::produce(art::Event& event) {
  int const packet_size(16); // in bytes

  event_ = &event;                      // cache for printouts
  if (debugLevel_ > 0) print_("-- START");

   // Collection of StrawDigis for the event
  std::unique_ptr<mu2e::StrawDigiCollection> straw_digis(new mu2e::StrawDigiCollection);
  std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> straw_digi_adcs(new mu2e::StrawDigiADCWaveformCollection);

  // IntensityInfoTrackerHits
  std::unique_ptr<mu2e::IntensityInfoTrackerHits> intInfo(new mu2e::IntensityInfoTrackerHits);

  // FIXME! this is temporary
  // std::unique_ptr<mu2e::ProtonBunchTime> pbt(new mu2e::ProtonBunchTime);
  // pbt->pbtime_ = 0;
  // pbt->pbterr_ = 0;
  //  event.put(std::move(pbt));
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
        if (debugLevel_) {
          if (diagLevel_ > 0) print_(std::format("-- fragment number:{}",ifrag));
          
                                        // debugLevel bit0: print fragment
          if (diagLevel_ & 0x1) print_fragment(frag);
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
        int dtc_id = offlineDtcID(seh->source_dtc_id);
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

            if (debugLevel_) {
              print_(std::format("--- DTC:{} ROC:{} nhits:{}",dtc_id,link_id,nhits));
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
              int ch_id         = hit_data->StrawIndex & 0x7f;   // channel ID within the panel
              int mnid          = hit_data->StrawIndex >> 7;
              const TrkPanelMap_t* pm = panel_map_[dtc_id][link_id];  // DB here
              if (pm->mnid != mnid) {
                print_(std::format("ERROR: hit chid:{:04x} inconsistent with the dtc_id:{} and link_id:{}",
                                   hit_data->StrawIndex, dtc_id, link_id));
//-----------------------------------------------------------------------------
// in case of a single channel ID error no need to skip the rest of the ROC data -
// force geographical address and mark the produced digi
//-----------------------------------------------------------------------------
                digi_flag = mu2e::StrawDigiFlag::corrupted;
              }
              else if (hit_data->NumADCPackets != nADCPackets_) {
                int np = hit_data->NumADCPackets;
                print_(std::format("ERROR: wrong NADCpackets:{} , expected:{}, STOP PROCESSING HITS",
                                   np,nADCPackets_));
                break;
              }
//-----------------------------------------------------------------------------
// convert channel_id into a strawID
//-----------------------------------------------------------------------------
              mu2e::StrawId sid(pm->plane, pm->panel, ch_id);
              
              mu2e::TrkTypes::TDCValues tdc = {hit_data->TDC0(), hit_data->TDC1()};
              mu2e::TrkTypes::TOTValues tot = {hit_data->TOT0, hit_data->TOT1};
              mu2e::TrkTypes::ADCValue  pmp = hit_data->PMP;
              if (debugLevel_ and debugBit_[1]) {
                if (header_printed == 0) {
                                        // print header
                  std::cout << " offset sid_data  mnid ch_id panel_id sid_ofln      TDC0    TDC1   TOT0 TOT1   PMP\n";
                  header_printed = 1;
                }
                std::cout << std::format("0x{:04x} 0x{:04x} MN{:03d} {:5}  {:3}  0x:{:04x}  {:9} {:9} {:2}  {:2} {:5}\n",
                                         offset,hit_data->StrawIndex,mnid,ch_id,pm->panel,sid.straw(),hit_data->TDC0(),
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
    if (debugLevel_ and (diagLevel_ & 0x2)) {
//-----------------------------------------------------------------------------
// print waveforms - before moving, that invalidates the pointer...
// make sure that the case of 2 packets prints in one line, the rest is less important
//-----------------------------------------------------------------------------
      print_(std::format("--- waveforms: n:{}",straw_digi_adcs->size()));
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

  if (debugLevel_) print_("-- END");
}

}

// ======================================================================

DEFINE_ART_MODULE(mu2e::StrawDigisFromArtdaqFragments)

// ======================================================================
