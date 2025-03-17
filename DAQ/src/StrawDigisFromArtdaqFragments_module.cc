// ======================================================================
// PM
// StrawDigisFromArtdaqFragments:  add tracker data products to the event
// each diagnostic printout has a level >= 0
// b) the higher it is, the less important is the printout
// c) in a job with diagLevel_ set, only printouts with level <= diagLevel_ are enabled
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

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


namespace art {
  class StrawDigisFromArtdaqFragments;
}
// ======================================================================

class art::StrawDigisFromArtdaqFragments : public EDProducer {

public:

  struct Config {
    fhicl::Atom<int> diagLevel    {fhicl::Name("diagLevel"    ), fhicl::Comment("diagnostic severity level, default = 0" ), 0};
    fhicl::Atom<int> debugLevel   {fhicl::Name("debugLevel"   ), fhicl::Comment("debug level, default = 0"               ), 0};
    fhicl::Atom<int> saveWaveforms{fhicl::Name("saveWaveforms"), fhicl::Comment("save StrawDigiADCWaveforms, default = 1"), 1};
  };

  // --- C'tor/d'tor:
  explicit StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config);
  virtual ~StrawDigisFromArtdaqFragments() {}

  int          panelID(uint16_t ChannelID, int DtcID, int LinkID);

  void         print_(const std::string&  Message, int DiagLevel = -1,
                      const std::source_location& location = std::source_location::current());

  void         print_fragment(const artdaq::Fragment* Frag);

    // --- overloaded functions of the art producer
  virtual void produce (art::Event& ArtEvent) override;
  virtual void beginRun(art::Run&   ArtRun  ) override;

private:
                                        // talk-to parameters
  int       diagLevel_    ;
  int       debugLevel_   ;
  int       saveWaveforms_;
                                        // the rest
  int       nADCPackets_{-1};           // N(ADC packets per hit)
  int       nSamples_   {-1};           // N(ADC samples per hit)
  int       np_per_hit_ {-1};           // N(data packets per hits)

  int       panel_map_[36][6];          // panel_map_[idtc][ilink] = panel MN number

  const art::Event* event_;

};

// ======================================================================
art::StrawDigisFromArtdaqFragments::StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_    (config().diagLevel    ()),
    debugLevel_   (config().debugLevel   ()),
    saveWaveforms_(config().saveWaveforms())
{
  produces<mu2e::StrawDigiCollection>();
  if (saveWaveforms_) produces<mu2e::StrawDigiADCWaveformCollection>();

  produces<mu2e::IntensityInfoTrackerHits>();
  // FIXME!
  //  produces<mu2e::ProtonBunchTime>();
}


//-----------------------------------------------------------------------------
void art::StrawDigisFromArtdaqFragments::print_(const std::string& Message, int DiagLevel,
                                                const std::source_location& location) {
  if (DiagLevel > diagLevel_) return;
  std::cout << std::format(" event:{}:{}:{}",event_->run(),event_->subRun(),event_->event())
            << " " << location.file_name() << ":" << location.line()
    //            << location.function_name()
            << ": " << Message << std::endl;
}

//-----------------------------------------------------------------------------
// HEX print of a fragment, the data has to be in 2-byte words
//-----------------------------------------------------------------------------
void art::StrawDigisFromArtdaqFragments::print_fragment(const artdaq::Fragment* Frag) {
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
// to be written - needs DB access
// for now - all good
// return panel ID part of the geographical (offline) channel ID - first channel ID of the panel
//-----------------------------------------------------------------------------
int art::StrawDigisFromArtdaqFragments::panelID(uint16_t MnID, int DtcID, int LinkID) {
  int panel_id(-1), idtc(-1);

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
//-----------------------------------------------------------------------------
// for the moment, assign unique panel ID's just to get through,
//-----------------------------------------------------------------------------
  panel_id = idtc*6+LinkID;

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
  return panel_id;
}

//-----------------------------------------------------------------------------
void art::StrawDigisFromArtdaqFragments::beginRun(art::Run&  ArtRun) {
  // fill panel_map_ for a given run - should come from the database

  for (int idtc=0; idtc<36; idtc++) {
    for (int ilink=0; ilink<6; ilink++) {
      panel_map_[idtc][ilink] = idtc*6+ilink;
    }
  }

}
// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void art::StrawDigisFromArtdaqFragments::produce(Event& event) {
  int const packet_size(16); // in bytes

  event_ = &event;                      // cache for printouts
  print_("-- START",1);

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
        // if (contf.fragment_type() != mu2e::FragmentType::DTCEVT) {
        //   break;
        // }

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

        print_(std::format("-- fragment number:{}",ifrag),1);
        if (debugLevel_ & 0x1) {
                                        // debug: print fragment
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

            print_(std::format("--- DTC:{} ROC:{} nhits:{}",dtc_id,link_id,nhits),1);
            for (int ihit=0; ihit<nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
              mu2e::TrackerDataDecoder::TrackerDataPacket* hit_data ;
              int offset          = (ihit*np_per_hit_+1)*packet_size;   // in bytes
              hit_data = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+offset);
//-----------------------------------------------------------------------------
// at this point, check consistency between the channel_id, dtc_id and link_id for a given run
// panel ID is a derivative of the DTC ID and the link iD
// mn_id - 'MinnesotaID' of the panel
//-----------------------------------------------------------------------------
              mu2e::StrawDigiFlag digi_flag;
              int ch_id    = hit_data->StrawIndex & 0x7f;   // channel ID within the panel
              int mn_id    = hit_data->StrawIndex >> 7;
              int panel_id = panelID(mn_id,dtc_id,link_id);
              if (panel_id < 0) {
                print_(std::format("ERROR: hit chid:{:04x} inconsistent with the dtc_id:{} and link_id:{}",
                                   hit_data->StrawIndex, dtc_id, link_id));
//-----------------------------------------------------------------------------
// in case of a single channel ID error no need to skip the rest of the ROC data -
// force geographical address and mark the produced digi
//-----------------------------------------------------------------------------
                panel_id  = panel_map_[dtc_id][link_id];
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
              uint16_t straw_index = (panel_id << 7) | ch_id;
              mu2e::StrawId sid(straw_index);
              mu2e::TrkTypes::TDCValues tdc = {hit_data->TDC0(), hit_data->TDC1()};
              mu2e::TrkTypes::TOTValues tot = {hit_data->TOT0, hit_data->TOT1};
              mu2e::TrkTypes::ADCValue  pmp = hit_data->PMP;

              print_(std::format("offset:0x{:04x} sid:{:5} times: {:9} {:9} TOT:{:2}:{:2} pmp:{}",
                                 offset, straw_index,hit_data->TDC0(),hit_data->TDC1(),tot[0],tot[1],pmp),1);

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
    if (debugLevel_ & 0x2) {
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

  print_("-- END",1);
}



// ======================================================================

DEFINE_ART_MODULE(art::StrawDigisFromArtdaqFragments)

// ======================================================================
