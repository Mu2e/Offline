// ======================================================================
// clang-format off
// PM
// StrawDigisFromArtdaqFragments:  add tracker data products to the event
// debugMode_ > 0 : enables diagnostic printouts
// debugBit_[0]: print raw fragments
// debugBit_[1]: digis
// debugBit_[2]: waveforms
// ======================================================================
#include <string>

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalSequence.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/Decoders/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_RocDataHeaderPacket.h"
#include "artdaq-core-mu2e/Overlays/DTC_Packets/DTC_EventHeader.h"

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>

#include <iostream>
#include <format>

#include <regex>
#include <string>
#include <format>

#include <map>
#include <memory>

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"

// #define TRACEMF_USE_VERBATIM 1

// #include "TRACE/tracemf.h"
// #define TRACE_NAME "StrawDigisFromArtdaqFragments"


namespace mu2e {
  class StrawDigisFromArtdaqFragments;
}
// ======================================================================

class mu2e::StrawDigisFromArtdaqFragments : public art::EDProducer {

public:

  struct Config {
    fhicl::Atom<int>             diagLevel        {fhicl::Name("diagLevel"        ), fhicl::Comment("diagnostic severity level, default:0"      )};
    fhicl::Atom<int>             debugMode        {fhicl::Name("debugMode"        ), fhicl::Comment("debug mode, default:0"                     )};
    fhicl::Sequence<std::string> debugBits        {fhicl::Name("debugBits"        ), fhicl::Comment("debug bits"                                )};
    fhicl::Atom<bool>            saveWaveforms    {fhicl::Name("saveWaveforms"    ), fhicl::Comment("save StrawDigiADCWaveforms, default:true"  )};
    fhicl::Atom<bool>            missingDTCHeaders{fhicl::Name("missingDTCHeaders"), fhicl::Comment("true for runs <= 107246, default:false"    )};
    fhicl::Atom<bool>            keyOnMnid        {fhicl::Name("keyOnMnid"        ), fhicl::Comment("true if need to key on MnID, default:false")};
    fhicl::Atom<bool>            allowOfflineFallbackWhenPanelMapMissing{
      fhicl::Name("allowOfflineFallbackWhenPanelMapMissing"),
      fhicl::Comment("If TrackerPanelMap lookup fails, decode using offline StrawId(dtc,link,straw)")};
    fhicl::Atom<bool>            forceOfflineAddressing{
      fhicl::Name("forceOfflineAddressing"),
      fhicl::Comment("Ignore TrackerPanelMap/mnid and decode StrawId directly as (dtc,link,straw)")};

  };

  // --- C'tor/d'tor:
  explicit StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config);
  virtual ~StrawDigisFromArtdaqFragments() {}

  void         print_(const std::string&  Message,
                      const std::source_location& location = std::source_location::current());

  void         print_fragment(const artdaq::Fragment* Frag);

    // --- overloaded functions of the art producer
  virtual void produce (art::Event& ArtEvent) override;

  enum {kNDebugBits = 100};

private:
                                        // talk-to parameters
  int       diagLevel_    ;
  int       debugMode_   ;

  std::vector<std::string> debugBits_;
  int                      debugBit_[kNDebugBits];
  bool      saveWaveforms_;
  bool      missingDTCHeaders_;
  bool      keyOnMnid_;
  bool      allowOfflineFallbackWhenPanelMapMissing_;
  bool      forceOfflineAddressing_;
                                        // the rest
  int       nADCPackets_{-1};           // N(ADC packets per hit)
  int       nSamples_   {-1};           // N(ADC samples per hit)
  int       np_per_hit_ {-1};           // N(data packets per hits)

  const art::Event*        event_;
                                                // for now, IDTC=2*nodename+PCIE_ADDR
  ProditionsHandle<TrackerPanelMap> _tpm_h;
  const TrackerPanelMap*            _trackerPanelMap;

  };

// ======================================================================
mu2e::StrawDigisFromArtdaqFragments::StrawDigisFromArtdaqFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer   {config},
    diagLevel_        (config().diagLevel    ()),
    debugMode_        (config().debugMode    ()),
    debugBits_        (config().debugBits    ()),
    saveWaveforms_    (config().saveWaveforms()),
    missingDTCHeaders_(config().missingDTCHeaders()),
    keyOnMnid_        (config().keyOnMnid()),
    allowOfflineFallbackWhenPanelMapMissing_(config().allowOfflineFallbackWhenPanelMapMissing()),
    forceOfflineAddressing_(config().forceOfflineAddressing()),
    event_            (nullptr)
{
  produces<mu2e::StrawDigiCollection>();
  if (saveWaveforms_) produces<mu2e::StrawDigiADCWaveformCollection>();

  produces<mu2e::IntensityInfoTrackerHits>();

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

    print_(std::format("StrawDigisFromArtdaqFragments: bit={:4d} is set to {}",index,debugBit_[index]));
  }
}


std::vector<std::string> splitString(const std::string& str, const std::string& delimiter) {
    std::vector<std::string> result;
    std::regex re(delimiter);
    std::sregex_token_iterator it(str.begin(), str.end(), re, -1);
    std::sregex_token_iterator end;
    while (it != end) {
        result.push_back(*it++);
    }
    return result;
}

//-----------------------------------------------------------------------------
void mu2e::StrawDigisFromArtdaqFragments::print_(const std::string& Message, const std::source_location& location) {

  std::string s;
  if (event_) s = std::format("event: {}:{}:{} ",event_->run(),event_->subRun(),event_->event());

  std::vector<std::string> ss = splitString(location.file_name(),"/");

  mf::LogVerbatim("MAKE_SD")
     << s << ss.back() << ":" << location.line()
     //            << location.function_name()
     << " : " << Message;

  // std::cout << s << ss.back() << ":" << location.line()
  //   //            << location.function_name()
  //      << " : " << Message << std::endl;
}

//-----------------------------------------------------------------------------
// HEX print of a fragment, the Mu2e data come in 2-byte words
//-----------------------------------------------------------------------------
void mu2e::StrawDigisFromArtdaqFragments::print_fragment(const artdaq::Fragment* Frag) {
  ushort* buf = (ushort*) (Frag->dataBegin());
  int nw      = Frag->dataSizeBytes()/2;
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

// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void mu2e::StrawDigisFromArtdaqFragments::produce(art::Event& event) {
  int const packet_size(16); // in bytes

  if (debugMode_ > 0) print_("-- START");

  event_ = &event;                      // cache to print events

  _trackerPanelMap = &_tpm_h.get(event.id());

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

  if (debugMode_ > 0) {
    std::string msg = std::format("n_fragment_collections:{}",fragmentHandles.size());
    print_(msg);
  }

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
      int n_fragments = handle->size();

      if (debugMode_) {
        print_(std::format("-- next fragment collection with n_fragments:{}",n_fragments));
      }

      for (int ifrag=0; ifrag<n_fragments; ifrag++) {
        const artdaq::Fragment* frag = &handle->at(ifrag);

        if (debugMode_ and (debugBit_[0] > 0)) {
          print_(std::format("-- fragment number:{} version:{} timestamp:{} data_size:{} type:{} DTC_SubEventHeader.size:{}",
                             ifrag,frag->version(),frag->timestamp(),frag->dataSizeBytes(),
                             frag->typeString(),sizeof(DTCLib::DTC_SubEventHeader)));
          print_fragment(frag);
        }
//-----------------------------------------------------------------------------
// skip CFO fragment (type = 12)
//-----------------------------------------------------------------------------
        if (frag->type() == mu2e::FragmentType::CFO)        continue;
        uint8_t* fdata = (uint8_t*) (frag->dataBegin());
        if (not missingDTCHeaders_) {
//-----------------------------------------------------------------------------
// skip fragments with the payload size less than the DTC header size
// do it only for the current data format (runs > 107236)
//-----------------------------------------------------------------------------
          if (frag->dataSizeBytes() <= sizeof(DTCLib::DTC_SubEventHeader)) {
            std::string msg = std::format("ERROR: fragment:{} data size:{} < DTC_SubEventHeader.size:{}. SKIP FRAGMENT",
                                          ifrag,frag->dataSizeBytes(),sizeof(DTCLib::DTC_SubEventHeader));
            print_(msg);
            continue;
          }
          fdata += sizeof(DTCLib::DTC_EventHeader);
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
        uint32_t dtc_id = seh->source_dtc_id;
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
// 2026-03-20: don't store waveforms if only one packet per hit
//-----------------------------------------------------------------------------
              if (roc_data+packet_size >= last_address) {
                print_(std::format("ERROR: dtc_id:{} roc_data:{} last_address:{} , SKIPPING",
                                   dtc_id, (void*) roc_data, (void*) last_address));
                break;
              }
              mu2e::TrackerDataDecoder::TrackerDataPacket* h0;
              h0           = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+packet_size);
//-----------------------------------------------------------------------------
// don't expect zero ADC packets, the error seen was that one of the ROCs just
// doesn't send the second packet at all... work under that assumption
//-----------------------------------------------------------------------------
              if (h0->NumADCPackets == 0) {
                print_(std::format("ERROR: dtc_id:{} link_id:{} N(ADC packets) = 0, skip ROC data",
                                   dtc_id,(int) rdh->linkID));

                roc_data += (rdh->packetCount+1)*packet_size;
                continue;
              }
              nADCPackets_ = h0->NumADCPackets;
              nSamples_    = 3+12*nADCPackets_;
              np_per_hit_  = nADCPackets_+1;
            }
            uint32_t link_id = rdh->linkID;
            nhits            = rdh->packetCount/(nADCPackets_+1);
//-----------------------------------------------------------------------------
// there should not be more than 255 hits per ROC, if nhits>255 it is a corruption,
// stop processing of the event
//-----------------------------------------------------------------------------
            if (nhits > 255) {
              print_(std::format("ERROR: nhits:{}, skip DTC",nhits));
              break;
            }

            const TrkPanelMap::Row* tpm(nullptr);
            if (!forceOfflineAddressing_ && not keyOnMnid_) {
              tpm = _trackerPanelMap->panel_map_by_online_ind(dtc_id,link_id);
              if (tpm == nullptr) {
                if (!allowOfflineFallbackWhenPanelMapMissing_) {
//-----------------------------------------------------------------------------
// either DTC ID or link ID are corrupted. Haven't seen that so far, switch to the next ROC anyway
//-----------------------------------------------------------------------------
                  print_(std::format("ERROR: either dtc_id:{} or link_id:{} is corrupted, skip ROC data",
                                     dtc_id,link_id));

                  roc_data += (nhits*np_per_hit_+1)*packet_size;
                  continue;
                }
                print_(std::format("WARNING: no panel map for dtc_id:{} link_id:{}, using offline fallback",
                                   dtc_id,link_id));
              }
            }

            if (debugMode_) {
              print_(std::format("-- DTC:{} ROC:{} nhits:{}",dtc_id,link_id,nhits));
            }

            for (int ihit=0; ihit<nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
              mu2e::TrackerDataDecoder::TrackerDataPacket* hit_data ;

              int offset = (ihit*np_per_hit_+1)*packet_size;   // in bytes
              hit_data   = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+offset);
              if (roc_data+offset >= last_address) {
                print_(std::format("ERROR: dtc_id:{} link_id:{} roc_data:{} offset:{} last_address:{} , SKIPPING",
                                   dtc_id, link_id, (void*) roc_data, offset, (void*) last_address));
                break;
              }
//-----------------------------------------------------------------------------
// at this point, check consistency between the channel_id, dtc_id and link_id for a given run
// panel ID is a derivative of the DTC ID and the link iD
// mn_id - 'MinnesotaID' of the panel
//-----------------------------------------------------------------------------
              mu2e::StrawDigiFlag digi_flag;
              uint16_t channel = static_cast<uint16_t>(hit_data->StrawIndex);
              uint16_t chid   = mu2e::StrawId(channel).straw(); // channel ID within the panel

              if (chid >= StrawId::_nstraws) {
                if (debugBit_[52] == 0) print_(std::format("ERROR: hit with corrupted chid:{:04x} : straw:{} / dtc_id:{} link_id:{}, SKIPPING",
                                                          hit_data->StrawIndex, chid, dtc_id, link_id));
                continue;
              }

              uint16_t mnid    = channel >> mu2e::StrawId::_panelsft;

              if (!forceOfflineAddressing_ && keyOnMnid_) {
                tpm = _trackerPanelMap->panel_map_by_mnid(mnid);
                if (tpm == nullptr) {
                  if (!allowOfflineFallbackWhenPanelMapMissing_) {
//-----------------------------------------------------------------------------
// bad mnid. Likely, corrupted data block. For now, skip the hit data and proceed with the next hit
//-----------------------------------------------------------------------------
                    if (debugBit_[51] == 0) print_(std::format("ERROR: corrupted mnid:{}, skip hit data",mnid));
                    continue;
                  }
                  print_(std::format("WARNING: no panel map for mnid:{}, using offline fallback",mnid));
                }
              }
// in principle, could this could become an 'else if'
              if (!forceOfflineAddressing_ && tpm != nullptr && tpm->mnid() != mnid) {
                print_(std::format("ERROR: mnid:{:3d} tpm->mnid():{:3d} hit chid:{:04x} inconsistent with the dtc_id:{:2d} and link_id:{}",
                                   mnid,tpm->mnid(),hit_data->StrawIndex, dtc_id, link_id));
//-----------------------------------------------------------------------------
// in case of a single channel ID error no need to skip the rest of the ROC data -
// force geographical address and mark the produced digi
//-----------------------------------------------------------------------------
                digi_flag = mu2e::StrawDigiFlag::corrupted;
              }

              if (hit_data->NumADCPackets != nADCPackets_) {
                int np = hit_data->NumADCPackets;
                print_(std::format("ERROR: wrong NADCpackets:{} , expected:{}, GO TO THE NEXT ROC",
                                   np,nADCPackets_));
                break;
              }
//-----------------------------------------------------------------------------
// convert channel_id into a strawID
//-----------------------------------------------------------------------------
              mu2e::StrawId sid = (forceOfflineAddressing_ || tpm == nullptr)
                ? mu2e::StrawId(dtc_id, link_id, chid)
                : mu2e::StrawId(tpm->uniquePlane(), tpm->panel(), chid);

              mu2e::TrkTypes::TDCValues tdc = {hit_data->TDC0(), hit_data->TDC1()};
              mu2e::TrkTypes::TOTValues tot = {hit_data->TOT0  , hit_data->TOT1  };
              mu2e::TrkTypes::ADCValue  pmp = hit_data->PMP;
              if (debugMode_ and debugBit_[1]) {
                if (header_printed == 0) {
                                        // print header
                  std::cout << "index offset sid_data  mnID  plane panel    straw      TDC0       TDC1  TOT0  TOT1   PMP\n";
                  header_printed = 1;
                }

                auto const planeForPrint = (tpm != nullptr) ? tpm->uniquePlane() : dtc_id;
                auto const panelForPrint = (tpm != nullptr) ? tpm->panel() : link_id;

                int ind = straw_digis->size();
                std::cout << std::format("{:5} 0x{:04x}   0x{:04x} MN{:03d}   {:3} {:3}      0x:{:04x}  {:9} {:9}   {:2}   {:2}  {:5}\n",
                                         ind,offset,hit_data->StrawIndex,mnid,planeForPrint,panelForPrint,sid.straw(),hit_data->TDC0(),
                                         hit_data->TDC1(),tot[0],tot[1],pmp);
              }

              straw_digis->emplace_back(sid, tdc, tot, pmp);
//-----------------------------------------------------------------------------
// an if could be more disruptive
//-----------------------------------------------------------------------------
              auto digi = straw_digis->back();
              digi.digiFlag() = digi_flag;
//------------------------------------------------------------------------------
// the corresponding waveform, store only if at least the second packet is present (nSamples_ = 15 or more)
//-----------------------------------------------------------------------------
              if (saveWaveforms_) {
                if (nSamples_ <= 3) {
                  print_(std::format("ERROR: nSamples:{}, do not store waveforms",nSamples_));
                }
                else {
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

  if (debugMode_) print_("-- END");
}

// ======================================================================

DEFINE_ART_MODULE(mu2e::StrawDigisFromArtdaqFragments)

// ======================================================================
