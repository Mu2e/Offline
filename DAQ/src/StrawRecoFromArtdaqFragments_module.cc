// ======================================================================
//
// StrawRecoFromArtdaqFragments:  add tracker data products to the event
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Data/TrackerDataDecoder.hh"
#include "artdaq-core-mu2e/Overlays/FragmentType.hh"

#include "Offline/DataProducts/inc/TrkTypes.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

#include <artdaq-core/Data/Fragment.hh>
#include <artdaq-core/Data/ContainerFragment.hh>

#include <iostream>

#include <string>

#include <memory>

#define TRACEMF_USE_VERBATIM 1

#include "TRACE/tracemf.h"
#define TRACE_NAME "StrawRecoFromArtdaqFragments"

namespace art {
  class StrawRecoFromArtdaqFragments;
}
// ======================================================================

class art::StrawRecoFromArtdaqFragments : public EDProducer {

public:

  struct RocDataHeaderPacket_t {            // 8 16-byte words in total
                                            // 16-bit word 0
    uint16_t            byteCount    : 16;
                                            // 16-bit word 1
    uint16_t            unused       : 4;
    uint16_t            packetType   : 4;
    uint16_t            linkID       : 3;
    uint16_t            DtcErrors    : 4;
    uint16_t            valid        : 1;
                                            // 16-bit word 2
    uint16_t            packetCount  : 11;
    uint16_t            unused2      : 2;
    uint16_t            subsystemID  : 3;
                                            // 16-bit words 3-5
    uint16_t            eventTag[3];
                                            // 16-bit word 6
    uint8_t             status       : 8;
    uint8_t             version      : 8;
                                            // 16-bit word 7
    uint8_t             dtcID        : 8;
    uint8_t             onSpill      : 1;
    uint8_t             subrun       : 2;
    uint8_t             eventMode    : 5;
                                            // decoding status

    int                 empty     () { return (status & 0x01) == 0; }
    int                 invalid_dr() { return (status & 0x02); }
    int                 corrupt   () { return (status & 0x04); }
    int                 timeout   () { return (status & 0x08); }
    int                 overflow  () { return (status & 0x10); }

    int                 error_code() { return (status & 0x1e); }
  };

  struct Config {
    fhicl::Atom<int> diagLevel    {fhicl::Name("diagLevel"    ), fhicl::Comment("diagnostic level"          )};
    fhicl::Atom<int> saveWaveforms{fhicl::Name("saveWaveforms"), fhicl::Comment("save StrawDigiADCWaveforms")};
    fhicl::Atom<int> nADCPackets  {fhicl::Name("nADCPackets"  ), fhicl::Comment("N(ADC packets per hit)"    )};
  };

  // --- C'tor/d'tor:
  explicit StrawRecoFromArtdaqFragments(const art::EDProducer::Table<Config>& config);
  virtual ~StrawRecoFromArtdaqFragments() {}

  // --- Production:
  virtual void produce(Event&);

private:

  int       diagLevel_;
  int       saveWaveforms_;
  int       nADCPackets_;

}; // StrawRecoFromArtdaqFragments

// ======================================================================

art::StrawRecoFromArtdaqFragments::StrawRecoFromArtdaqFragments(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    diagLevel_    (config().diagLevel    ()),
    saveWaveforms_(config().saveWaveforms()),
    nADCPackets_  (config().nADCPackets  ())
{
  produces<mu2e::StrawDigiCollection>();
  if (saveWaveforms_) produces<mu2e::StrawDigiADCWaveformCollection>();

  produces<mu2e::IntensityInfoTrackerHits>();
  // FIXME!
  produces<mu2e::ProtonBunchTime>();
}

// ----------------------------------------------------------------------
// runs on tracker Artdaq fragments
//-----------------------------------------------------------------------------
void art::StrawRecoFromArtdaqFragments::produce(Event& event) {
  int const packet_size(16); // in bytes

  TLOG(TLVL_DEBUG) << "--- START event:" << event.run() << ":" << event.subRun() << ":" << event.event();

   // Collection of StrawDigis for the event
  std::unique_ptr<mu2e::StrawDigiCollection> straw_digis(new mu2e::StrawDigiCollection);
  std::unique_ptr<mu2e::StrawDigiADCWaveformCollection> straw_digi_adcs(
                 new mu2e::StrawDigiADCWaveformCollection);

  // IntensityInfoTrackerHits
  std::unique_ptr<mu2e::IntensityInfoTrackerHits> intInfo(new mu2e::IntensityInfoTrackerHits);

  // FIXME! this is temporary

  std::unique_ptr<mu2e::ProtonBunchTime> pbt(new mu2e::ProtonBunchTime);
  pbt->pbtime_ = 0;
  pbt->pbterr_ = 0;
  event.put(std::move(pbt));

  int np_per_hit  = 1+nADCPackets_;
  int nsamples    = 3+12*nADCPackets_;

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
        const auto& frag = handle->at(ifrag);
        uint8_t* fdata  = (uint8_t*) (frag.dataBegin());
        TLOG(TLVL_DEBUG) << "-- fragment number:" << ifrag;
//-----------------------------------------------------------------------------
// skip non-tracker fragments
//-----------------------------------------------------------------------------
        DTCLib::DTC_SubEventHeader* seh = (DTCLib::DTC_SubEventHeader*) fdata;
        if (seh->source_subsystem != DTCLib::DTC_Subsystem::DTC_Subsystem_Tracker) continue;
        int dtc_id  = seh->source_dtc_id;
//-----------------------------------------------------------------------------
// this is a tracker DTC fragment, loop over the ROCs
//-----------------------------------------------------------------------------
        ushort*  buf          = (ushort*) fdata;
        int      nbytes       = buf[0]; // frag.dataSizeBytes() includes extra 0x20
        uint8_t* roc_data     = fdata+sizeof(*seh);
        uint8_t* last_address = fdata+nbytes;     //

        while (roc_data < last_address) {
          RocDataHeaderPacket_t* rdh = (RocDataHeaderPacket_t*) roc_data;

          int nhits = rdh->packetCount/(nADCPackets_+1); // nADCPackets_ defaults to 1

          int link_id = rdh->linkID;
          TLOG(TLVL_DEBUG) << "DTC:" << dtc_id << " ROC:" << link_id
                           << " offset:" << std::hex
                           << " nhits:" << std::dec << nhits;

          for (int ihit=0; ihit<nhits; ihit++) {
//-----------------------------------------------------------------------------
// first packet, 16 bytes, or 8 ushort's is the data header packet
//-----------------------------------------------------------------------------
            mu2e::TrackerDataDecoder::TrackerDataPacket* hit ;
            int offset          = (ihit*np_per_hit+1)*packet_size;   // in bytes
            hit = (mu2e::TrackerDataDecoder::TrackerDataPacket*) (roc_data+offset);
//-----------------------------------------------------------------------------
// unpack hit and fill the StrawDigiCollection
// HV side has an extra bit in teh straw ID
//-----------------------------------------------------------------------------
            uint16_t straw_id = hit->StrawIndex;
            if (straw_id >= 0x80) straw_id = straw_id - 0x80;

            mu2e::StrawId sid(hit->StrawIndex);
            mu2e::TrkTypes::TDCValues tdc = {hit->TDC0(), hit->TDC1()};
            mu2e::TrkTypes::TOTValues tot = {hit->TOT0, hit->TOT1};
            mu2e::TrkTypes::ADCValue  pmp = hit->PMP;

            TLOG(TLVL_DEBUG+1) << "times:" << std::setw(10) << hit->TDC0() << " " << hit->TDC1()
                               << " TOT:" << std::setw(4) << hit->TOT0 << " " << hit->TOT1
                               << " pmp:" << hit->PMP;

            straw_digis->emplace_back(sid, tdc, tot, pmp);
//------------------------------------------------------------------------------
// the corresponding waveform
//-----------------------------------------------------------------------------
            if (saveWaveforms_) {
              std::vector<uint16_t> wf(nsamples);

              wf[0] = hit->ADC00;
              wf[1] = hit->ADC01();
              wf[2] = hit->ADC02;

              auto npackets   = 0;
              auto idx        = 2;
              auto adc_packet = (mu2e::TrackerDataDecoder::TrackerADCPacket*)((char*) hit + 16); // the packet size is 16 bytes
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
          roc_data += (nhits*np_per_hit+1)*packet_size;
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
    if (diagLevel_ > 1) {
      // printout better be matched to the chID's - but OK for now
      for (auto& wf : *straw_digi_adcs) {
        int loc = 0;
        for (int i=0; i<nsamples; i++) {
          std::cout << std::format(" 0x{:4x}",wf.samples()[i]);
          if (loc >= 15) {
            std::cout << "\n";
            loc = 0;
          }
        }
        if (loc > 0) std::cout << "\n";
      }
    }
    event.put(std::move(straw_digi_adcs));
  }

  TLOG(TLVL_DEBUG) << "--- END";
}



// ======================================================================

DEFINE_ART_MODULE(art::StrawRecoFromArtdaqFragments)

// ======================================================================
