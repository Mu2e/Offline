# DQMHelpers

Shared DQM histogram helpers used by both offline art analyzers and the
otsdaq online monitor. The helper classes own histogram booking and filling;
art modules own I/O, visualization, and histogram shipping.

This package has no GeometryService or ProditionsService dependency, so the
same fill path can run in the DAQ process.

## CRVDigiDQM

`mu2e::CRVDigiDQM` books and fills every live digi histogram from
`otsdaq-mu2e-crv` `CrvDQM_module.cc` (`mu2e/ots_ops`), plus optional
`ValCrvDigi` inclusive plots (`BarId`, `SiPM`, `ADC`) gated by
`Config::fillInclusive`.

Typical use from an art module:

```cpp
CRVDigiDQM dqm(config);          // constructor / beginJob
dqm.Book(tfs->mkdir("CRVDigiDQM"));
dqm.Fill(*digis, *status);       // analyze, once per event
dqm.WriteGraphs();               // endJob (TGraph is not auto-saved)
```

`CRVReco/src/CrvDQMcollector_module.cc` uses this helper for all per-event
digi histograms and the CRVId rate maps (`crvDigiRates_ROC*`, 2D
`crvDigiRates`, `crvDigisPerChannel`). The collector still owns reco/PE/
coincidence products. Per-sector `crvDigisPerChannelAndEvent_CRVsector*`
is filled in the collector (and optionally `CRVDigiDQMAnalyzer`) `endJob`
from the helper's offline-channel counts; GeometryService supplies sector
names. The helper itself has no GeometryService or Proditions dependency.

`Config::fillCrvIdRates` (default true) books the two rate maps and the
detector-wide 1D vs offline channel. `WriteGraphs()` scales those three by
`1/nEvents`. VST `h1_channels` / `h2_channels` are unchanged (folded 25-slot
map).

Event-window tags for the time-series plots come from
`CrvStatus::GetEventWindowTag()`. If the status collection is empty (typical
MC), occupancy / ADC / TDC histograms are still filled and the EWT graphs,
rolling occupancy, and MicroBunchStatus plots are skipped.

Reco pulses and coincidences are out of scope; those belong in a future
`CRVRecoDQM` helper in this package.

## CRVStatusDQM

`mu2e::CRVStatusDQM` books and fills ROC-firmware health histograms from
`CrvStatus` / `CrvDAQerror`. It does **not** fold `roc==4` onto `roc==2`
(status is per DTC link). The required online plot is TH2 `errorBitsVsRoc`
(firmware bits 24–31 vs `dtcId*6+linkId`).

Typical use:

```cpp
CRVStatusDQM dqm(config);
dqm.Book(tfs->mkdir("CRVStatusDQM"));
dqm.Fill(*status, *daqErrors);   // analyze
dqm.EndSubRun(run, subrun);      // endSubRun
dqm.WriteGraphs();               // endJob
```

Empty `CrvStatus` (typical MC) still counts `nEvents` and skips ROC/latency
fills. `lastEventRocs()` supplies the five artdaq LastPoint scalars
(`TriggerCount`, `EventWindowTag`, `ActiveFEBCount`, `MicroBunchStatus`,
`WordCount`) with names `CRV.DTC<n>.ROC<m>.*`.

```text
mu2e -c Offline/DQMHelpers/fcl/CRVStatusDQM.fcl -s <digi.art>
```

A missing `CrvStatus` product is tolerated (empty collection). A missing
`CrvDAQerror` product skips unpack-error histograms.

## Offline analyzer

`CRVDigiDQMAnalyzer` is a thin wrapper:

```text
mu2e -c Offline/DQMHelpers/fcl/CRVDigiDQM.fcl -s <digi.art>
```

A missing `CrvStatus` product is tolerated (empty collection). A missing
`CrvDigi` product skips the event.

`fillSectorOccupancy: false` by default so the stock FCL needs no CRV
geometry. Set it true (as in `dqm_data_vst.fcl`) to book
`crvDigisPerChannelAndEvent_CRVsector*` from CosmicRayShield sector names
and fill them in `endJob`. The collector still skips `notConnected`
channels via Proditions; the analyzer does not.

## otsdaq CrvDQM

`otsdaq-mu2e-crv` `CrvDQM_module.cc` (`feature/CRVDigiDQM`, based on
`mu2e/ots_ops`) constructs a `mu2e::CRVDigiDQM` member with
`Config::fillInclusive = false`, calls `Book`/`Fill`/`WriteGraphs`, and
links `Offline::DQMHelpers`. The module still owns `h1_dummy`,
`HistoSender`, `THttpServer`/`TCanvas`/`CrvDQMStyle`, rate-log counters,
and the per-FEB timing summary canvases.

The DAQ build must pick up an Offline that contains `DQMHelpers` (the
`off_dqm` checkout, or a later Offline tag). Constant-fraction timing
lives in `Offline/DQMHelpers/inc/CRVCFTime.hh`; the local
`ArtModules/CrvCFTime.hh` is unused by `CrvDQM`.

Binning FHiCL defaults in `CRVDigiDQM::Config` match the online
`ps.get(...)` defaults so a cutover on the same file is comparable.

## otsdaq CrvStatusMetrics

`otsdaq-mu2e-crv` `CrvStatusMetrics_module.cc` constructs a
`mu2e::CRVStatusDQM` member, prefers `CrvStatus`/`CrvDAQerror` products,
and still publishes the five LastPoint series. If the status product is
absent it falls back to DTC-fragment decode for LastPoint only.
`HistoSender` ships `errorBitsVsRoc` (and the other status hists) when
`sendHists: true`. File-mode FCL: `fcl/RunCrvStatusDQM_vst_raw.fcl`.

