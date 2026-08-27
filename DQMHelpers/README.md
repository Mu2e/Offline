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
`1/nEvents`.

### Readout geography (`kppReadout`)

The helper carries two channel-ID conventions, and `Config::kppReadout`
(default true) picks which one is live. They are mutually exclusive — the
detector is either cabled the KPP way or it is not.

`kppReadout: true` is the KPP cabling.
DTC link 3 is read out as ROC 4 and folded onto
ROC 2, FEBs are numbered `(roc-1)*25 + feb`, and `h1_channels` /
`h1_channelsLastEwt` / `h2_channels` are booked over the resulting 33 FEB
slots.

`kppReadout: false` is the full CRV, which does not exist yet. The ROC fold is turned off
and the occupancy trio is not booked rather than resized, since 394 of 432 FEBs
would fall past the 2112-bin axis and `crvDigisPerChannel` / `crvDigiRates`
already cover the full detector correctly binned.

`h1_channelsLastEwt`, the rolling EWT-window occupancy, has no full-CRV
equivalent and would need rebasing onto the CRVId convention. Irrelevant while
KPP is the only CRV that exists.

`dtVsFeb` is booked in both modes; its x-axis follows the same geography
(`nFebIdBins()`) so it stays legible rather than reserving 450 bins for FEBs
that do not exist.

### Inter-FEB sync (`dtVsFeb`)

FEBs drifting out of sync with each other is the error case these timing plots
exist to catch. Because a slip is a per-FEB offset, every pairwise difference is
just `dt(i,j) = d_j - d_i` — the N x N pair matrix carries only N independent
numbers. `dtVsFeb` therefore stores those N numbers directly: x is the full
`globalFebId`, y is that FEB's first constant-fraction hit time minus the median
of the *other* FEBs' first hit times in the same event. A slipped FEB is a
displaced vertical stripe; a whole ROC slipping is a block of adjacent stripes,
since `globalFebId` is ROC-ordered.

The reference excludes the FEB being filled so that one bad FEB does not drag the
reference and smear its partners. Events with fewer than two FEBs are skipped
(no reference exists), and with exactly two FEBs both entries are `+/-` the pair
difference — attribution then comes from aggregating over events, since the bad
FEB is displaced against *every* partner while each partner is displaced only in
the events it shares with the bad one.

### Axis-coverage diagnostics

Both axes can hide entries, so the raw values are counted before `Fill` and
exposed as accessors — `maxFebIdSeen()`, `nFebIdOutOfAxis()`, `nDtOutOfRange()`,
`maxAbsDtSeen()`. Reading them beats reading an overflow bin, since
`nDtOutOfRange` and `maxAbsDtSeen` are measured on the true `dt` and so register
a slip of any size. Per-FEB attribution comes from the histograms below.

**Off-axis FEB — one warning per job.** The first time a digi arrives from a FEB
outside the axis, one `mf::LogWarning` names the ROC, FEB and `globalFebId`, then
latches. The axis size comes from `Config::kppReadout`, so this is a fixed
configuration error: it cannot change between runs, and repeating it would add
nothing. The helper needs no run-boundary hook. This is the only thing it logs.

### Per-FEB desync counters

An off-scale `dt` is a physics observation, not a misconfiguration, so it is
counted rather than logged. Two `TH1F`s, both on the same `globalFebId` axis as
`dtVsFeb`, counting events in which that FEB had `abs(dt) > dtVsFebRange`:

| Histogram | Period |
|---|---|
| `dtOutOfRangePerFeb` | whole job / file — the offline view |
| `dtOutOfRangePerFebLastEwt` | rolling `channelsWindowEwts` — the online view |

A FEB whose clock has slipped shows a bar standing above its neighbours. The
rolling twin exists so a slip that starts now is not diluted by hours of earlier
good data; it shares the deque-and-`AddBinContent` mechanism used by
`h1_channelsLastEwt`, so like it it carries bin content without entry counts and
only updates on events that carry an EWT (i.e. not MC).

A FEB contributes at most one entry per event, and only when it has a digi
passing `cfTime` (>= 3 samples, `peak - adcs[0] >= minAmplitude`, a leading-edge
crossing) *and* at least one other FEB in the event does too — a lone FEB has no
reference. So the count carries an occupancy term, which matters only while the
healthy baseline is non-zero; see the caveat under `dtVsFebRange` above.

This replaces the former per-FEB-pair `timing_feb/dt_febXX_febYY` histograms.
Intra-FEB `timing_fpga/` histograms are unchanged.

Event-window tags for the time-series plots come from
`CrvStatus::GetEventWindowTag()`. If the status collection is empty (typical
MC), occupancy / ADC / TDC histograms are still filled and the EWT graphs,
rolling occupancy, and MicroBunchStatus plots are skipped.

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
geometry. Set it true to book
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

