# DQMHelpers

Shared DQM histogram helpers used by both offline art analyzers and the
otsdaq online monitor. The helper classes own histogram booking and filling;
art modules own I/O, visualization, and histogram shipping.

This package has no GeometryService or ProditionsService dependency, so the
same fill path can run in the DAQ process.

## Output file layout

`TFileService` writes **one file per art job** and gives **each module its own
directory**, named after the module label. That directory is automatic and not
optional — it is what keeps two modules from colliding in one file.

`outputTag` adds a *second* level inside the module directory. Use it only when
the module books histograms of its own alongside the helper's, so the two owners
stay separated:

| Module | `outputTag` | Why |
|---|---|---|
| `CRVStatusDQMAnalyzer` | `""` | books nothing itself — a second level would just repeat the module label |
| `CRVDigiDQMAnalyzer` | `""` | as above; set it if `fillSectorOccupancy` is on and you want the per-sector hists kept apart |
| `CrvDQMcollector` | `"CRVDigiDQM"` | **keep** — the module also books `crvPEsMPV*`, `crvPedestals*`, `crvCalibConstants*` and the `crvMetaData` tree at module level |

With `CRVDQM.fcl` (module labels `CRVDigiDQM` / `CRVStatusDQM`, empty `outputTag`):

```text
crvDQM.root
  CRVDigiDQM/
    nEvents
    h1_digisPerEvt, h1_peakAdc, h1_tdc
    h1_channels, h2_channels
    crvDigiRates, crvDigiRates_ROC*, crvDigisPerChannel
    dtVsFeb, dtOutOfRangePerFeb
    timing_fpga/
    BarId, SiPM, ADC                  (fillInclusive)
  CRVStatusDQM/
    nEvents
    nRocHeaders, errorBits, errorBitsVsRoc
    linkLatency, linkLatency_dtc*_roc*
    rocCensus, eventHasError, eventHasDaqError
    daqErrorCode, ewtMismatch
    errorsPerSubrun, meanLatencyPerSubrun
    ...
```

`CrvDQMcollector` keeps PE / pedestal / calib histograms at module level, so the
helper's occupancy maps live one directory deeper:

```text
<collector module label>/
  crvPEsMPV*, crvPedestals*, crvCalibConstants*, crvMetaData
  CRVDigiDQM/
    nEvents
    crvDigiRates, crvDigiRates_ROC*, crvDigisPerChannel
    crvDigisPerChannelAndEvent_CRVsector*   (if BookSectorOccupancy was called)
    ...
```

### Combining per-job files

`hadd` of per-job DQM files is the intended combine step (run, week, …). Occupancy
maps are stored as **raw counts**. After `hadd`,

```text
rate = content / nEvents->GetBinContent(1)
```

TGraphs vs EWT or subrun, and `*LastEwt` rolling snapshots, are live-monitor
objects (`fillLivePlots: true` in the online modules). They are not booked by
the offline analyzers or the collector, so they are not in the files that get
`hadd`'d.

### Per-sector occupancy

`crvDigisPerChannelAndEvent_CRVsector*` is owned by the helper. Sector names and the
channel->sector map are geometry-derived, so the caller injects them via
`BookSectorOccupancy(...)` and the helper stays free of GeometryService. A
negative sector entry skips that channel, which is how `CrvDQMcollector` drops
its Proditions `notConnected` channels without the helper knowing about
Proditions. `WriteGraphs()` fills them once per file from `nDigisOffline` as a
per-channel rate distribution. After `hadd`, rebuild that distribution from
`crvDigisPerChannel` / `nEvents` rather than using the combined hist.

An empty `outputTag` books directly in the module directory. The module labels in
`CRVDQM.fcl` are therefore `CRVDigiDQM` / `CRVStatusDQM`, so the directory name is
the same either way.

`CRVDQM.fcl` is the only FCL in this package: digi and status always run together
in one job, into one file. Running one side alone is a matter of dropping the
other module from `physics.ana`, not a separate config to keep in step.

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
dqm.WriteGraphs();               // endJob (sector occupancy + persist TGraphs if booked)
```

`CRVReco/src/CrvDQMcollector_module.cc` uses this helper for all per-event
digi histograms and the CRVId occupancy maps (`crvDigiRates_ROC*`, 2D
`crvDigiRates`, `crvDigisPerChannel`). The collector still owns reco/PE/
coincidence products. Per-sector `crvDigisPerChannelAndEvent_CRVsector*` is
owned by the helper; the modules only inject the geometry-derived sector names
and channel->sector map. The helper itself has no GeometryService or Proditions
dependency.

`Config::fillCrvIdRates` (default true) books the two occupancy maps and the
detector-wide 1D vs offline channel. Those three are raw counts; divide by
`nEvents` after `hadd`.

### Readout geography (`kppReadout`)

`Config::kppReadout` (default true) sizes the FEB axes. It selects detector
size only — it does **no** ROC remapping.

`kppReadout: true` sizes for KPP, the extracted CRV and the only one built so
far: ROC 1-2, FEBs numbered `(roc-1)*25 + feb`, and `h1_channels` /
`h2_channels` booked over 33 FEB slots (`h1_channelsLastEwt` only when
`fillLivePlots` is true). Every existing dataset wants this mode, which is why
it is the default here and in `prolog_v12.fcl`.

`kppReadout: false` is the full CRV, which does not exist yet — a seam for when
it does, not a mode anything runs in today. The occupancy trio is **not booked**
rather than resized, since 394 of 432 FEBs would fall past the 2112-bin axis and
`crvDigisPerChannel` / `crvDigiRates` already cover the full detector correctly
binned.

**No ROC fold here.** DTC link 3 is read out as ROC 4 and folded onto ROC 2 by
`CrvDigisFromArtdaqFragmentsFEBII` (`useROC4asROC2`) at unpack time, so digis
reach the helper with the correct ROC already. Re-folding would duplicate a
DAQ-level mapping and silently mask an unpacker configured without it; instead a
ROC 4 arriving here lands off-axis and raises the warning below.

`h1_channelsLastEwt`, the rolling EWT-window occupancy, is booked only with
`fillLivePlots` and has no full-CRV equivalent. Irrelevant while KPP is the
only CRV that exists.

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
exposed as accessors — `maxFebIdSeen()`, `nFebIdOutOfAxis()`, `nCrvIdOutOfRange()`,
`nDtOutOfRange()`, `maxAbsDtSeen()`. Reading them beats reading an overflow bin,
since `nDtOutOfRange` and `maxAbsDtSeen` are measured on the true `dt` and so
register a slip of any size. Per-FEB attribution comes from the histograms below.

**Off-axis FEB — one warning per job.** The first time a digi arrives from a FEB
outside the axis, one `mf::LogWarning` names the ROC, FEB and `globalFebId`, then
latches. The axis size comes from `Config::kppReadout`, so this is a fixed
configuration error: it cannot change between runs, and repeating it would add
nothing. The helper needs no run-boundary hook.

**CRVId-range ROC/FEB — one warning per job.** Occupancy maps (`crvDigiRates*`) skip a
digi whose ROC/FEB/channel is outside `CRVId`. Occupancy and timing still fill.
`nCrvIdOutOfRange()` counts the skips. This is not a hard failure: the helper
runs in the DAQ process.

This is the only logging the helper does.

### Per-FEB desync counters

An off-scale `dt` is a physics observation, not a misconfiguration, so it is
counted rather than logged. `dtOutOfRangePerFeb` is on the same `globalFebId`
axis as `dtVsFeb`, counting events in which that FEB had `abs(dt) > dtVsFebRange`
over the whole job. With `fillLivePlots`, `dtOutOfRangePerFebLastEwt` is the
rolling `channelsWindowEwts` twin for the online monitor.

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
`CrvStatus::GetEventWindowTag()`. Those plots are booked only when
`fillLivePlots` is true. If the status collection is empty (typical MC),
occupancy / ADC / TDC histograms are still filled and the EWT graphs, rolling
occupancy, and MicroBunchStatus plots are skipped.

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

Empty `CrvStatus` (typical MC, product present but empty) still counts `nEvents`
and skips ROC/latency fills. A **missing** product (wrong tag) is the same fill
path plus one `LogWarning` per job from the analyzer — it is not a throw, so MC
jobs that omit status still run. `lastEventRocs()` supplies the five artdaq
LastPoint scalars
(`TriggerCount`, `EventWindowTag`, `ActiveFEBCount`, `MicroBunchStatus`,
`WordCount`) with names `CRV.DTC<n>.ROC<m>.*`.

Per-link latency histograms are booked only for `(dtcId, linkId)` that fall on
the ROC axis (`linkId < nROCPerDTC` and `dtcId*nROCPerDTC+linkId < nROC`).
Corrupt IDs still fill the inclusive `linkLatency` histogram. Latency is a
DTC-link field, so it is filled whether or not `HasROCHeader()` is true.

```text
mu2e -c Offline/DQMHelpers/fcl/CRVDQM.fcl -s <digi.art>
```

A missing `CrvStatus` product is tolerated (empty collection, one warning per
job). A missing `CrvDAQerror` product skips unpack-error histograms (one
warning per job).

## Offline analyzer

`CRVDigiDQMAnalyzer` is a thin wrapper:

```text
mu2e -c Offline/DQMHelpers/fcl/CRVDQM.fcl -s <digi.art>
```

A missing `CrvStatus` product is tolerated (empty collection, one warning per
job). A missing `CrvDigi` product skips the event (one warning per job).

`fillSectorOccupancy: false` by default so the stock FCL needs no CRV
geometry. Set it true to book
`crvDigisPerChannelAndEvent_CRVsector*` from CosmicRayShield sector names
and fill them in `endJob`. The collector still skips `notConnected`
channels via Proditions; the analyzer does not.

## otsdaq CrvDQM

`otsdaq-mu2e-crv` `CrvDQM_module.cc` (`feature/CRVDigiDQM`, based on
`mu2e/ots_ops`) constructs a `mu2e::CRVDigiDQM` member with
`Config::fillInclusive = false` and `Config::fillLivePlots = true`, calls
`Book`/`Fill`/`WriteGraphs`, and links `Offline::DQMHelpers`. The module still
owns `h1_dummy`, `HistoSender`, `THttpServer`/`TCanvas`/`CrvDQMStyle`, rate-log
counters, and the per-FEB timing summary canvases.

The DAQ build must pick up an Offline that contains `DQMHelpers` (the
`off_dqm` checkout, or a later Offline tag). Constant-fraction timing
lives in `Offline/DQMHelpers/inc/CRVCFTime.hh`; the local
`ArtModules/CrvCFTime.hh` is unused by `CrvDQM`.

Binning FHiCL defaults in `CRVDigiDQM::Config` match the online
`ps.get(...)` defaults so a cutover on the same file is comparable.

## otsdaq CrvStatusMetrics

`otsdaq-mu2e-crv` `CrvStatusMetrics_module.cc` constructs a
`mu2e::CRVStatusDQM` member with `fillLivePlots = true`, prefers
`CrvStatus`/`CrvDAQerror` products, and still publishes the five LastPoint
series. If the status product is absent it falls back to DTC-fragment decode
for LastPoint only. `HistoSender` ships `errorBitsVsRoc` (and the other status
hists) when `sendHists: true`. File-mode FCL: `fcl/RunCrvStatusDQM_vst_raw.fcl`.

