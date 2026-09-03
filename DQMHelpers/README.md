# DQMHelpers

Shared DQM histogram helpers used by both Offline and the otsdaq online
monitor. The helper classes own histogram booking and filling; art modules own
I/O, visualization, and histogram shipping.

This package is helper **classes only** — no art plugins, no FHiCL, and so no
GeometryService or ProditionsService dependency anywhere in its build. That is
what lets the same fill path run in the DAQ process. Two consumers:

| Consumer | Uses | Where |
|---|---|---|
| `CrvDQMcollector` | all three helpers | `Offline/CRVReco/src/CrvDQMcollector_module.cc` |
| `CrvDQM` / `CrvStatusMetrics` | digi and status helpers | `otsdaq-mu2e-crv` |

Anything a helper needs from geometry or Proditions is **injected** by the
consumer — the sector map, the channel-to-sector map, the CRV envelope — so the
dependency stays on the caller's side of the boundary.

## Output file layout

`TFileService` writes **one file per art job** and gives **each module its own
directory**, named after the module label. That directory is automatic and not
optional — it is what keeps two modules from colliding in one file.

`CrvDQMcollector` books three helpers into one file, so two of them take a
subdirectory and the third keeps the module directory:

| Setting | Default | Why |
|---|---|---|
| `crvDigiDQMDir` | `"CRVDigiDQM"` | its own level, so it does not mix with what the module books |
| `crvStatusDQMDir` | `"CRVStatusDQM"` | as above |
| `crvRecoDQMDir` | `""` | **keep empty** — `crvPEsMPV*` and `crvCoincidencesClusters` have always sat in the module directory and the KPP extractor reads them there |

The module itself books `crvPedestals*`, `crvCalibConstants*` and the
`crvMetaData` tree at module level, alongside the reco helper's histograms.

```text
<collector module label>/
  crvPedestals*, crvCalibConstants*, crvMetaData        (module)
  nEvents, nEventsWithCoincidenceClusters               (CRVRecoDQM)
  crvCoincidencesClusters                               (CRVRecoDQM)
  crvPEsMPV, crvPEsMPV_ROC*, crvPEsMPV_CRVsector*       (CRVRecoDQM)
  NPulses, NPulse2, BarIdr, SiPMr, PEr, PEHeight,       (CRVRecoDQM,
  PulseTime, PulseTime2, chi2, logchi2, LeadingTime,     fillInclusive)
  LeadingTime2, NClus, NPc, PEc, tc, t2c, X, Y, Z
  CRVDigiDQM/
    nEvents
    h1_digisPerEvt, h1_peakAdc, h1_tdc
    h1_channels, h2_channels
    crvDigiRates, crvDigiRates_ROC*, crvDigisPerChannel
    crvDigisPerChannelAndEvent_CRVsector*   (if BookSectorOccupancy was called)
    dtVsFeb, dtOutOfRangePerFeb
    timing_fpga/
    BarId, SiPM, ADC                        (fillInclusive)
  CRVStatusDQM/
    nEvents
    nRocHeaders, errorBits, errorBitsVsRoc
    linkLatency, linkLatency_dtc*_roc*
    rocCensus, eventHasError, eventHasDaqError
    daqErrorCode, ewtMismatch
    errorsPerSubrun, meanLatencyPerSubrun
```

### Combining per-job files

`hadd` of per-job DQM files is the intended combine step (run, week, …). Occupancy
maps are stored as **raw counts**. After `hadd`,

```text
rate = content / nEvents->GetBinContent(1)
```

TGraphs vs EWT or subrun, and `*LastEwt` rolling snapshots, are live-monitor
objects (`fillLivePlots: true` in the online modules). The collector leaves
`fillLivePlots` false, so they are not in the files that get `hadd`'d.

### Per-sector occupancy

`crvDigisPerChannelAndEvent_CRVsector*` is owned by the helper. Sector names and the
channel->sector map are geometry-derived, so the caller injects them via
`BookSectorOccupancy(...)` and the helper stays free of GeometryService. A
negative sector entry skips that channel, which is how `CrvDQMcollector` drops
its Proditions `notConnected` channels without the helper knowing about
Proditions. `WriteGraphs()` fills them once per file from `nDigisOffline` as a
per-channel rate distribution. After `hadd`, rebuild that distribution from
`crvDigisPerChannel` / `nEvents` rather than using the combined hist.

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
`crvDigiRates`, `crvDigisPerChannel`); reco pulses and coincidences go through
`CRVRecoDQM` below, and the collector keeps only the Proditions-derived
pedestal / calibration histograms and the metadata tree. Per-sector
`crvDigisPerChannelAndEvent_CRVsector*` is
owned by the helper; the modules only inject the geometry-derived sector names
and channel->sector map. The helper itself has no GeometryService or Proditions
dependency.

`Config::fillCrvIdRates` (default true) books the two occupancy maps and the
detector-wide 1D vs offline channel. Those three are raw counts; divide by
`nEvents` after `hadd`.

### Readout geography (`kppReadout`)

`Config::kppReadout` (default true) sizes the FEB axes. It selects detector
size only — it does no ROC remapping.

`kppReadout: true` sizes for KPP, the extracted CRV and the only one built so
far: ROC 1-2, FEBs numbered `(roc-1)*25 + feb`, and `h1_channels` /
`h2_channels` booked over 33 FEB slots (`h1_channelsLastEwt` only when
`fillLivePlots` is true). Every existing dataset wants this mode, which is why
it is the default here and in `prolog_v12.fcl`.

`kppReadout: false` is the full CRV, which does not exist yet — a seam for when
it does, not a mode anything runs in today. The occupancy trio is not booked
rather than resized, since 394 of 432 FEBs would fall past the 2112-bin axis and
`crvDigisPerChannel` / `crvDigiRates` already cover the full detector correctly
binned.

`h1_channelsLastEwt`, the rolling EWT-window occupancy, is booked only with
`fillLivePlots` and has no full-CRV equivalent.

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
`nCrvIdOutOfRange()` counts the skips.

### Per-FEB desync counters

`dtOutOfRangePerFeb` is on the same `globalFebId` axis as `dtVsFeb`,
counting events in which that FEB had `abs(dt) > dtVsFebRange` over the whole job.
With `fillLivePlots`, `dtOutOfRangePerFebLastEwt` is the
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

## CRVRecoDQM

`mu2e::CRVRecoDQM` books and fills reco-pulse and coincidence histograms: 
`crvCoincidencesClusters`, `crvPEsMPV_CRVsector*`, `crvPEsMPV_ROC*`, and the 2D `crvPEsMPV`.

```cpp
CRVRecoDQM dqm(config);              // constructor / beginJob
dqm.Book(tfs->mkdir("CRVRecoDQM"));  // or *tfs, to book in the module directory
dqm.BookSectorMPV(sectorNames, channelToSector);   // once, optional
dqm.Fill(*clusters, *recoPulses);    // analyze, once per event
dqm.WriteGraphs();                   // endJob: runs the fits, fills the MPV maps
```

Only pulses that belong to a coincidence cluster enter the PE spectra and the
MPV maps. The inclusive plots below see every pulse,
which is why they take the second `Fill` argument.

### Inclusive per-event plots (`fillInclusive`, default true)

```text
reco pulses (whole CrvRecoPulseCollection)
  NPulses, NPulse2, BarIdr, SiPMr, PEr, PEHeight,
  PulseTime, PulseTime2, chi2, logchi2, LeadingTime, LeadingTime2
coincidence clusters
  NClus, NPc, PEc, tc, t2c, X, Y, Z
```

The axes that follow the *detector* are configurable:

| Config | Applies to | Default |
|---|---|---|
| `nBinsTime`, `minTime`, `maxTime` | `PulseTime`, `LeadingTime`, `tc` | 100, 0, 2000 ns |
| `nBinsTime2`, `minTime2`, `maxTime2` | `PulseTime2`, `LeadingTime2`, `t2c` | 100, 0, 100 µs |
| `nBinsPos`, `minX`/`maxX`, `minY`/`maxY`, `minZ`/`maxZ` | `X`, `Y`, `Z` | 100, the `DqmCrv` full-CRV ranges |

Each time quantity is booked twice — a short view and a full-window view. The
three short axes share one definition and the three long axes another, so pulse,
leading-edge and cluster times stay overlayable; `X`/`Y`/`Z` share a bin count
but not a range.

The `X`/`Y`/`Z` entries are the **fallback only** — a caller with geometry
should inject the envelope instead, as `CrvDQMcollector` does; see *Cluster
position axes* below.

`NClus` gets an entry on **every** event, zero included, so it is a rate. The
reco-pulse block only fills through the two-argument `Fill`; a caller with no
`CrvRecoPulse` product still gets everything the clusters drive.

### No channel map

`CrvRecoPulse` carries `ROC` / `FEB` / `FEBchannel` copied from the digi it was
reconstructed from, and both `CrvDigitizer` (MC) and the artdaq unpacker (data)
fill those from `CRVOrdinal`. So the online-indexed plots need no Proditions. 
An ID outside the `CRVId` ranges is counted by `nOnlineIdOutOfRange()` and raises
one `mf::LogWarning`; the per-sector MPV still gets that pulse.

`crvPEsMPV_ROC*` and the 2D `crvPEsMPV` index the same online channel through
`onlineChannelIndex()` / `rocChannel()` / `febPort()`, and their axes are sized
from the same constants, so the index and the axis cannot drift apart.

### Per-sector MPV

`crvPEsMPV_CRVsector*` is owned by the helper. As with `BookSectorOccupancy`,
sector names and the channel->sector map are geometry-derived and are injected by
the caller, and a negative sector entry skips the channel — which is how
`CrvDQMcollector` drops its Proditions `notConnected` channels without the helper
knowing about Proditions. A caller with no Proditions keeps those channels;
they enter at MPV zero.

### The fit

`WriteGraphs()` runs the Landau(x)Gauss fit (`$ROOTSYS/tutorials/fit/langaus.C`)
on every channel spectrum and fills the MPV maps from it. A channel too sparse to
fit, or one whose fit lands on the range edge, enters as **zero** rather than
being skipped — a hole in the distribution at zero is what makes a dead or
dying channel visible. `nFits()`, `nFitsSucceeded()` and `meanFitChi2PerNdf()`
separate "MPV is zero" from "the fit never ran"; the analyzer prints them at
`diagLevel > 0`.

### Per-channel spectra

The `crvPEs_channel*` / `crvPEsROC_channel*` spectra are the fit input, not a DQM
product, and across all of `CRVId` there are ~22 k offline and ~27 k online of
them. They are booked **on first fill**, so a KPP-sized geometry pays for the few
thousand channels it actually reads out rather than for the whole detector — the
collector used to allocate every one of them in `beginRun`.

By default the helper owns them outright (`SetDirectory(nullptr)`) and nothing
writes them. Set `writePerChannelPE` to put the ones that were filled in a
`perChannelPE/` subdirectory, when a channel's MPV goes bad and someone needs to
see why.

### hadd

`nEvents`, `nEventsWithCoincidenceClusters` and `crvCoincidencesClusters` are raw
counts and add correctly. The MPV objects **do not**: `crvPEsMPV_ROC*` and the 2D
`crvPEsMPV` are bar charts whose bin content *is* the MPV, so `hadd` sums MPVs
rather than averaging them, and `crvPEsMPV_CRVsector*` combines two files' fits
of partial statistics rather than one fit of the sum. To combine runs, re-run the
job over the combined input, or refit from `perChannelPE/`.

## CRVStatusDQM

`mu2e::CRVStatusDQM` books and fills ROC-firmware health histograms from
`CrvStatus` / `CrvDAQerror`. 

Typical use:

```cpp
CRVStatusDQM dqm(config);
dqm.Book(tfs->mkdir("CRVStatusDQM"));
dqm.Fill(*status, *daqErrors);   // analyze
dqm.EndSubRun(run, subrun);      // endSubRun
dqm.WriteGraphs();               // endJob
```

Empty `CrvStatus` (typical MC, product present but empty) still counts `nEvents`
and skips ROC/latency fills. A missing product (wrong tag) is the same fill
path plus one `LogWarning` per job from the analyzer — it is not a throw, so MC
jobs that omit status still run. `lastEventRocs()` supplies the five artdaq
LastPoint scalars
(`TriggerCount`, `EventWindowTag`, `ActiveFEBCount`, `MicroBunchStatus`,
`WordCount`) with names `CRV.DTC<n>.ROC<m>.*`.

Per-link latency histograms are booked only for `(dtcId, linkId)` that fall on
the ROC axis (`linkId < nROCPerDTC` and `dtcId*nROCPerDTC+linkId < nROC`).
Corrupt IDs still fill the inclusive `linkLatency` histogram. Latency is a
DTC-link field, so it is filled whether or not `HasROCHeader()` is true.

A missing `CrvStatus` product is tolerated (empty collection, one warning per
job). A missing `CrvDAQerror` product skips unpack-error histograms (one
warning per job).

## Cluster position axes

`X`, `Y` and `Z` are the one part of the inclusive block whose right range
depends on where the detector physically is, so they are not booked by
`Book()`. A caller holding geometry injects the CRV envelope:

```cpp
dqm.BookPositionAxes(crvMin, crvMax);   // min/max corner, Mu2e coordinates
dqm.BookPositionAxes();                 // or: use the Config ranges
```

`CrvDQMcollector` builds the envelope in `beginRun` from the counters:

```cpp
for (const auto &bar : CRS->getAllCRSScintillatorBars())
  // crvMin[i] = min(crvMin[i], bar->getPosition()[i] - bar->getHalfLengths()[i]), etc.
```

On run 124155 the counter envelope gives x [-7204, -604], y [4203, 4743],
z [21116, 24887] with a 5 % margin, and every cluster lands in range. The
`Config` ranges remain as the no-geometry fallback; they are the full-CRV
`DqmCrv` values and are entirely overflow in y and z on extracted geometry.

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
`CrvStatus`/`CrvDAQerror` products, and publishes the five LastPoint
series. If the status product is absent it falls back to DTC-fragment decode
for LastPoint only. `HistoSender` ships `errorBitsVsRoc` (and the other status
hists) when `sendHists: true`. File-mode FCL: `fcl/RunCrvStatusDQM_vst_raw.fcl`.

