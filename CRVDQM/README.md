# CRVDQM

The CRV clients of the [`DQMHelpers`](../DQMHelpers/README.md) core: the CRV
histogram sets, their fixed Run 1 binning, and the code that fills them. Read
the `DQMHelpers` README first — the three rules there (binning is code, the
histogram set is code, changing either is a version bump) are what this
package implements for the CRV.

```
CRVDQM/
  inc/CRVDQMRun1.hh       frozen CRV layout: numbering, configurations (constants only)
  inc/CRVDQMLayout.hh     geometry/conditions -> what a CRV client needs
  inc/CRVDigiDQM.hh       digis: occupancy, ADC, TDC, timing
  inc/CRVDigiDQMConfig.hh FHiCL of CRVDigiDQM: hists + the configured dtPairs
  inc/CRVStatusDQM.hh     ROC status packets and DAQ errors
  inc/CRVRecoDQM.hh       reco pulses and coincidence clusters, PE/MPV maps
  inc/CRVCFTime.hh        constant-fraction timing of a digi waveform
  fcl/prolog.fcl          CRVDQM.* client tables, offline and online
  data/CRVDQM_binning_v<N>_<configuration>.txt   reference catalogues, one per configuration
```

The art modules live outside Offline: offline ones in
[Mu2e/DQM](https://github.com/Mu2e/DQM) (`DqmCrvDigi`, `DqmCrvStatus`,
`DqmCrvReco`), online ones (`CrvDQM`, `CrvStatusMetrics`) in
[Mu2e/otsdaq-mu2e-dqm](https://github.com/Mu2e/otsdaq-mu2e-dqm), alongside the
calorimeter and tracker online DQM. Both link `Offline::CRVDQM`.

## The frozen layout

`CRVDQMRun1.hh` holds the frozen layout: dense CRVId numbering (`febId` =
`(ROC-1)*24+(FEB-1)`, 432 ports; `onlineChannel` = `febId*64+FEBchannel`,
27 648), the constant-fraction timing constants (`kCFFraction` 0.20,
`kCFMinAmplitude` 10 ADC), the partner-timing selection (`kDtMinAmplitude`
200 ADC, `kDtCoincWindow` 20 ns, `kDtMinLayers` 3), and the two supported
configurations:

| Configuration | Geometry `crs.name` | Sectors | Readout |
|---|---|---|---|
| `run1a` | `run1a_v01` | T1, T2 | 18 ROCs, 432 FEB ports |
| `extracted` (KPP) | `extracted` (v02, v03), `extracted_v04` | EX, T1, T2, M1–M8 | 2 ROCs; ROC 1 ports 1–24, ROC 2 ports 2–7 (global FEB ID 0–30). ROC 2 is physical ROC 4 on DTC link 3 (unpacker `useROC4asROC2`), so DTC 0 links 0–3 are booked for status |

A client is constructed for one configuration (`CRVDigiDQM(hists, cfg)`, the
DQMHelpers *variant*) and books for it alone: the global-FEB, online-channel,
FPGA-pair and DTC-link axes are sized by the layout's `nROC`/`nFebs`
(`CRVDQMRun1::globalFebAxis(cfg)` etc.), the offline-channel axis of
`crvDigisPerChannel` and `crvPEsVsChannel` by `nOfflineChannels` (4 × the bar
count of the configuration's largest geometry: 5120 for `run1a_v01`, 2304 for
`extracted_v04`), the per-ROC and per-link sets run to
`nROC`, and the per-sector occupancy, per-sector PE MPV and cluster position
axes `X/Y/Z_<configuration>` (the counter envelope ± 5% of its span, 100
bins) exist for that configuration only. The output is stamped
`variant=<configuration>`. A digi or status block from a ROC, port or link the
configuration does not read out, or with an offline channel past its axis, is
counted in the diagnostics (`onlineIdOutsideConfiguration`, `linkOutOfRange`,
`offlineChannelOutsideConfiguration`; a geometry with more channels than that
axis, `geometryChannelsOutsideConfiguration`) and left out; a geometry
that names a different configuration at run start is refused by
`SetConfiguration` (`configurationMismatch`) rather than filled into the wrong
axes. The online modules take the configuration from FHiCL
(`configuration : "extracted"`, their default), the offline modules default to
`run1a`.

## The clients

| Client | Fills | Needs from the module |
|---|---|---|
| `CRVDigiDQM` | occupancy per channel/FEB, peak ADC, TDC, intra-FEB `dtFpgaPairs`, partner-FEB `dtPartner_*`, coincidence multiplicity, configured `dtPair_*` | `Fill(digis, status)`; `SetConfiguration`, `SetFebTopology`; `dtPairs` at construction |
| `CRVStatusDQM` | ROC headers, latency (and the invalid-latency counts), firmware error bits, port flags, per-link views, DAQ error codes | `Fill(status[, daqErrors])` |
| `CRVRecoDQM` | clusters, reco pulses, `crvPEsVsChannel`, sector-pair `crvCoincidenceDt_*`, and at end of job the MPV maps | `Fill(clusters[, pulses])`; `SetConfiguration` |

`CRVDQMLayout.hh` turns geometry and conditions into those injections, and is
shared by the offline and online modules:

```cpp
const int cfg = CRVDQMLayout::configuration(*crs);      // throws on unknown geometry
dqm.SetConfiguration(cfg, CRVDQMLayout::channelToSector(*crs, cfg, &sipmStatus));
CRVDQMLayout::febTopology(*crs, channelMap, topology, channelToLayer);
dqm.SetFebTopology(topology, channelToLayer);
```

A geometry that is not in the table, or a sector with no histograms, throws:
that is a layout change and needs a version bump, not a silent empty plot.
Redo the calls when the run changes, since the channel map can. The
configuration passed here must be the one the client was constructed for.

**Partner-FEB timing.** One FEB reads two of the four layers on one side of one
module, so a muon fires a handful of geometrically related FEBs. dt is measured
only inside a local coincidence group (≥ 3 of a module group's 4 layers within
20 ns, the group being one module or two adjacent ones in a sector), never
against the detector as a whole, and never across sectors — a muon entering one
side and leaving the other gives two genuinely separated traversals. Each
ordered partner pair fills one of four `dtPartner_<class>` maps (same module
same side / other side, adjacent module same side / other side), x = global FEB ID.
A slipped FEB is one displaced column. `layersPerGroup`, `groupsPerEvent`,
`sectorsPerEvent` and `febNoGroup` make the selection auditable; a large
multi-sector event is a cosmic air shower, which is physics, not a fault.

**Configured FPGA pairs.** `dtPairs` in the digi table names a short list
of (ROC, FEB, FPGA) pairs A, B; each books one TH1F `dtPair_ROC<r>FEB<f>FPGA<p>_ROC<r>FEB<f>FPGA<p>`
of t_B − t_A over every digi pair on the two FPGAs in an event (±250 ns,
0.5 ns bins). No selection is applied to the digis: the time is the
constant-fraction time with no amplitude cut, or, when the waveform has no
crossing, the time of its largest sample (counted as `dtPairPeakTime`). A
and B are any two FPGAs of the readout. A = B pairs one FPGA's different
channels, as `dtFpgaPairs` does, and then t_B − t_A is the later digi in the
collection minus the earlier. The times are uncalibrated, so a pair on two
FEBs sits at their relative delay. These are the one part of the set chosen
in FHiCL: a non-empty list is appended, sorted, to the stamped variant
(`variant=extracted+dtPairs:ROC1FEB3FPGA0_ROC1FEB3FPGA1,...`), so files with
different lists never merge. The reference catalogues are written with no
pairs.

**Sector-pair timing.** `crvCoincidenceDt_<configuration>_<A>_<B>`, one per
unordered pair of the configuration's sectors (A before B in its list:
55 for `extracted`, 1 for `run1a`), holds t_B − t_A of the average hit times
of two clusters of one event in those sectors (±100 ns, 0.5 ns bins). A
cluster's sector comes from its pulses' channels through the injected sector
map. The coincidence finder groups sectors by type (M1-M8 share one, as do
Run 1A T1 and T2), so a cluster can have pulses in two sectors; its time mixes
them, and it is left out and counted in the `clusterSpansSectors` diagnostic.
In KPP run 124155 that is 407 of 37,307 clusters, all on the muon taggers
(17% of their 2445 clusters); the per-tagger peaks move by under 3 ns. Cluster times are built from calibrated pulse times, so a pair spanning
two ROCs is the ROC-to-ROC comparison with the calibration applied; the peak
also holds the flight time between the sectors.

**DAQ errors.** `daqErrorCode` has one bin per `CrvDAQerror` code with room
for codes added upstream (32 bins), labelled from
`CrvDAQerrorCodeDetail::names()`, so a new code appears without a binning
change. `eventHasDaqError` flags an event with any of them.

## FHiCL

The client tables are in `fcl/prolog.fcl`, which includes the generic presets:
`CRVDQM.{Digi,Status,Reco}` offline, `CRVDQM.{DigiOnline,StatusOnline,
StatusLinksOnline}` online. They are configuration-agnostic, so the same tables
serve Run 1A/1B and KPP.

```
#include "Offline/CRVDQM/fcl/prolog.fcl"
crvDigiDQM : { module_type : DqmCrvDigi  dqm : @local::CRVDQM.Digi }
physics.analyzers.crvDigiDQM.dqm.hists : @local::DQMHelpers.Hists.PerSubrun
```

## Changing binning

Follow the procedure in the `DQMHelpers` README. The CRV reference catalogues
are `data/CRVDQM_binning_v<N>_<configuration>.txt`, one per configuration:
two comment lines, then the dump that a module writes when its
`catalogueFile` is set, run against that configuration's geometry. Version 2
introduced the per-configuration readout axes; `v1` is the pre-variant set
with every configuration's family booked. Version 3 (`CRVDigiDQM`,
`CRVRecoDQM`; `CRVStatusDQM` stays at 2) sizes the offline-channel axis of
`crvDigisPerChannel` and `crvPEsVsChannel` by configuration. Version 4
(`CRVRecoDQM` only) adds the `crvCoincidenceDt_*` sector pairs.

## Testing

The harness lives in the analysis area (`analysis/crv_efficiency`), not in
Offline:

- `scripts/dqm_refactor_phase0.sh` — steps `v1` (clients), `v3` (FHiCL
  tables), `v4` (DQM-repo modules, hadd);
- `scripts/dqm_phase2_compare.py` — compare two DQM files histogram by
  histogram (`--exact` when both use the current numbering);
- `scripts/dqm_hadd_check.py` — `hadd(A, B) == A + B`, same axes, stamps agree.

The plan and every decision behind the DQMHelpers/CRVDQM split are in
`docs/crv_dqm_helpers_refactor_plan.md` of the analysis workspace.
