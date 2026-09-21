# CRVDQM

The CRV clients of the [`DQMHelpers`](../DQMHelpers/README.md) core: the CRV
histogram sets, their fixed Run 1 binning, and the code that fills them. Read
the `DQMHelpers` README first — the three rules there (binning is code, the
histogram set is code, changing either is a version bump) are what this
package implements for the CRV.

```
CRVDQM/
  inc/CRVDQMRun1.hh       frozen Run 1 CRV layout (constants only)
  inc/CRVDQMLayout.hh     geometry/conditions -> what a CRV client needs
  inc/CRVDigiDQM.hh       digis: occupancy, ADC, TDC, timing
  inc/CRVStatusDQM.hh     ROC status packets and DAQ errors
  inc/CRVRecoDQM.hh       reco pulses and coincidence clusters, PE/MPV maps
  inc/CRVCFTime.hh        constant-fraction timing of a digi waveform
  fcl/prolog.fcl          CRVDQM.* client tables, offline and online
  data/CRVDQM_binning_v1.txt   reference catalogue of the CRV histogram set
```

The art modules live outside Offline: offline ones in
[Mu2e/DQM](https://github.com/Mu2e/DQM) (`DqmCrvDigi`, `DqmCrvStatus`,
`DqmCrvReco`), online ones (`CrvDQM`, `CrvStatusMetrics`) in
[Mu2e/otsdaq-mu2e-dqm](https://github.com/Mu2e/otsdaq-mu2e-dqm), alongside the
calorimeter and tracker online DQM. Both link `Offline::CRVDQM`.

## The frozen layout

`CRVDQMRun1.hh` holds the frozen layout: dense CRVId numbering (`febPort` =
`(ROC-1)*24+(FEB-1)`, 432 ports; `onlineChannel` = `febPort*64+FEBchannel`,
27 648), the constant-fraction timing constants (`kCFFraction` 0.20,
`kCFMinAmplitude` 10 ADC), the partner-timing selection (`kDtMinAmplitude`
200 ADC, `kDtCoincWindow` 20 ns, `kDtMinLayers` 3), and the two supported
configurations:

| Configuration | Geometry `crs.name` | Sectors |
|---|---|---|
| `run1a` | `run1a_v01` | T1, T2 |
| `extracted` (KPP) | `extracted` (v02, v03), `extracted_v04` | EX, T1, T2, M1–M8 |

Anything whose meaning depends on the configuration is booked once per
configuration — per-sector occupancy, per-sector PE MPV, and the cluster
position axes `X/Y/Z_<configuration>` (the counter envelope ± 5% of its span,
100 bins). Every job books every family and fills the one its geometry names.

## The clients

| Client | Fills | Needs from the module |
|---|---|---|
| `CRVDigiDQM` | occupancy per channel/FEB, peak ADC, TDC, intra-FEB `dtFpgaPairs`, partner-FEB `dtPartner_*`, coincidence multiplicity | `Fill(digis, status)`; `SetConfiguration`, `SetFebTopology` |
| `CRVStatusDQM` | ROC headers, latency (and the invalid-latency counts), firmware error bits, port flags, per-link views, DAQ error codes | `Fill(status[, daqErrors])` |
| `CRVRecoDQM` | clusters, reco pulses, `crvPEsVsChannel`, and at end of job the MPV maps | `Fill(clusters[, pulses])`; `SetConfiguration` |

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
Redo the calls when the run changes, since the channel map can.

**Partner-FEB timing.** One FEB reads two of the four layers on one side of one
module, so a muon fires a handful of geometrically related FEBs. dt is measured
only inside a local coincidence group (≥ 3 of a module group's 4 layers within
20 ns, the group being one module or two adjacent ones in a sector), never
against the detector as a whole, and never across sectors — a muon entering one
side and leaving the other gives two genuinely separated traversals. Each
ordered partner pair fills one of four `dtPartner_<class>` maps (same module
same side / other side, adjacent module same side / other side), x = FEB port.
A slipped FEB is one displaced column. `layersPerGroup`, `groupsPerEvent`,
`sectorsPerEvent` and `febNoGroup` make the selection auditable; a large
multi-sector event is a cosmic air shower, which is physics, not a fault.

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

Follow the procedure in the `DQMHelpers` README. The CRV reference catalogue
is `data/CRVDQM_binning_v<N>.txt`: two comment lines, then the dump that a
module writes when its `catalogueFile` is set.

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
