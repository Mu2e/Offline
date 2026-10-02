# DQMHelpers

Subdetector-agnostic histogramming for DQM. The classes here have no art
module, no `Event`, no service, no Proditions lookup and no detector
dependency: a module reads the event and calls a subdetector *client* built on
this core, so the same client runs in an offline art job and inside the otsdaq
online monitor.

```
DQMHelpers/
  inc/DQMAxis.hh          one axis, as a constant
  inc/DQMHistSet.hh       the histogram set: booking, copies, publishing
  inc/DQMHist.hh          DQMH1<T> / DQMH2<T>, the fill handles
  inc/DQMClient.hh        base class of a subdetector client
  inc/DQMSeries.hh        capped TGraph time series (online only)
  inc/DQMDiagnostics.hh   counters, warn-once, end-of-job summary
  inc/DQMHistSetConfig.hh validated FHiCL for the histogram set
  inc/DQMStyle.hh         ROOT styling, shared by online and offline displays
  fcl/prolog.fcl          generic presets
```

Subdetector clients live in their own packages and link `Offline::DQMHelpers`; nothing
detector-specific belongs here, including on the link line. The CRV clients
are in [`CRVDQM`](../CRVDQM/README.md), the worked example of everything
below.

## Three rules

1. **Binning is hard-coded.** An axis is a `static constexpr DQMAxis` in the client
   header. `book1`/`book2` take nothing else. No FHiCL path reaches a binning.
   This is to enable merging/comparison across runs.
2. **The histogram set is code.** A client books everything in `book()`, up
   front, for the whole detector. Booking later, or booking one name twice,
   throws. FHiCL chooses extra *copies* of a histogram and what is published.
3. **Changing either is a version bump.** Edit the constant, bump the client's
   `kBinningVersion`, and regenerate the catalogue file (see
   [Changing binning](#changing-binning)).

The point of all three is merging: files from different runs and different
jobs of the same client `hadd` bin for bin, and a metrics tool can refuse to
combine two binning versions because every client stamps
`dqmBinningVersion` into its output directory.

## Adding a subdetector

A subdetector's DQM goes in **its own Offline package** (`Offline/<Subdetector>DQM`),
a sibling of `Offline/DQMHelpers`, never inside it. `Offline/CRVDQM` is the complete, working
example, and every step below points at the CRV file to copy from. The steps
use a hypothetical tracker package, `Offline/TrkDQM`, laid out like
`Offline/CRVDQM`:

```
Offline/TrkDQM/                     compare Offline/CRVDQM/
  CMakeLists.txt                    Offline/CRVDQM/CMakeLists.txt
  inc/TrkDigiDQM.hh                 Offline/CRVDQM/inc/CRVDigiDQM.hh
  src/TrkDigiDQM.cc                 Offline/CRVDQM/src/CRVDigiDQM.cc
  src/SConscript                    Offline/CRVDQM/src/SConscript
  fcl/prolog.fcl                    Offline/CRVDQM/fcl/prolog.fcl
  data/TrkDQM_binning_v1.txt        Offline/CRVDQM/data/CRVDQM_binning_v1.txt
  README.md                         Offline/CRVDQM/README.md
```

**1. Write the client** in `Offline/TrkDQM/inc` and `Offline/TrkDQM/src`.
Axes are constants; `book()` books; `Fill()` starts with `beginEvent()`. This
is a whole client:

```cpp
// Offline/TrkDQM/inc/TrkDigiDQM.hh
class TrkDigiDQM : public DQMClient {
 public:
  static constexpr int kBinningVersion = 1;          // bump on any axis change
  static constexpr DQMAxis kNDigis = DQMAxis::Counts(0, 100);
  static constexpr DQMAxis kTdc{100, 0., 80.e3};

  explicit TrkDigiDQM(const DQMHistSet::Config& hists = {}) :
      DQMClient("TrkDigiDQM", kBinningVersion, hists) {}

  void Fill(const StrawDigiCollection& digis) {
    beginEvent();                        // counts the event, advances windows
    nDigis_.Fill(digis.size());
    for (const auto& d : digis) tdc_.Fill(d.TDC(StrawEnd::cal));
  }

 private:
  void book() override {
    nDigis_ = hists().book1<TH1F>("nDigis", "Digis / event;N;Events", kNDigis);
    tdc_    = hists().book1<TH1F>("tdc", "Digi TDC;TDC;Digis", kTdc);
  }
  DQMH1<TH1F> nDigis_, tdc_;
};
```

The CRV clients show the same pattern at full size:
- `Offline/CRVDQM/inc/CRVStatusDQM.hh` and `src/CRVStatusDQM.cc`: the simplest
  CRV client. It needs no layout from the module, only the status collection
  (and optionally the DAQ errors), and declares its axes at the top of the
  header.
- `Offline/CRVDQM/inc/CRVDigiDQM.hh` and `src/CRVDigiDQM.cc`: a client that
  also takes the detector layout from the module (`SetConfiguration`,
  `SetFebTopology`) and books one histogram family per configuration.
- `Offline/CRVDQM/inc/CRVRecoDQM.hh` and `src/CRVRecoDQM.cc`: end-of-job fit
  results booked with `bookSummary1`/`bookSummary2`.
- `Offline/CRVDQM/inc/CRVDQMRun1.hh`: frozen layout constants kept apart from
  the clients, so several clients share one numbering.

**2. Add it to the build**, in both build systems. `Offline/CRVDQM/CMakeLists.txt`
is the template; change the sources and the subdetector libraries:

```cmake
# Offline/CRVDQM/CMakeLists.txt (comments added here)
cet_make_library(
    SOURCE
      src/CRVDQMRun1.cc
      src/CRVDigiDQM.cc
      src/CRVRecoDQM.cc
      src/CRVStatusDQM.cc
    LIBRARIES PUBLIC
      Offline::DQMHelpers              # always
      Offline::CRVConditions           # the rest are what the CRV clients use
      Offline::CosmicRayShieldGeom
      Offline::RecoDataProducts
      Offline::DataProducts
      ROOT::Hist
      ROOT::MathCore
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CRVDQM/fcl)
```

The `install_fhicl` line is what makes
`#include "Offline/CRVDQM/fcl/prolog.fcl"` resolve in a job. For the tracker
it would read `SUBDIRNAME Offline/TrkDQM/fcl`.

`Offline/CRVDQM/src/SConscript` is the scons template. scons does not pass
libraries on transitively, so name everything the `.cc` files use, starting
with `mu2e_DQMHelpers`:

```python
# Offline/CRVDQM/src/SConscript (excerpt)
helper.make_mainlib([
    'mu2e_DQMHelpers',
    'mu2e_CRVConditions',
    'mu2e_CosmicRayShieldGeom',
    'mu2e_RecoDataProducts',
    'mu2e_DataProducts',
    'artdaq-core-mu2e_Overlays',
    'CLHEP',
    # ... art, fhiclcpp, cetlib, rootlibs: copy the rest of the list
])
```

Finally, register the package in Offline's top-level `Offline/CMakeLists.txt`,
in alphabetical order. For the CRV that is:

```cmake
add_subdirectory(CRVConfig)
add_subdirectory(CRVDQM)
add_subdirectory(CRVFilters)
```

scons needs no registration: it finds `Offline/TrkDQM/src/SConscript` itself.

**3. Write the modules.** They live outside Offline: offline ones in
[Mu2e/DQM](https://github.com/Mu2e/DQM) (`src/`), online ones in
[Mu2e/otsdaq-mu2e-dqm](https://github.com/Mu2e/otsdaq-mu2e-dqm)
(`otsdaq-mu2e-dqm/ArtModules/`). Every module drives its client the same way:

```cpp
fhicl::Table<DQMClientFhicl> dqm{fhicl::Name("dqm")};   // its whole FHiCL schema
TrkDigiDQM dqm_{toConfig(conf().dqm().hists())};

beginJob:    dqm_.Book(*tfs);                 // or tfs->mkdir("...")
beginSubRun: dqm_.BeginSubRun(sr.run(), sr.subRun());
analyze:     dqm_.Fill(*event.getValidHandle<StrawDigiCollection>(tag_));
endSubRun:   dqm_.EndSubRun();
endJob:      dqm_.EndJob();                   // series written, diagnostics printed
beginRun:    dqm_.ResetForNewRun();           // online only
publish:     for (auto& [group, copies] : dqm_.hists().publishedCopies()) ...
```

The CRV modules to copy from:
- Mu2e/DQM `src/DqmCrvStatus_module.cc`: the shortest, a complete offline
  module in under 100 lines.
- Mu2e/DQM `src/DqmCrvDigi_module.cc`: recomputes the layout on each new run
  and passes it to the client.
- otsdaq-mu2e-dqm `ArtModules/CrvDQM_module.cc`: online, with publishing,
  HistoSender and the THttpServer display.

In the module's own build files, link `Offline::TrkDQM` (CMake) and
`mu2e_TrkDQM` (scons) as well as `Offline::DQMHelpers` / `mu2e_DQMHelpers`,
since the module includes both packages' headers. The `DqmCrv*` blocks in
Mu2e/DQM `src/CMakeLists.txt` and `src/SConscript` show both.

**4. Add FHiCL tables** in `Offline/TrkDQM/fcl/prolog.fcl`, the new
subdetector package. That file includes the generic presets from
`Offline/DQMHelpers/fcl/prolog.fcl` (this package) and builds one table per
client out of them. For the CRV:

```
# Offline/CRVDQM/fcl/prolog.fcl (excerpt)
#include "Offline/DQMHelpers/fcl/prolog.fcl"

BEGIN_PROLOG

CRVDQM : {
  # offline: one job-integrated copy of everything
  Digi   : { hists : @local::DQMHelpers.Hists.JobOnly }
  Status : { hists : @local::DQMHelpers.Hists.JobOnly }

  # online: the generic online base plus rules choosing what is published
  DigiOnline : {
    hists : {
      @table::DQMHelpers.Hists.Online
      rules : [
        { match : "h1_peakAdc"  modes : [ "job" ]  publish : true },
        { match : "dtFpgaPairs" modes : [ "job" ]  publish : true  group : "timing_fpga" }
      ]
    }
  }
}

END_PROLOG
```

A job then includes the subdetector prolog and hands each module its table.
Adapted from Mu2e/DQM `fcl/crvDQM.fcl`:

```
#include "Offline/CRVDQM/fcl/prolog.fcl"

physics.analyzers.crvStatusDQM : {
  module_type : DqmCrvStatus
  statusTag   : "CrvDigi"
  daqErrorTag : "CrvDigi"
  dqm         : @local::CRVDQM.Status
}
# per-subrun copies instead, overriding below the include:
physics.analyzers.crvStatusDQM.dqm.hists : @local::DQMHelpers.Hists.PerSubrun
```

A subdetector that needs no named tables can skip `Offline/TrkDQM/fcl` and use
the presets in `Offline/DQMHelpers/fcl/prolog.fcl` directly, e.g.
`dqm : { hists : @local::DQMHelpers.Hists.JobOnly }`.

**5. Write the reference catalogue**, `Offline/TrkDQM/data/TrkDQM_binning_v1.txt`:
run a module with `catalogueFile` set and keep the dump, as in
[Changing binning](#changing-binning).

## Core API

### `DQMAxis`

```cpp
constexpr DQMAxis(int nBins, double lo, double hi);
static constexpr DQMAxis Counts(int first, int last);     // one bin per integer
static constexpr DQMAxis Symmetric(double range, double binWidth);  // e.g. a dt axis
std::string describe() const;                             // "n,lo,hi"
```

### `DQMHistSet`

Booking (only in `book()`):

```cpp
DQMH1<H> book1<H>(path, title, const DQMAxis& x);
DQMH2<H> book2<H>(path, title, const DQMAxis& x, const DQMAxis& y);
// A summary filled once at end of job (a fit result, a rate): no rule applies,
// because segment copies of it would all be identical. Still in the catalogue.
DQMH1<H> bookSummary1<H>(path, title, x);
DQMH2<H> bookSummary2<H>(path, title, x, y);
```

`path` may name a subdirectory (`"timing/dtFpgaPairs"`); rules match the whole
path. The returned handle fans one `Fill` out to every copy the rules asked
for, and converts to `H*` (the job copy), so accessors keep working:

```cpp
DQMH1<TH1F> h = hists().book1<TH1F>(...);
h.Fill(x);  h.Fill(x, weight);  h.ForEach([](TH1F* p) { ... });  // styling, labels
TH1F* job = h;                                                   // the job copy
```

The job copy always exists, so an accessor never returns null; a rule's `job`
mode decides only whether it is *written*.

Reading and lifecycle (the client base class drives most of this):
`copies(path)`, `allCopies()`, `live(path)`, `publishedCopies()`,
`Advance(event, ewt)`, `BeginSubRun`, `EndSubRun`, `RefreshLabels`,
`Finalize`, `ResetContents`, `WriteCatalogue(ostream&)`.

### `DQMClient`

Owns the set, the series and the diagnostics.

| Module calls | Client overrides (protected) |
|---|---|
| `Book(dir)` — books, then freezes | `book()` — required |
| `BeginSubRun(run, subrun)` / `EndSubRun()` | `endSubRun()` |
| `EndJob()` — hook, finalize, write series, print diagnostics | `endJob()` |
| `ResetForNewRun()` — online | `resetForNewRun()` |
| `Fill(...)` — the client's own typed method | `beginEvent(clock)` first in `Fill` |

Also `hists()`, `series()`, `diag()`, `nEvents()`, `run()`, `subrun()`,
`booked()`, `name()`, `binningVersion()`. `Book()` writes a
`dqmBinningVersion` `TNamed` into the client's directory and books `nEvents`
for you — do not book your own.

### `DQMDiagnostics`

With fixed axes, an out-of-range value lands in an overflow bin; this is how
it gets *reported*.

```cpp
diag().Count("clusterOutsideEnvelope", "a cluster position is outside ...");
diag().Note("maxAbsDt", dt);   // largest value seen for a key
```

`Count` warns the first time per key and totals are printed by `EndJob()`.
Use it for conditions the data can produce. A condition that means the job is
misconfigured should throw instead.

### `DQMSeries` (online only)

Graphs do not merge — `hadd` concatenates their points — so they are outside
the fixed-set rule and offline presets leave them off.

```cpp
DQMSeries& g = series().book("graphs/g_rate", "Rate;EWT;Hz");  // path, title, cap
g.Add(x, y);  g.AddIfChanged(x, y);  g.Step(x, y);             // step: status words
```

With `liveSeries` off, `book()` returns an inert series and every `Add` is a
no-op, so a client fills unconditionally. `EndJob()` writes them.

### `DQMStyle`

The ROOT translation of `mu2e.mplstyle`. It lives here, not in the online
package, so a histogram looks the same on the shifter's page and in whatever
offline job later draws the merged file — the fixed binning exists so those two
can be compared, and they are easier to compare when they are drawn alike.

```cpp
DQMStyle::SetStyle();                      // process-wide; once, by whoever draws
DQMStyle::FormatHist(h, "blue");           // black (default), blue, green, red
DQMStyle::FormatHist2D(h2);                // axis offsets, stats off
DQMStyle::FormatGraph(g, "red");
```

Nothing here books, fills or interprets a histogram, and no client calls it:
styling is the display's job, so the module that owns the canvases calls it.
`SetStyle()` ends in `gROOT->ForceStyle()` and so is process-wide — in a job
that monitors several subdetectors, the last caller wins, which is an argument
for this one style rather than a style per package.

## FHiCL

A client's whole schema is its `hists` table (`DQMClientFhicl` when it has no
other knobs). There is no binning in it and no way to drop a histogram.

```
hists : {
  annotateTitles : true       # append the range each copy covers to its title
  stampMetadata  : true       # add the dqmSegment TNamed to every copy
  subrunDir      : "bySubrun"
  segmentDir     : "segments"
  liveSeries     : false      # book the client's graphs (online only)
  rules : [                   # omit for job-only, which is what hadd wants
    { match    : "nDigis"     # glob over the directory-qualified path
      modes    : [ "job", "window" ]     # any of job, subrun, window
      jobPersist : true                  # write the job copy
      subrun   : { keep : -1  persist : true  persistLive : false }
      window   : { span : 50000  unit : "ewt"   # event | ewt | subrun
                   subdivisions : 10  keep : 4
                   persist : false  persistLive : false }
      liveName : "nDigisLastEwt"         # name of the rolling copy
      publish  : true                    # hand its copies to the consumer
      group    : "occupancy"             # collect them under one label
      archiveGroup : "occupancy_history" # label for _prevN / per-subrun copies
    }
  ]
}
```

First matching rule wins. `publish`/`group` are what `publishedCopies()`
returns; the registry does not know what a consumer does with them, which is
what lets DQMHelpers build in the DAQ process without knowing HistoSender
exists.

Presets, from `fcl/prolog.fcl`:

| `@local::DQMHelpers.Hists.` | What it gives |
|---|---|
| `JobOnly` | one job-integrated copy of everything (offline default) |
| `PerSubrun` | plus a copy per subrun, named by run and subrun |
| `Online` | annotated titles, no stamp, `liveSeries`; add `rules` for publishing |

```
#include "Offline/DQMHelpers/fcl/prolog.fcl"
physics.analyzers.myDQM.dqm.hists : @local::DQMHelpers.Hists.PerSubrun
# or splice and extend:
physics.analyzers.myDQM.dqm.hists : { @table::DQMHelpers.Hists.Online  rules : [ ... ] }
```

### What the copies are called

| Copy | Name | Written when |
|---|---|---|
| job | `<name>` | `jobPersist` (default true) |
| subrun, in progress | `<name>_sub` | `subrun.persistLive` |
| subrun, finished | `bySubrun/<name>_rNNNNNN_sNNNNNN` | `subrun.persist` |
| window, rolling | `liveName`, else `<name>_last` | `window.persistLive` |
| window, previous spans | `segments/<name>_prev1 … _prevK` | `window.persist` |

`_prev1` is always the most recent completed span, so a GUI can subscribe to a
fixed name. Only the job copy and the per-subrun archives are `hadd`-safe: the
rest overlap in time or are rebuilt in place. With `annotateTitles`, every copy
says what it covers, e.g. `[last 50000 EWT: 1234-5678, 4321 events]`, and with
`stampMetadata` the same is machine-readable in a `dqmSegment` `TNamed`.

## Changing binning

1. Edit the `DQMAxis` constant (or add/remove a histogram) in the client header.
2. Bump that client's `kBinningVersion`.
3. Regenerate the catalogue: run a job with the module's `catalogueFile` set
   and copy the result over the client package's reference file (for the CRV,
   `CRVDQM/data/CRVDQM_binning_v<N>.txt`).
4. Say so in the PR: merging tools and the DQM database see a new version, and
   old files no longer combine with new ones.

`WriteCatalogue` prints one line per histogram — path, class, axes — so
`diff` against the reference file is the check that a job's set is the
expected one.
