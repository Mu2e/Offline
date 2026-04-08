# CalPatRec Code Review

Detailed review of the CalPatRec package covering all modules, algorithms, data
structures, and utilities. Issues are organized by severity.

---

## Critical Issues

### 1. Buffer overflow in `DeltaCandidate::removeSeed()` — unbounded array walk

**File:** `src/DeltaCandidate.cc:49-55`

```cpp
if (Station == fFirstStation) {
  while (fSeed[++fFirstStation] == nullptr) {}
}
if (Station == fLastStation) {
  while (fSeed[--fLastStation] == nullptr) {}
}
```

`fSeed` has size `kNStations` (18). If all remaining seeds are null after the
removed one, the `while` loop increments/decrements past the array bounds
without any bound check. `++fFirstStation` can exceed 17;
`--fLastStation` can go below 0. This is undefined behavior and can cause a
crash or memory corruption.

---

### 2. Stale pointer after loop in `CalHelixFinderAlg` — wrong face used

**File:** `src/CalHelixFinderAlg.cc:2485`

After the loop at lines 2378-2479 (`for (int f=SeedIndex.face; f>=0; --f)`),
`facez` points to the last face the loop touched — typically face 0 — **not**
`ibest.face`:

```cpp
if (ibest.panel >= 0){
  panelz = &facez->panelZs[ibest.panel];   // BUG: facez ≠ ibest.face
```

Fix: insert `facez = &Helix._oTracker[ibest.face];` before using `facez`.

---

### 3. Division-by-zero in `DeltaSeed` accessor methods

**File:** `inc/DeltaSeed.hh:91-93, 106, 110-111`

```cpp
float Chi2ParN ()  { return fChi2Par /fNHits; }
float Chi2PerpN()  { return fChi2Perp/fNHits; }
float Chi2TotN ()  { return (fChi2Par+fChi2Perp)/fNHits; }
float EDep    ()   { return fSumEDep/fNStrawHits ; }
float TMean   ()   { return fSumT/fNHits;  }
float Chi2Time()   { return (fSumT2/fNHits-(fSumT/fNHits)*(fSumT/fNHits))/fSigT2; }
```

None of these check for zero denominators. `fNHits`, `fNStrawHits`, and
`fSigT2` can all be zero, producing division-by-zero (FPE or inf/NaN that
propagates).

The same pattern appears in:
- `inc/DeltaCandidate.hh:64` — `EDep()` divides by `fNStrawHits`
- `inc/ProtonCandidate.hh:83-87` — `eDep()`, `xMean()`, `yMean()`, `FBest()`
  all divide without checks

---

### 4. Division by zero — zero determinant in `DeltaSeed::CalculateCogAndChi2()`

**File:** `src/DeltaSeed.cc:165-168`

```cpp
double d  = fSnx2*fSny2-fSnxy*fSnxy;
double xc = (fSnyr*fSnx2-fSnxr*fSnxy)/d;
double yc = (fSnyr*fSnxy-fSnxr*fSny2)/d;
```

If the determinant `d` is zero (collinear geometry), division by zero occurs.
Same pattern in `src/DeltaCandidate.cc` at lines 67-69, 137-139, and 215-217.

---

### 5. Division by zero — `rho` calculation in `DeltaCandidate`

**File:** `src/DeltaCandidate.cc:74-76` (also 144-146 and 222-224)

```cpp
double rho = sqrt(xc*xc+yc*yc);
fNx        = xc/rho;
fNy        = yc/rho;
```

If the center-of-mass is at the origin, `rho == 0` and division by zero occurs.

---

### 6. Null pointer dereference in `PrefetchData_module`

**File:** `src/PrefetchData_module.cc:233-237`

```cpp
if (_shcol) {                                          // only checks _shcol
  ...
  const StrawHitPosition& shp = _shpcol->at(ish);     // dereferences _shpcol!
```

`_shpcol` is only set when `_fetchStrawHitPositions` is true (line 192-194),
but it is dereferenced unconditionally inside the `if (_shcol)` block. If
`_fetchStrawHits` is true but `_fetchStrawHitPositions` is false, this is a
null pointer dereference.

---

### 7. Division by zero in `PhiClusterFinder::initCluster()`

**File:** `src/PhiClusterFinder_module.cc:447-451`

```cpp
tacc/=weight;
tacc2/=weight;
xacc/=weight;
yacc/=weight;
zacc/=weight;
```

No check that `weight > 0` before division. If `_strawHitIdxs` is empty, weight
is 0. Note that `CalLineTimePeakFinder` (line 326) properly checks
`if(weight > 0.)` before dividing.

---

## High Severity Issues

### 8. Dead code — `isHitUsed()` always returns 0

**File:** `src/CalHelixFinderAlg.cc:124-134`

```cpp
int CalHelixFinderAlg::isHitUsed(int index) {
  return 0;          // Everything below is unreachable
  // ...commented-out logic...
}
```

The function unconditionally returns 0. The rest is dead code.

---

### 9. Memory leak — `DeltaFinderAlg` allocated with `new`, never deleted

**File:** `src/DeltaFinder_module.cc:170`

```cpp
_finder = new DeltaFinderAlg(config().finderParameters,&_data);
```

`_finder` is a raw pointer allocated with `new` but the module has no
destructor or `delete` call. Should use `std::unique_ptr`.

---

### 10. Out-of-bounds array access after detecting invalid indices

**File:** `src/DeltaFinderAlg.cc:732-743`

```cpp
if (_printErrors) {
  if ((os < 0) || (os >= kNStations)) printf("ERROR: ...\n");
  if ((of < 0) || (of >= kNFaces   )) printf("ERROR: ...\n");
  if ((op < 0) || (op >= kNPanelsPerFace)) printf("ERROR: ...\n");
}
// BUG: continues to use invalid indices unconditionally:
FaceZ_t* fz = &_data->fFaceData[os][of];
```

Invalid indices are detected and printed but execution continues with the
out-of-bounds access. The code also has a `FIXME` comment acknowledging this.

---

### 11. Incomplete state update in `DeltaSeed::ReplaceFirstHit()`

**File:** `src/DeltaSeed.cc:131-145`

`ReplaceFirstHit` updates `fSumEDep`, `fSumT`, `fSumT2` but does **not**
update the coordinate-sum accumulators `fSnx2`, `fSnxy`, `fSny2`, `fSnxr`,
`fSnyr`. Those stale values are subsequently used in `CalculateCogAndChi2()`
(line 165), producing incorrect geometry.

---

### 12. Division by zero — time-vs-z linear fit

**File:** `src/DeltaCandidate.cc:97`

```cpp
fDtDz = (tzm-tm*zm)/(z2m-zm*zm);
```

Denominator `z2m - zm*zm` is zero when all hits are at the same z-coordinate
(single station). Same issue at lines 167 (`AddSeed`), 242
(`MergeDeltaCandidate`), and in `src/ProtonCandidate.cc:340`.

---

### 13. Unchecked weight used as divisor in `CalHelixFinderAlg`

**File:** `src/CalHelixFinderAlg.cc:676, 802`

```cpp
weight = calculatePhiWeight(*hit, helCenter, radius, 0, PhiZInfo.banner);
err    = 1./sqrt(weight);
```

If `weight` is zero, `sqrt(0)=0` and `1./0.` is infinity.

---

### 14. Division by zero — `_dfdz` in helix calculations

**File:** `src/CalHelixFinderAlg.cc:1519, 2688, 2695`

```cpp
float lambda(1./Helix._dfdz);             // line 1519
float tollMax = fabs(2.*M_PI/dfdz);       // line 2688
```

`_dfdz` can be zero, producing infinity/NaN. No validation before division.

---

### 15. Division by zero in `defineHelixParams`

**File:** `src/CalHelixFinderAlg.cc:86, 98`

```cpp
pvec[HelixTraj::omegaIndex]  = amsign/radius;
pvec[HelixTraj::tanDipIndex] = amsign/(radius*Helix._dfdz);
```

No validation of `radius` or `_dfdz` before these divisions.

---

### 16. Division by zero in `TZClusterFinder`

**File:** `src/TZClusterFinder_module.cc:353, 816`

```cpp
comboHit.hWeight = 1/(hit->timeVar());
```

No validation that `timeVar()` is non-zero.

---

## Medium Severity Issues

### 17. Memory leak — ROOT objects in `TZClusterFinder`

**File:** `src/TZClusterFinder_module.cc:212, 385, 404, 422`

```cpp
_c1 = new TCanvas("_c1", "t vs. z", 900, 900);
```

The destructor at line 224 is empty. `TCanvas`, `TGraph`, and `TF1` objects
allocated with `new` in `plotTZ()` are never deleted.

---

### 18. Signed/unsigned cast issues in loop variables

**File:** `src/TZClusterFinder_module.cc:545, 555, 667`

```cpp
for (int i=(int)_f.cHits.size()-1; i>=0; i--)
size_t start = (size_t)_f.startIndex;          // int(-1) → SIZE_MAX
```

Casting `size_t` to `int` can overflow; casting `-1` to `size_t` creates
`SIZE_MAX`.

---

### 19. Missing bounds check in `ComboHitFilter`

**File:** `src/ComboHitFilter_module.cc:157`

```cpp
int ind = ch->indexArray().at(0);
```

Assumes `indexArray()` is non-empty. No size check before `.at(0)`.

---

### 20. Assertions disabled in release builds

**File:** `src/DeltaSeed.cc:133, 163`

```cpp
assert(fSFace[0] == Hd->fZFace);
assert(fNHits > 2);
```

`assert` is compiled out in optimized builds. If these conditions fail in
production, silent data corruption follows instead of a clean abort.

---

### 21. `printf`-based error handling in `CalHelixFinder_module`

**File:** `src/CalHelixFinder_module.cc:174-208`

```cpp
} else {
  _chcol = 0;
  printf(" >>> ERROR in CalHelixFinder::findData: ...\n");
}
```

Uses `printf` instead of `mf::LogError`. Sets `_chcol = 0` (should be
`nullptr`) and continues, risking null pointer dereference downstream.

---

### 22. Vector erase during iteration in `TZClusterFinder::combineChunks()`

**File:** `src/TZClusterFinder_module.cc:614-651`

```cpp
for (size_t i=0; i<_f.chunks.size()-1; i++) {
  ...
  _f.chunks.erase(_f.chunks.begin()+chunkTwoIdx);
```

Erasing elements from a vector while iterating with indices can skip elements
or access wrong indices if `chunkTwoIdx <= i`.

---

## Low Severity / Code Quality

### 23. Hardcoded magic numbers

- `src/CalHelixFinderAlg.cc:852` — `float weight_cl = 784.;`
- `src/DeltaFinderAlg.cc:589` — `float max_d2(20*20); // FIXME`
- `src/DeltaFinder_module.cc:334` — `if (pc->nHitsTot() >= 15)`
- `src/DeltaFinderAlg.cc:47` — `if (fabs(tpc-t) > 30)` (undocumented units)
- `src/CalTimePeakFinder_module.cc:263` — `TrkT0(cl_time, 0.1)` (magic error)

### 24. Extensive dead / commented-out code

- `src/CalHelixFinderData.cc:22-67` — commented-out copy constructor
- `src/CalHelixFinderAlg.cc:124-134` — dead code after early return
- `src/CalTimePeakFinder_module.cc:93-95` — commented-out member variables
- `src/PrefetchData_module.cc:144-145, 196-199, 252-258, 260-271` — multiple
  dead blocks

### 25. `NULL` vs `nullptr` inconsistency

- `src/CalHelixFinderData.cc:16` — `_helix = NULL;`
- `inc/DeltaCandidate.hh:61` — `return (fSeed[I] != NULL);`

Modern C++ should use `nullptr` everywhere.

### 26. Implicit `size_t` to `int` conversions

- `src/PhiClusterFinder_module.cc:217, 231, 431` — `int nh = ordchcol.size();`

Risk of overflow on very large collections.

---

## Summary

| Severity | Count | Key categories |
|----------|-------|----------------|
| Critical | 7 | Buffer overflow, stale pointer, division by zero, null deref |
| High | 9 | Memory leak, dead code, missing bounds checks, unvalidated divisions |
| Medium | 6 | ROOT memory leaks, signed/unsigned casts, assertions, vector erase |
| Low | 4 | Magic numbers, dead code, style inconsistencies |
| **Total** | **26** | |

The most impactful issues to fix first are **#1** (buffer overflow in
`removeSeed`), **#2** (stale `facez` pointer), **#6** (null deref in
PrefetchData), and **#3/#4/#5** (the family of division-by-zero issues in
DeltaSeed/DeltaCandidate).
