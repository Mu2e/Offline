#pragma once

// Validation plots for the VDResampler generate modules (VDResamplerGenerateFromModel
// and VDResamplerGenerateMix): for every model dimension the GENERATED distribution is
// compared against the MOTHER (source MC truth) distribution, and four sparse-bin-robust
// agreement metrics (W1 / JSD / TV / KS) plus chi2/NDF are reported.
//
// Two histogram families are compared:
//   * the basis' TRANSFORMED coordinates -- what the model actually learned
//       V1        : xTrans, yTrans, tTrans,      prTrans,   pphiTrans,      pzTrans
//       V2 plain  : xTrans, yTrans, tTrans,      pTotTrans, urTrans,        uphiTrans
//       V2 asinh  : xTrans, yTrans, tTrans,      pTotTrans, urTransAsinh,   uphiTransAsinh
//       V3        : xTrans, yTrans, tTransAsinh, pTotTrans, urTransAsinh,   uphiTransAsinh
//   * the recovered PHYSICAL quantities px, py, pz, x, y, t -- what downstream sees.
// There is no z distribution: every generated sample lands on the VDz0 plane by
// construction, so z is a constant and carries nothing to compare.
// Histograms are BOOKED at a fine native resolution so nothing is thrown away before the
// statistics are known.
//
// Adaptive rebinning (the reason binning is decided late): a generation job's sample
// count is not known until it ends, and the mother file's count is unrelated to it --
// either side may run from a few thousand to O(1e7) entries. Binning fine enough for 1e7
// leaves a few-thousand-entry comparison as noise. So finalize() -- where BOTH counts are
// known -- coarsens every histogram by an integer TH1::Rebin factor chosen from the
// smaller of the two, and only then computes the metrics. Mother and generated always
// take the SAME factor, so their bin edges stay aligned (required for the metrics and for
// any ratio drawn from these hists later).
//
// Only ONE generated distribution is compared per dimension, because a generation job
// commits to a single network via useEMANetworkIfAvailable.
//
// The mother distribution is filled by re-reading the source ROOT file (the same file and
// accepted-hit selection the PtotResampler uses, via forEachAcceptedHitRoot) and pushing
// each hit through forwardTransformSample for the model's basis. The generated
// distribution is filled one sample at a time from the module's produce().
//
// Every dimension -- transformed AND physical -- additionally gets a ready-made comparison
// CANVAS ("c_<suffix>") written beside the histograms: the two distributions overlaid on top
// (log-y where the dynamic range needs it), and the generated/mother ratio underneath with
// the agreement metrics printed in its legend. The canvases are built from the same rebinned,
// normalized histograms the metrics are computed from, so a number and the shape behind it
// are always read together without re-deriving anything downstream.
//
// Alongside the per-dimension 1D comparison, a set of 2D CORRELATION views is booked for
// both sides, so the generated and the training pictures can be read side by side:
//   <pz>_vs_r, pt_vs_r, pr_vs_r, pphi_vs_r, pt_vs_pz, y_vs_x, <pz>_vs_t
// where <pz> is pTot for the V2/V3 pTotal+slopes bases -- there the momentum structure lives
// in pTotal rather than pz -- and pz for V1.
// For photons three further time-vs-energy views zoom on the STM lines of interest
// (1809, 844 and 347 keV, each +/- 15 keV).
// The 2D views carry NO metrics -- W1/JSD/TV/KS and chi2/NDF are defined here for the 1D
// marginals only -- and are NOT rebinned. They are made comparable by eye in two steps:
// each is normalized to unit volume over its side's WHOLE sample (flow bins included, so a
// zoomed window is not rescaled by the fraction of the sample that happens to fall in it),
// and both sides of a pair are then given the same stored z-range, so equal colour means
// equal density. Palette is left to the drawing script (gStyle is a draw-time global that
// does not persist in the file); kInferno is the intended one.
//
// Header-only; depends on ROOT (TH1D/TH2D/TTree/TDirectory) and the VDResampler transforms.
// Yongyi Wu, Aug. 2026

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "art_root_io/TFileDirectory.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPad.h"
#include "TTree.h"

#include "Offline/STMMC/inc/VDResamplerGenerateCommon.hh"
#include "Offline/STMMC/inc/VDResamplerPtotResampler.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"

namespace mu2e {
namespace VDResampler {

// ---------------------------------------------------------------------------
// DistributionMetrics — shape-agreement metrics between a test and a reference
//   histogram, computed DIRECTLY from the (plotted, post-rebin) bin contents so they
//   reflect exactly what is drawn. These four are far more robust to sparse /
//   low-population bins than chi2/NDF, whose (obs-exp)^2/sigma^2 term blows up where
//   counts are few, so they carry the verdict and chi2/NDF is reported only alongside.
//     W1  : Wasserstein-1 (earth-mover) distance, in x-axis UNITS -- its scale is
//           per-variable (large for wide axes like time, small for slopes).
//     JSD : Jensen-Shannon distance in [0,1], symmetric, sqrt of the JS divergence.
//     TV  : total variation 0.5*sum|p-q| in [0,1] -- fraction of prob. mass misplaced.
//     KS  : Kolmogorov-Smirnov DISTANCE in [0,1] = max|CDF_p - CDF_q| (not a p-value).
//   All four are 0 for a perfect match; chi2/NDF ~ 1 for one.
// ---------------------------------------------------------------------------
struct DistributionMetrics {
    double w1   = std::nan("");
    double jsd  = std::nan("");
    double tv   = std::nan("");
    double ks   = std::nan("");
    double chi2 = std::nan("");
    int    ndf  = 0;
    double chi2ndf() const { return ndf > 0 ? chi2 / ndf : std::nan(""); }
};

// Compute the metrics of `hist` (test/generated) against `ref` (mother/truth). Both are
// re-normalized to unit probability mass here, so any prior scaling of either is
// irrelevant -- only the SHAPES are compared. chi2 uses the stored bin errors, falling
// back to the reference content as the Poisson variance when both errors vanish.
inline DistributionMetrics computeDistributionMetrics(const TH1D& hist, const TH1D& ref,
                                                      const std::string& moduleName)
{
    DistributionMetrics m;
    const int nb = std::min(hist.GetNbinsX(), ref.GetNbinsX());
    if (hist.GetNbinsX() != ref.GetNbinsX())
        mf::LogWarning(moduleName)
            << "computeDistributionMetrics: bin-count mismatch (" << hist.GetNbinsX()
            << " vs " << ref.GetNbinsX() << "); comparing the first " << nb << " bins.";
    if (nb <= 0) return m;

    std::vector<double> p(nb), q(nb);
    double sp = 0.0, sq = 0.0;
    for (int b = 1; b <= nb; ++b) {
        p[b - 1] = std::max(hist.GetBinContent(b), 0.0);   // test
        q[b - 1] = std::max(ref.GetBinContent(b),  0.0);   // reference
        sp += p[b - 1];
        sq += q[b - 1];
    }
    if (sp <= 0.0 || sq <= 0.0) return m;   // one side is empty; metrics stay NaN
    for (int i = 0; i < nb; ++i) { p[i] /= sp; q[i] /= sq; }

    // total variation
    double tv = 0.0;
    for (int i = 0; i < nb; ++i) tv += std::abs(p[i] - q[i]);
    m.tv = 0.5 * tv;

    // Jensen-Shannon distance: base-2 divergence in [0,1], distance = sqrt(divergence)
    double div = 0.0;
    for (int i = 0; i < nb; ++i) {
        const double mid = 0.5 * (p[i] + q[i]);
        if (mid <= 0.0) continue;
        if (p[i] > 0.0) div += 0.5 * p[i] * std::log2(p[i] / mid);
        if (q[i] > 0.0) div += 0.5 * q[i] * std::log2(q[i] / mid);
    }
    m.jsd = std::sqrt(std::max(0.0, div));

    // single CDF sweep for both Wasserstein-1 and the KS statistic
    const TAxis* axis = ref.GetXaxis();
    double cp = 0.0, cq = 0.0, w1 = 0.0, ks = 0.0;
    for (int i = 0; i < nb; ++i) {
        cp += p[i];
        cq += q[i];
        const double diff = std::abs(cp - cq);
        ks = std::max(ks, diff);                        // KS = sup |CDF_p - CDF_q|
        w1 += diff * axis->GetBinWidth(i + 1);          // W1 = integral |CDF_p - CDF_q| dx
    }
    m.w1 = w1;
    m.ks = ks;

    // chi2/NDF on the contents as stored
    double chi2 = 0.0;
    int ndf = 0;
    for (int b = 1; b <= nb; ++b) {
        const double obs = hist.GetBinContent(b);
        const double exp = ref.GetBinContent(b);
        double err2 = hist.GetBinError(b) * hist.GetBinError(b)
                    + ref.GetBinError(b)  * ref.GetBinError(b);
        if (err2 == 0.0) {
            if (exp <= 0.0) continue;   // empty in both; contributes nothing
            err2 = exp;                 // fall back to the Poisson variance
        }
        chi2 += (obs - exp) * (obs - exp) / err2;
        ++ndf;
    }
    m.chi2 = chi2;
    m.ndf  = ndf;
    return m;
}

// ---------------------------------------------------------------------------
// Adaptive rebinning
//   targetBinCount: how many bins a comparison of `n` entries can support. ~2*sqrt(n)
//     keeps the mean occupancy climbing with statistics (~sqrt(n)/2 entries per bin)
//     instead of pinning it, so a 1e7-entry comparison keeps fine structure that a fixed
//     target would erase, while a few-thousand-entry one collapses to a readable ~100
//     bins. Clamped to [kMinTargetBins, kMaxTargetBins].
//   rebinFactor: the TH1::Rebin group size bringing `nbins` closest to `target`. It must
//     be an exact DIVISOR of nbins (ROOT makes a non-divisor a fatal error), so the
//     achieved bin count only approximates the target -- how closely depends on how many
//     divisors the native binning happens to have.
// ---------------------------------------------------------------------------
constexpr int kMinTargetBins = 50;
constexpr int kMaxTargetBins = 4000;

inline int targetBinCount(long n) {
    if (n <= 0) return kMinTargetBins;
    const double t = 2.0 * std::sqrt(static_cast<double>(n));
    return std::max(kMinTargetBins, std::min(kMaxTargetBins, static_cast<int>(std::lround(t))));
}

// TH1::Rebin merges `ngroup` adjacent bins and REQUIRES ngroup to divide nbins exactly. 
// Otherwise it raises a fatal error. The factor cannot simply be nbins/target
// rounded: that ratio is almost never a divisor. Search the actual divisors of nbins 
// instead and take the one whose resulting bin count lands closest to the target.
//
// Closeness is judged on the RATIO rather than the difference, because bin counts are
// compared multiplicatively here.
inline int rebinFactor(int nbins, int target) {
    if (target <= 0 || nbins <= 0) return 1;
    if (target >= nbins) return 1;   // already at or below the target; nothing to merge

    int bestFactor = 1;
    double bestScore = -1.0;
    // Divisors come in pairs (f, nbins/f), so scanning to sqrt(nbins) covers all of them.
    for (int f = 1; f * f <= nbins; ++f) {
        if (nbins % f != 0) continue;
        for (int candidate : {f, nbins / f}) {
            const int resulting = nbins / candidate;
            const double ratio = (resulting > target)
                ? static_cast<double>(resulting) / target
                : static_cast<double>(target) / resulting;
            if (bestScore < 0.0 || ratio < bestScore) {
                bestScore = ratio;
                bestFactor = candidate;
            }
        }
    }
    return bestFactor;
}

// ---------------------------------------------------------------------------
// Per-dimension histogram specification: name suffix, x-axis title and native binning.
// `transformed` marks the model's own coordinates (as opposed to the recovered physical
// quantities), which decides how a sample is routed when filling. `binWidthUnit` is the
// unit shown in the "counts per <width><unit>" y-axis label, which finalize() writes
// once the rebin factor is known.
// ---------------------------------------------------------------------------
struct DimSpec {
    std::string suffix;         // e.g. "pTotTrans", "px"
    std::string xTitle;         // x-axis title, e.g. "px [MeV]"
    std::string binWidthUnit;   // "" (dimensionless), " MeV", "mm", "ns"
    int    nbins = 0;
    double lo = 0.0, hi = 0.0;
    bool   transformed = false;
    int    slot = 0;            // transformed: index into (xTrans,yTrans,tTrans,m0,m1,m2)
                                // physical:    index into (px,py,pz,x,y,t)
};

// Physical-quantity slots, in the order fillGenerated takes them.
enum PhysicalSlot { kPx = 0, kPy = 1, kPz = 2, kX = 3, kY = 4, kT = 5 };

// TransformedDimStats / TransformedStatsBySlot are defined in VDResamplerGenerateCommon.hh
// (included above), where collectTransformedStats builds them from the loaded models. They
// live there rather than here because this header includes that one, not the other way
// round.

// ---------------------------------------------------------------------------
// Axis sizing from the training statistics.
//
// The static ranges below were hand-picked to be safely wide for every sample, which makes
// them far too wide for most: a distribution occupying 0.2 units of a +/-5 axis lands in a
// handful of bins and the comparison sees nothing. When the model reports its training
// statistics, the axis is instead sized to where the population actually is, BEFORE any
// sample has been generated.
//
// The range is the union of two requirements:
//   * the observed data extent [min, max] — nothing the training set contained should fall
//     outside the axis, since the model can reproduce it;
//   * a mean +/- kSigmaSpan sigma band — so a heavy-tailed dimension whose min/max are set
//     by a couple of outliers still devotes most of the axis to the bulk.
// Taking the union (widest of the two) means neither a compact-but-outliered distribution
// nor a broad one gets clipped. A kRangePadding margin is added so the extremes are not
// exactly on the boundary, and the result is snapped outward to a round number so the axis
// reads sensibly and the bin edges are not arbitrary.
//
// Generated samples CAN fall outside the training extent — the model is not bounded by its
// training data — so the padded range is deliberately wider than [min,max], and anything
// beyond it still lands in the flow bins and is counted by the metrics' normalization.
// ---------------------------------------------------------------------------
constexpr double kSigmaSpan     = 6.0;   // half-width of the bulk band, in sigma
constexpr double kRangePadding  = 0.10;  // fractional margin added to the chosen half-width
constexpr int    kStatsBinCount = 2400;  // native bins for a stats-sized axis (highly composite)

// Snap `x` outward to a multiple of `step` (down for the low edge, up for the high edge),
// so axis limits read as round numbers rather than -4.8317. Only ever widens.
inline double snapOutward(double x, double step, bool up) {
    if (!(step > 0.0) || !std::isfinite(x)) return x;
    return step * (up ? std::ceil(x / step) : std::floor(x / step));
}

// A 1/2/5 x 10^k step that divides `width` into roughly kSnapDivisions parts — the
// granularity the axis edges are snapped to. Tied to the width so a 0.2-wide axis snaps to
// 0.01 and a 20-wide one to 1, instead of one absolute step suiting neither.
constexpr int kSnapDivisions = 20;
inline double niceStep(double width) {
    if (!(width > 0.0)) return 0.0;
    const double raw = width / kSnapDivisions;
    const double decade = std::pow(10.0, std::floor(std::log10(raw)));
    const double scaled = raw / decade;        // in [1,10)
    for (const double s : {1.0, 2.0, 5.0})
        if (scaled <= s) return s * decade;
    return 10.0 * decade;
}

// Apply `stats` to a spec's range, leaving nbins/lo/hi untouched when they are unusable.
// The spec's static range is the fallback, so a model without statistics behaves exactly
// as before.
//
// The two edges are computed and snapped INDEPENDENTLY rather than as a half-width about a
// centre: these distributions are frequently asymmetric (log-time sits well off zero, pTot
// is skewed), and forcing a symmetric range about the midpoint would pad the short side by
// as much as the long one — for a -4..6 dimension that inflates a 10-wide axis to 20, worse
// than the static fallback it replaces.
inline void applyStatsToSpec(DimSpec& spec, const TransformedDimStats& stats) {
    if (!stats.valid) return;
    if (!std::isfinite(stats.min) || !std::isfinite(stats.max)) return;
    if (!std::isfinite(stats.mean) || !std::isfinite(stats.stdev)) return;
    if (!(stats.max > stats.min)) return;   // degenerate/constant dimension

    // The observed extent is the authority: the model was trained on exactly this data, so
    // the axis must cover it, and there is no reason to reach beyond it. The sigma band only
    // TIGHTENS the axis, never widens it — for a dimension whose extent is set by a handful
    // of far outliers, clipping to the bulk band keeps the bins where the population is (the
    // outliers then land in the flow bins, still counted by the normalization). Clamping the
    // band to [min,max] is what keeps a large-sigma dimension from producing an axis wider
    // than the static fallback it replaces.
    double lo = stats.min;
    double hi = stats.max;
    if (stats.stdev > 0.0) {
        lo = std::max(lo, stats.mean - kSigmaSpan * stats.stdev);
        hi = std::min(hi, stats.mean + kSigmaSpan * stats.stdev);
        // A degenerate band (sigma tiny next to the extent, or a mean far outside it) would
        // invert or collapse the range; fall back to the full extent in that case.
        if (!(hi > lo)) { lo = stats.min; hi = stats.max; }
    }

    const double width = hi - lo;
    if (!(width > 0.0)) return;

    // Pad each side, then snap outward so the edges are round numbers.
    const double pad  = kRangePadding * width;
    const double step = niceStep(width);
    spec.lo    = snapOutward(lo - pad, step, /*up=*/false);
    spec.hi    = snapOutward(hi + pad, step, /*up=*/true);
    spec.nbins = kStatsBinCount;
}

// The six transformed dimensions for `basis`. Slot order is
// (xTrans, yTrans, tTrans, m0, m1, m2) — the TRANSFORM's order, not the model's internal
// (t,x,y,...) data layout.
//
// The ranges here are STATIC FALLBACKS, wide enough for any sample. When `stats` carries a
// slot's training statistics the axis is re-sized to the actual population (see
// applyStatsToSpec); slots without statistics keep the fallback.
inline std::vector<DimSpec> transformedDimSpecs(MomentumBasis basis,
                                                const TransformedStatsBySlot* stats = nullptr) {
    const bool asinhSlopes = basisUsesAsinhSlopes(basis);
    const bool asinhTime   = basisUsesAsinhTime(basis);
    const bool isV1        = (basis == MomentumBasis::V1_CylindricalTransformed);

    std::vector<DimSpec> specs;
    specs.push_back({"xTrans", "xTrans", "", 2500, -5., 5., true, 0});
    specs.push_back({"yTrans", "yTrans", "", 2500, -5., 5., true, 1});

    // V3 asinh-tames the log-time, so its time slot holds a DIFFERENT quantity than the
    // plain log-time the other bases store; it gets a distinct name accordingly.
    if (asinhTime)
        specs.push_back({"tTransAsinh", "tTransAsinh", "", 2000, -6., 6., true, 2});
    else
        specs.push_back({"tTrans", "tTrans", "", 2400, -4., 8., true, 2});

    if (isV1) {
        specs.push_back({"prTrans",   "prTrans",   "", 3000, -2., 4., true, 3});
        specs.push_back({"pphiTrans", "pphiTrans", "", 3000, -3., 3., true, 4});
        specs.push_back({"pzTrans",   "pzTrans",   "", 2000, -6., 6., true, 5});
    } else {
        // Slopes are ~0 for a forward beam; the asinh variants compress the wide-angle tails.
        specs.push_back({"pTotTrans", "pTotTrans", "", 2400, -6., 6., true, 3});
        specs.push_back({asinhSlopes ? "urTransAsinh"   : "urTrans",
                         asinhSlopes ? "urTransAsinh"   : "urTrans",   "", 2000, -6., 6., true, 4});
        specs.push_back({asinhSlopes ? "uphiTransAsinh" : "uphiTrans",
                         asinhSlopes ? "uphiTransAsinh" : "uphiTrans", "", 2000, -6., 6., true, 5});
    }

    if (stats)
        for (DimSpec& spec : specs)
            applyStatsToSpec(spec, (*stats)[spec.slot]);
    return specs;
}

// The six recovered physical quantities. Two windows are widened for neutrons: the momentum
// one, whose spectrum extends well past the few-MeV range that covers the other species, and
// the time one, since neutrons are slow enough to keep arriving long after the prompt species
// have. Both keep the other species' bin width, so only the range changes, not the resolution.
inline std::vector<DimSpec> physicalDimSpecs(int pdgId) {
    const bool isNeutron = (pdgId == 2112);
    const double pmax  = isNeutron ? 50. : 8.;
    const int    pnbin = 1000;
    const double tmax  = isNeutron ? 20000. : 5000.;
    const int    tnbin = isNeutron ? 2000 : 500;   // 10 ns per bin either way

    std::vector<DimSpec> specs;
    specs.push_back({"px", "px [MeV]", " MeV", 2 * pnbin, -pmax, pmax, false, kPx});
    specs.push_back({"py", "py [MeV]", " MeV", 2 * pnbin, -pmax, pmax, false, kPy});
    specs.push_back({"pz", "pz [MeV]", " MeV", pnbin, 0., pmax, false, kPz});
    specs.push_back({"x",  "x [mm]", "mm", 800, -3904. - 2000., -3904. + 2000., false, kX});
    specs.push_back({"y",  "y [mm]", "mm", 800, -2000., 2000., false, kY});
    specs.push_back({"t",  "t [ns]", "ns", tnbin, 0., tmax, false, kT});
    return specs;
}

// Full per-dimension comparison list: transformed coordinates + physical quantities.
// `stats` (optional) re-sizes the TRANSFORMED axes to the model's training population; the
// physical axes are always the static ranges, since those are in detector units the model
// statistics say nothing about.
inline std::vector<DimSpec> perDimensionSpecs(MomentumBasis basis, int pdgId,
                                              const TransformedStatsBySlot* stats = nullptr) {
    std::vector<DimSpec> specs = transformedDimSpecs(basis, stats);
    const std::vector<DimSpec> phys = physicalDimSpecs(pdgId);
    specs.insert(specs.end(), phys.begin(), phys.end());
    return specs;
}

// ---------------------------------------------------------------------------
// 2D correlation views
//
// The quantities the 2D views plot are DERIVED from the six physical values both sides
// already produce; nothing new has to be read or generated. Derived2D holds them so the
// mother and generated fill paths compute them identically, from one place.
//   r    : in-plane radius about the detector axis, sqrt((x-x0)^2 + y^2), measured about
//          the ACTUAL x0/y0 the transform was built with rather than a nominal centre, so
//          it stays correct if the geometry ever moves.
//   pt   : transverse momentum sqrt(px^2 + py^2).
//   pr / pphi : the local cylindrical projection of (px,py) about (dx,dy) -- the same
//          decomposition the V1 basis uses, via the shared cartesianToLocalPolar.
//   Ek   : kinetic energy in keV, sqrt(|p|^2 + m^2) - m. keV (not MeV) because the STM
//          lines these views zoom on are quoted in keV. The rest mass is not looked up
//          here: it is passed into book() by the module, which already holds it from the
//          GlobalConstantsService ParticleDataList (pdt->particle(pdgId).mass()), so the
//          conditions DB stays the single source of truth for it.
//   pTot : |p|, the z-axis variable for the V2/V3 bases.
// ---------------------------------------------------------------------------
struct Derived2D {
    double r = 0.0, pt = 0.0, pr = 0.0, pphi = 0.0, pTot = 0.0, ek = 0.0;
};

inline Derived2D computeDerived2D(double x, double y, double px, double py, double pz,
                                  double x0, double y0, double massMeV)
{
    Derived2D d;
    const double dx = x - x0, dy = y - y0;
    d.r    = std::sqrt(dx * dx + dy * dy);
    d.pt   = std::sqrt(px * px + py * py);
    cartesianToLocalPolar(px, py, dx, dy, d.r, d.pr, d.pphi);
    d.pTot = std::sqrt(px * px + py * py + pz * pz);
    d.ek   = 1000. * (std::sqrt(d.pTot * d.pTot + massMeV * massMeV) - massMeV);
    return d;
}

// Which derived (or physical) quantity an axis of a 2D view plots.
enum class Var2D { kR, kPt, kPr, kPphi, kPz, kPTot, kX, kY, kT, kEk };

// One 2D view: name, the two variables, and each axis' binning.
struct Dim2DSpec {
    std::string name;        // named <y>_vs_<x>, e.g. "pt_vs_r" for pt on y against r on x
    std::string xTitle, yTitle;
    Var2D  xVar = Var2D::kR, yVar = Var2D::kPt;
    int    xbins = 0;  double xlo = 0.0, xhi = 0.0;
    int    ybins = 0;  double ylo = 0.0, yhi = 0.0;
};

// The 2D views for this basis and species. Unlike the 1D hists these are binned once, at
// book time, and never rebinned: a generation job's sample count is not known that early,
// so the binning is fixed at the finer end of what is useful -- an over-fine 2D view is
// still readable, while an over-coarse one has already thrown the structure away.
inline std::vector<Dim2DSpec> correlationDimSpecs(MomentumBasis basis, int pdgId) {
    const bool isNeutron = (pdgId == 2112);
    const bool isPhoton  = (pdgId == 22);
    const double pmax  = isNeutron ? 50. : 8.;
    const int    pnbin = 1000;
    const double pbin  = pmax / pnbin;             // the 1D momentum bin width
    const int    p2bin = pnbin / 5;                // 2D views use 5x coarser momentum bins

    // "<width> MeV", the per-bin width quoted in the momentum axis labels.
    auto perMeV = [](double w) {
        char buf[64];
        snprintf(buf, sizeof(buf), "%.03f MeV", w);
        return std::string(buf);
    };
    const std::string pPer = perMeV(5 * pbin);

    // The V2/V3 bases carry the momentum structure in pTotal, so the r- and t-correlated
    // views are drawn against pTot there and against pz for V1 -- and the histogram is
    // NAMED for whichever it holds, so a plot is never mislabelled.
    const bool isV2 = (basis != MomentumBasis::V1_CylindricalTransformed);
    const Var2D  zVar  = isV2 ? Var2D::kPTot : Var2D::kPz;
    const std::string zName = isV2 ? "pTot" : "pz";

    std::vector<Dim2DSpec> specs;
    specs.push_back({zName + "_vs_r", "r [mm] per 5mm", zName + " per " + pPer,
                     Var2D::kR, zVar, 400, 0., 2000., p2bin, 0., pmax});
    specs.push_back({"pt_vs_r", "r [mm] per 5mm", "pt per " + pPer,
                     Var2D::kR, Var2D::kPt, 400, 0., 2000., p2bin, 0., pmax});
    // pr and pphi are SIGNED projections, so unlike pt/pTot they need a two-sided y range:
    // a range starting at 0 would silently drop the inward- / backward-going half of the
    // population into the underflow.
    specs.push_back({"pr_vs_r", "r [mm] per 5mm", "pr per " + pPer,
                     Var2D::kR, Var2D::kPr, 400, 0., 2000., 2 * p2bin, -pmax, pmax});
    specs.push_back({"pphi_vs_r", "r [mm] per 5mm", "pphi per " + pPer,
                     Var2D::kR, Var2D::kPphi, 400, 0., 2000., 2 * p2bin, -pmax, pmax});
    specs.push_back({"pt_vs_pz", "pz per " + pPer, "pt per " + pPer,
                     Var2D::kPz, Var2D::kPt, p2bin, 0., pmax, p2bin, 0., pmax});
    specs.push_back({"y_vs_x", "x [mm] per 10mm", "y [mm] per 10mm",
                     Var2D::kX, Var2D::kY, 400, -3904. - 2000., -3904. + 2000.,
                     400, -2000., 2000.});
    // Neutrons arrive late and slow, so their time-momentum view needs a far wider window.
    if (isNeutron)
        specs.push_back({zName + "_vs_t", "t [ns] per 40ns", zName + " per " + perMeV(25 * pbin),
                         Var2D::kT, zVar, 500, 0., 20000., p2bin, 0., 5 * pmax});
    else
        specs.push_back({zName + "_vs_t", "t [ns] per 1ns", zName + " per " + pPer,
                         Var2D::kT, zVar, 500, 0., 500., p2bin, 0., pmax});

    // Photon line zooms: time vs energy in a +/- 15 keV window about each STM line of
    // interest. 1 keV in E gives 30 bins across the window -- enough occupancy per bin for
    // the 2D structure to show at realistic generated statistics; 10 ns in t covers the
    // full 0-5000 ns spill the 1D time plot uses, so a line's time profile is visible in
    // the same view.
    if (isPhoton) {
        const double kLines[] = {1809., 844., 347.};   // keV, in STM order of interest
        const double kHalf = 15.;                      // keV, half-window about each line
        const double kEkBin = 1.0;                     // keV per energy bin
        for (double line : kLines) {
            char nameBuf[64];
            snprintf(nameBuf, sizeof(nameBuf), "t_vs_Ek_%gkeV", line);
            specs.push_back({nameBuf,
                             "E [keV] per 1keV", "t [ns] per 10ns",
                             Var2D::kEk, Var2D::kT,
                             static_cast<int>(std::lround(2 * kHalf / kEkBin)),
                             line - kHalf, line + kHalf,
                             500, 0., 5000.});
        }
    }
    return specs;
}

// Value of `v` for one sample, from its physical values and the derived quantities.
inline double var2DValue(Var2D v, const Derived2D& d,
                         double x, double y, double t, double pz)
{
    switch (v) {
        case Var2D::kR:    return d.r;
        case Var2D::kPt:   return d.pt;
        case Var2D::kPr:   return d.pr;
        case Var2D::kPphi: return d.pphi;
        case Var2D::kPTot: return d.pTot;
        case Var2D::kEk:   return d.ek;
        case Var2D::kPz:   return pz;
        case Var2D::kX:    return x;
        case Var2D::kY:    return y;
        case Var2D::kT:    return t;
    }
    return 0.0;
}

// ---------------------------------------------------------------------------
// ValidationPlots — one comparison set (mother vs generated) for a single
//   (source, pdgId, basis) triple.
//
// Lifetime:
//   book(...)        once, at module construction: creates the hist pair for every
//                    dimension inside a per-set subdirectory of the TFileService file
//                    at native binning, and fills the mother side from the source file.
//   fillGenerated(g, x,y,t,px,py,pz)
//                    once per generated event, from produce().
//   finalize()       at endJob: rebin both sides by the statistics-driven factor,
//                    normalize, compute metrics, write the summary tree, log the table.
// ---------------------------------------------------------------------------
class ValidationPlots {
public:
    // Below this many generated entries the metrics are dominated by sampling noise
    // rather than by model bias, so finalize() flags the comparison as under-sampled.
    static constexpr long kLowStatisticsWarn = 1000;

    bool booked() const { return booked_; }
    const std::string& label() const { return label_; }
    long generatedEntries() const { return nGenerated_; }
    long motherEntries() const { return nMother_; }

    // Create the histogram pair for each dimension and fill the mother side from
    // `sourceFile`:`sourceTree`, selecting hits with the given VD id and pdgId (the same
    // accepted-hit selection the trainer and the PtotResampler use). `dir` is the
    // TFileDirectory the hists are created in — give each set its own subdirectory so
    // several (source, pdg) sets can coexist in one TFileService file. `ip` must be the
    // SAME transform parameters generation uses, otherwise mother and generated live in
    // different coordinate systems. `massMeV` is the species' rest mass, used only to turn
    // |p| into the kinetic energy the 2D line-zoom views plot; pass the module's own
    // pdt->particle(pdgId).mass() so the ParticleDataList stays its single source.
    // `stats` (optional, may be null) carries the model's per-slot training statistics in
    // TRANSFORM slot order; when given, the transformed axes are sized to that population
    // instead of the static fallback ranges.
    void book(const art::TFileDirectory& dir, const std::string& label,
              MomentumBasis basis, int pdgId, double massMeV,
              const std::string& sourceFile, const std::string& sourceTree,
              unsigned long virtualDetectorID, const InverseParams& ip,
              const std::string& moduleName,
              const TransformedStatsBySlot* stats = nullptr)
    {
        label_      = label;
        basis_      = basis;
        pdgId_      = pdgId;
        moduleName_ = moduleName;
        dir_        = std::make_unique<art::TFileDirectory>(dir);
        specs_      = perDimensionSpecs(basis, pdgId, stats);
        specs2D_    = correlationDimSpecs(basis, pdgId);
        // The 2D views' r and pr/pphi are measured about the SAME detector axis the
        // transform uses, so both sides share one origin and one mass.
        x0_   = ip.x0;
        y0_   = ip.y0;
        mass_ = massMeV;

        mother_.reserve(specs_.size());
        generated_.reserve(specs_.size());
        for (const DimSpec& s : specs_) {
            const std::string title = label + " " + s.suffix;
            mother_.push_back(dir.make<TH1D>((s.suffix + "_mother").c_str(),
                                             (title + " (mother)").c_str(),
                                             s.nbins, s.lo, s.hi));
            generated_.push_back(dir.make<TH1D>((s.suffix + "_generated").c_str(),
                                                (title + " (generated)").c_str(),
                                                s.nbins, s.lo, s.hi));
            mother_.back()->Sumw2();
            generated_.back()->Sumw2();
        }

        // 2D correlation views, one mother/generated pair each. Same naming convention as
        // the 1D pairs, so a drawing script finds them the same way.
        mother2D_.reserve(specs2D_.size());
        generated2D_.reserve(specs2D_.size());
        for (const Dim2DSpec& s : specs2D_) {
            const std::string title = label + " " + s.name;
            auto makeOne = [&](const char* role) {
                TH2D* h = dir.make<TH2D>((s.name + "_" + role).c_str(),
                                         (title + " (" + role + ")").c_str(),
                                         s.xbins, s.xlo, s.xhi, s.ybins, s.ylo, s.yhi);
                // Unlike the 1D hists, whose axis titles wait for the rebin factor, these are
                // never rebinned, so their titles are final at book time.
                h->GetXaxis()->SetTitle(s.xTitle.c_str());
                h->GetYaxis()->SetTitle(s.yTitle.c_str());
                h->Sumw2();
                return h;
            };
            mother2D_.push_back(makeOne("mother"));
            generated2D_.push_back(makeOne("generated"));
        }
        // Metric summary, one row per dimension, written by finalize().
        metricTree_ = dir.make<TTree>("validationMetrics",
                                      ("Generated-vs-mother agreement metrics: " + label).c_str());
        metricTree_->Branch("dimension",   &brDimension_);
        metricTree_->Branch("transformed", &brTransformed_, "transformed/O");
        metricTree_->Branch("pdgId",       &pdgId_,         "pdgId/I");
        metricTree_->Branch("W1",          &brW1_,          "W1/D");
        metricTree_->Branch("JSD",         &brJsd_,         "JSD/D");
        metricTree_->Branch("TV",          &brTv_,          "TV/D");
        metricTree_->Branch("KS",          &brKs_,          "KS/D");
        metricTree_->Branch("chi2",        &brChi2_,        "chi2/D");
        metricTree_->Branch("ndf",         &brNdf_,         "ndf/I");
        metricTree_->Branch("nbins",       &brNbins_,       "nbins/I");
        metricTree_->Branch("motherEntries",    &brMotherEntries_,    "motherEntries/L");
        metricTree_->Branch("generatedEntries", &brGeneratedEntries_, "generatedEntries/L");

        // Report which transformed axes were sized from the model statistics and which fell
        // back to the static range, with the resulting range — an axis that silently kept a
        // far-too-wide fallback is otherwise only discoverable by looking at the plots.
        if (stats) {
            std::ostringstream axes;
            axes << "ValidationPlots [" << label << "]: transformed axis ranges";
            for (const DimSpec& s : specs_) {
                if (!s.transformed) continue;
                const TransformedDimStats& st = (*stats)[s.slot];
                axes << "\n  " << s.suffix << ": [" << s.lo << ", " << s.hi << "] in "
                     << s.nbins << " bins";
                if (st.valid)
                    axes << "  (from training stats: mean=" << st.mean << " sigma=" << st.stdev
                         << " range=[" << st.min << ", " << st.max << "])";
                else
                    axes << "  (STATIC FALLBACK — no training statistics for this slot)";
            }
            mf::LogInfo(moduleName_) << axes.str();
        }

        fillMother(sourceFile, sourceTree, virtualDetectorID, ip);
        booked_ = true;
    }

    // Fill the generated side with one sample. `g` carries the transformed coordinates as
    // the model produced them; (x,y,t,px,py,pz) are the inverted physical values. Passing
    // both means the transformed plots show the model's own output rather than a
    // re-forward-transform of the inverse (which would hide any inversion asymmetry).
    // z is not taken: it is VDz0 for every sample and has no distribution.
    void fillGenerated(const GeneratedTransformed& g,
                       double x, double y, double t, double px, double py, double pz)
    {
        if (!booked_) return;
        // For the V2 resampler path m0 is set from the raw drawn pTotal by
        // generateTwoStage, so the transformed pTotal slot is meaningful in both paths.
        const double trans[6] = {g.xTrans, g.yTrans, g.tTrans, g.m0, g.m1, g.m2};
        const double phys[6]  = {px, py, pz, x, y, t};
        for (size_t i = 0; i < specs_.size(); ++i) {
            const DimSpec& s = specs_[i];
            generated_[i]->Fill(s.transformed ? trans[s.slot] : phys[s.slot]);
        }
        fill2D(generated2D_, x, y, t, px, py, pz);
        ++nGenerated_;
    }

    // Rebin both sides to a statistics-appropriate resolution, normalize to unit area,
    // compute every dimension's metrics, fill the summary tree and log the table.
    // Returns the per-dimension metrics in spec order; safe on an unbooked set (empty).
    std::vector<DistributionMetrics> finalize() {
        std::vector<DistributionMetrics> out;
        if (!booked_) return out;

        if (nGenerated_ == 0) {
            mf::LogWarning(moduleName_)
                << "ValidationPlots [" << label_ << "]: no generated samples were recorded; "
                << "the comparison is empty. (For GenerateMix this simply means the mix never "
                << "drew this source/particle.)";
            return out;
        }
        if (nGenerated_ < kLowStatisticsWarn)
            mf::LogWarning(moduleName_)
                << "ValidationPlots [" << label_ << "]: only " << nGenerated_
                << " generated sample(s) (< " << kLowStatisticsWarn << "); the metrics below are "
                << "dominated by sampling noise, not by model bias. Generate more events before "
                << "reading anything into them.";

        // Resolution is set by whichever side is thinner: fine bins the sparser sample
        // cannot populate only add noise to every metric. Both sides then take the SAME
        // integer factor, keeping their edges aligned.
        const long nMin = std::min(nMother_, nGenerated_);
        const int target = targetBinCount(nMin);

        // The whole table is accumulated into one stream and emitted as a SINGLE LogInfo:
        // a per-dimension message would stamp each line with its own timestamp/header and
        // scatter the rows, which is both noisier to read and slower to write.
        std::ostringstream report;
        report << "ValidationPlots [" << label_ << "] pdg " << pdgId_ << ", basis "
               << momentumBasisName(basis_) << ": generated " << nGenerated_
               << " sample(s) vs " << nMother_ << " mother hit(s); rebinning to ~" << target
               << " bins (driven by the smaller count, " << nMin << ").\n"
               << "Per-dimension agreement (0 is a perfect match for W1/JSD/TV/KS; "
               << "chi2/NDF ~ 1):";

        for (size_t i = 0; i < specs_.size(); ++i) {
            const DimSpec& s = specs_[i];
            const int factor = rebinFactor(s.nbins, target);
            if (factor > 1) {
                mother_[i]->Rebin(factor);
                generated_[i]->Rebin(factor);
            }
            // Axis titles are set only now: the y-axis states the ACHIEVED bin width, which
            // is not known until the rebin factor is. TH1::Rebin may drop a partial group of
            // trailing bins, so read the width back off the axis rather than recomputing it.
            const double width = mother_[i]->GetXaxis()->GetBinWidth(1);
            const std::string yTitle = binWidthLabel(width, s.binWidthUnit);
            for (TH1D* h : {mother_[i], generated_[i]}) {
                h->GetXaxis()->SetTitle(s.xTitle.c_str());
                h->GetYaxis()->SetTitle(yTitle.c_str());
            }

            // Unit-area normalization on BOTH sides: the mother file and the generated job
            // have unrelated statistics, so only the shapes are comparable.
            normalizeToUnitArea(*mother_[i]);
            normalizeToUnitArea(*generated_[i]);

            const DistributionMetrics m =
                computeDistributionMetrics(*generated_[i], *mother_[i], moduleName_);
            out.push_back(m);

            brDimension_        = s.suffix;
            brTransformed_      = s.transformed;
            brW1_               = m.w1;
            brJsd_              = m.jsd;
            brTv_               = m.tv;
            brKs_               = m.ks;
            brChi2_             = m.chi2;
            brNdf_              = m.ndf;
            brNbins_            = mother_[i]->GetNbinsX();
            brMotherEntries_    = nMother_;
            brGeneratedEntries_ = nGenerated_;
            metricTree_->Fill();

            report << "\n  " << (s.transformed ? "[trans] " : "[phys ] ") << s.suffix
                   << " (" << mother_[i]->GetNbinsX() << " bins): W1=" << m.w1
                   << " JSD=" << m.jsd << " TV=" << m.tv << " KS=" << m.ks
                   << " chi2/NDF=" << m.chi2 << "/" << m.ndf << "=" << m.chi2ndf();

            // Overlay + ratio canvas for THIS dimension. specs_ holds the transformed
            // coordinates followed by the physical ones, so every dimension of both families
            // gets one. Built here, after the rebin/normalize/metrics above, so the canvas
            // shows exactly the histograms the metrics were computed from.
            writeComparisonCanvas(s, *mother_[i], *generated_[i], m);
        }

        // The 2D views are normalized (to unit volume, on both sides) but not rebinned and
        // not scored: they are read by eye, and the rebin factor above is chosen for the
        // metrics of the 1D marginals, which do not apply here.
        //
        // Both sides of a pair then get the SAME z-range, so equal colour means equal
        // density when the two are put next to each other -- without this ROOT autoscales
        // each pad to its own max and the comparison is visually meaningless. Unlike a
        // palette (a draw-time gStyle global), SetMinimum/SetMaximum are STORED on the
        // histogram, so the matched range travels with the ROOT file to whatever draws it.
        // The floor is pinned at 0 rather than left to autoscale, so an empty bin reads as
        // the bottom of the scale in both plots.
        for (size_t i = 0; i < specs2D_.size(); ++i) {
            normalizeToUnitVolume(*mother2D_[i]);
            normalizeToUnitVolume(*generated2D_[i]);

            // GetMaximum() with no argument returns the largest bin content (the plotted
            // range only, which is what the colour scale has to cover).
            const double zmax = std::max(mother2D_[i]->GetMaximum(),
                                         generated2D_[i]->GetMaximum());
            for (TH2D* h : {mother2D_[i], generated2D_[i]}) {
                h->SetMinimum(0.0);
                // An all-empty pair would give zmax=0 and an inverted range; leave those to
                // ROOT rather than storing min==max.
                if (zmax > 0.0) h->SetMaximum(zmax);
            }
        }
        if (!specs2D_.empty())
            report << "\n  " << specs2D_.size() << " 2D correlation view(s) written "
                   << "(mother/generated pairs, unit-normalized, no metrics).";

        mf::LogInfo(moduleName_) << report.str();
        return out;
    }

private:
    // Fill every 2D view of one side with a single sample. The derived quantities are
    // computed once and shared across the views, and the same code serves both sides, so
    // mother and generated can never disagree on how a quantity was formed.
    void fill2D(std::vector<TH2D*>& hists,
                double x, double y, double t, double px, double py, double pz)
    {
        if (hists.empty()) return;
        const Derived2D d = computeDerived2D(x, y, px, py, pz, x0_, y0_, mass_);
        for (size_t i = 0; i < specs2D_.size(); ++i) {
            const Dim2DSpec& s = specs2D_[i];
            hists[i]->Fill(var2DValue(s.xVar, d, x, y, t, pz),
                           var2DValue(s.yVar, d, x, y, t, pz));
        }
    }

    // "counts per <width><unit> (normalized)", the width to 3 significant figures.
    static std::string binWidthLabel(double width, const std::string& unit) {
        char buf[128];
        snprintf(buf, sizeof(buf), "counts per %.3g%s (normalized)", width, unit.c_str());
        return std::string(buf);
    }

    // Dimensions drawn on a log-y scale: sharply peaked slope/angular distributions, and the
    // time and momentum distributions, which span several decades. On a linear scale the
    // peak flattens everything else into the axis and the tails -- where the model most often
    // disagrees -- become invisible.
    static bool useLogY(const std::string& suffix) {
        static const std::set<std::string> kLogY = {
            "prTrans", "pphiTrans", "urTrans", "uphiTrans",
            "urTransAsinh", "uphiTransAsinh",
            "t", "tTrans", "tTransAsinh",
            "px", "py", "pz"
        };
        return kLogY.count(suffix) > 0;
    }

    // Build the mother-vs-generated overlay canvas for one dimension and write it into the
    // plot set's directory.
    //
    // Layout is two stacked pads, matching the reference comparison plots:
    //   top    -- both distributions overlaid (normalized, so directly comparable), log-y for
    //             the wide-dynamic-range dimensions, with a legend naming the two.
    //   bottom -- generated/mother bin-by-bin ratio, flat at 1 for a perfect match, with the
    //             agreement metrics printed in the legend header so the number and the shape
    //             that produced it are read together.
    //
    // Called AFTER both sides have been rebinned and normalized, so the ratio divides
    // histograms that already share bin edges and total area.
    //
    // The clones are owned by the canvas (kCanDelete) rather than leaked or added to the
    // directory: the canvas is what gets written, and the stored hists must keep their own
    // identity in the file.
    void writeComparisonCanvas(const DimSpec& s, const TH1D& mother, const TH1D& generated,
                               const DistributionMetrics& m)
    {
        if (!dir_) return;

        // makeAndRegister (not make): a TCanvas does not attach itself to a TDirectory the
        // way a TH1 does, so registering is what puts it in the output file. It is then
        // written by TFileService at close -- calling Write() as well would store it twice.
        const std::string cname = "c_" + s.suffix;
        const std::string ctitle = label_ + " " + s.suffix;
        TCanvas* c = dir_->makeAndRegister<TCanvas>(cname.c_str(), ctitle.c_str());
        c->SetCanvasSize(800, 1000);
        c->Divide(1, 2);

        // --- top pad: the two distributions overlaid ---
        c->cd(1);
        gPad->SetLeftMargin(0.15);
        if (useLogY(s.suffix)) gPad->SetLogy(1);

        TH1D* hMother = static_cast<TH1D*>(mother.Clone((s.suffix + "_ov_mother").c_str()));
        TH1D* hGen    = static_cast<TH1D*>(generated.Clone((s.suffix + "_ov_generated").c_str()));
        for (TH1D* h : {hMother, hGen}) {
            h->SetDirectory(nullptr);   // owned by the canvas, not by the output directory
            h->SetTitle("");
            h->SetStats(0);
            h->SetMarkerStyle(1);
            h->SetBit(kCanDelete);
        }
        hMother->SetLineColor(kRed + 1);
        hMother->SetMarkerColor(kRed + 1);
        hGen->SetLineColor(kBlue - 4);
        hGen->SetMarkerColor(kBlue - 4);

        THStack* stack = new THStack(("hs_" + s.suffix).c_str(), ctitle.c_str());
        stack->Add(hMother);
        stack->Add(hGen);
        stack->Draw("nostack");
        stack->GetHistogram()->GetXaxis()->SetTitle(s.xTitle.c_str());
        stack->GetHistogram()->GetYaxis()->SetTitle(mother.GetYaxis()->GetTitle());
        stack->GetHistogram()->GetYaxis()->SetTitleOffset(1.2);
        stack->SetBit(kCanDelete);

        // Transparent and inside the frame, so the distribution behind stays readable.
        TLegend* leg = new TLegend(0.62, 0.72, 0.88, 0.88);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->AddEntry(hMother, "mother (MC truth)", "lp");
        leg->AddEntry(hGen,    "generated", "lp");
        leg->Draw("SAME");
        leg->SetBit(kCanDelete);

        // --- bottom pad: generated / mother ratio ---
        c->cd(2);
        gPad->SetLeftMargin(0.15);
        gPad->SetLogy(0);

        TH1D* ratio = static_cast<TH1D*>(generated.Clone((s.suffix + "_ratio").c_str()));
        ratio->SetDirectory(nullptr);
        ratio->Divide(&mother);
        ratio->SetTitle("");
        ratio->SetStats(0);
        ratio->SetLineColor(kBlue - 4);
        ratio->SetMarkerColor(kBlue - 4);
        ratio->SetMarkerStyle(1);
        ratio->GetXaxis()->SetTitle(s.xTitle.c_str());
        ratio->GetYaxis()->SetTitle("generated / mother");
        ratio->GetYaxis()->SetTitleOffset(1.2);
        // Ratios rarely exceed ~3; a fixed window keeps the flat-at-1 reference readable
        // instead of letting one wild low-statistics bin set the scale.
        ratio->SetMinimum(0.0);
        ratio->SetMaximum(3.0);
        ratio->Draw("HIST");
        ratio->SetBit(kCanDelete);

        // Metrics in the legend header. Green (TLatex #color[8] == kGreen+2) marks the value
        // each takes for a PERFECT match: W1/JSD/TV/KS are all distances (0 = identical),
        // chi2/NDF ~ 1. Showing the scale next to the number saves looking up which way is
        // good for each metric.
        static const char* const kMetricHeader =
            "W1[#color[8]{#bf{0}},#infty) / JSD[#color[8]{#bf{0}},1] / "
            "TV[#color[8]{#bf{0}},1] / KS[#color[8]{#bf{0}},1] / "
            "#chi^{2}/NDF(#color[8]{#bf{~1}})";
        TLegend* legRatio = new TLegend(0.16, 0.74, 0.89, 0.88, kMetricHeader);
        legRatio->SetFillStyle(0);
        legRatio->SetBorderSize(0);
        char entry[512];
        snprintf(entry, sizeof(entry), "%.3g / %.3f / %.3f / %.2g / %.2f",
                 m.w1, m.jsd, m.tv, m.ks, m.chi2ndf());
        legRatio->AddEntry(ratio, entry, "lp");
        legRatio->Draw("SAME");
        legRatio->SetBit(kCanDelete);

        c->Update();
    }

    // Fill the mother histograms from the source ROOT file, pushing every accepted hit
    // through the SAME forward transform training used, so the transformed mother hists
    // are in the model's own coordinates.
    void fillMother(const std::string& sourceFile, const std::string& sourceTree,
                    unsigned long virtualDetectorID, const InverseParams& ip)
    {
        PzFallbackStats pzStats;
        forEachAcceptedHitRoot(
            sourceFile, sourceTree, virtualDetectorID, pdgId_, moduleName_,
            [&](double x, double y, double z, double t, double px, double py, double pz) {
                double xTrans, yTrans, tTrans, m0, m1, m2;
                forwardTransformSample(x, y, z, t, px, py, pz,
                                       ip.x0, ip.y0, ip.t0, ip.tScale, ip.p0, ip.VDr, ip.VDz0,
                                       xTrans, yTrans, tTrans, m0, m1, m2, basis_, &pzStats);
                const double trans[6] = {xTrans, yTrans, tTrans, m0, m1, m2};
                // The mother's physical coordinates must be the ones the transform defines:
                // the model only ever saw hits extrapolated onto the VDz0 plane, so comparing
                // against the raw step positions would show an offset that is an artifact of
                // that projection rather than a model defect. Recover them by inverting the
                // just-computed transform, which is exactly what generation does.
                double xm, ym, zm, tm, pxm, pym, pzm;
                invertGeneratedSample(xTrans, yTrans, tTrans, m0, m1, m2,
                                      ip.x0, ip.y0, ip.t0, ip.tScale, ip.p0, ip.VDr, ip.VDz0,
                                      xm, ym, zm, tm, pxm, pym, pzm, basis_);
                const double phys[6] = {pxm, pym, pzm, xm, ym, tm};
                for (size_t i = 0; i < specs_.size(); ++i) {
                    const DimSpec& s = specs_[i];
                    mother_[i]->Fill(s.transformed ? trans[s.slot] : phys[s.slot]);
                }
                // The 2D views take the SAME inverted physical values as the 1D ones (see
                // above): comparing generated samples against the raw step quantities would
                // show the VDz0 projection rather than a model defect.
                fill2D(mother2D_, xm, ym, tm, pxm, pym, pzm);
                ++nMother_;
            });

        if (nMother_ == 0)
            throw cet::exception(moduleName_)
                << "ValidationPlots [" << label_ << "]: no accepted hits in " << sourceFile
                << ":" << sourceTree << " (VDid=" << virtualDetectorID << ", pdgId=" << pdgId_
                << ", pz>0), so there is no mother distribution to compare against.";

        if (pzStats.count > 0)
            mf::LogWarning(moduleName_)
                << "ValidationPlots [" << label_ << "]: pz fell below the safety floor for "
                << pzStats.count << " of " << nMother_ << " mother hit(s) while building the "
                << "transformed reference; those slopes were computed against the floor.";

        mf::LogInfo(moduleName_)
            << "ValidationPlots [" << label_ << "]: filled the mother distribution with "
            << nMother_ << " hit(s) from " << sourceFile << ":" << sourceTree << ".";
    }

    // Scale to unit integral (underflow/overflow excluded, matching the plotted range).
    static void normalizeToUnitArea(TH1D& h) {
        const double integral = h.Integral();
        if (integral > 0.0) h.Scale(1.0 / integral);
    }

    // The 2D counterpart, normalizing over the WHOLE sample rather than the plotted range:
    // the flow bins are included in the integral (0..n+1 on both axes), so a bin's content
    // is its fraction of ALL entries that side recorded. This matters most for the photon
    // line zooms, whose +/- 15 keV window holds a tiny and DIFFERENT fraction of each side's
    // sample -- normalizing per-window would divide the two by unrelated denominators and
    // silently rescale any real difference in how much yield each side puts in the line.
    static void normalizeToUnitVolume(TH2D& h) {
        const double integral = h.Integral(0, h.GetNbinsX() + 1, 0, h.GetNbinsY() + 1);
        if (integral > 0.0) h.Scale(1.0 / integral);
    }

    bool booked_ = false;
    std::string label_;
    std::string moduleName_ = "VDResamplerValidationPlots";
    MomentumBasis basis_ = MomentumBasis::V1_CylindricalTransformed;
    int pdgId_ = 0;

    std::vector<DimSpec> specs_;
    // Owned by the TFileService directory they were made in, which is why these are raw
    // pointers and are never deleted here.
    std::vector<TH1D*> mother_;
    std::vector<TH1D*> generated_;

    // 2D correlation views, in specs2D_ order; owned by the TFileDirectory as above.
    std::vector<Dim2DSpec> specs2D_;
    std::vector<TH2D*> mother2D_;
    std::vector<TH2D*> generated2D_;

    // The directory the plot set was booked in, kept so finalize() can write the overlay
    // canvases into the same place. art::TFileDirectory has no default constructor, hence
    // the unique_ptr rather than a bare member.
    std::unique_ptr<art::TFileDirectory> dir_;
    // Detector-axis origin and rest mass the derived 2D quantities are formed with; taken
    // from the SAME InverseParams the transform uses and the module's ParticleDataList
    // mass, so both sides agree.
    double x0_ = 0.0, y0_ = 0.0, mass_ = 0.0;

    TTree* metricTree_ = nullptr;

    long nMother_ = 0;
    long nGenerated_ = 0;

    // metricTree_ branch buffers
    std::string brDimension_;
    bool   brTransformed_ = false;
    double brW1_ = 0.0, brJsd_ = 0.0, brTv_ = 0.0, brKs_ = 0.0, brChi2_ = 0.0;
    int    brNdf_ = 0, brNbins_ = 0;
    long   brMotherEntries_ = 0, brGeneratedEntries_ = 0;
};

} // namespace VDResampler
} // namespace mu2e
