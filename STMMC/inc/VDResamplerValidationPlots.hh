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
// Header-only; depends on ROOT (TH1D/TTree/TDirectory) and the VDResampler transforms.
// Yongyi Wu, Aug. 2026

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

#include "art_root_io/TFileDirectory.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1D.h"
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
//   rebinFactor: the integer TH1::Rebin group size bringing `nbins` closest to `target`.
//     TH1::Rebin merges an integer number of bins, so the achieved count only
//     approximates the target.
// ---------------------------------------------------------------------------
constexpr int kMinTargetBins = 50;
constexpr int kMaxTargetBins = 4000;

inline int targetBinCount(long n) {
    if (n <= 0) return kMinTargetBins;
    const double t = 2.0 * std::sqrt(static_cast<double>(n));
    return std::max(kMinTargetBins, std::min(kMaxTargetBins, static_cast<int>(std::lround(t))));
}

inline int rebinFactor(int nbins, int target) {
    if (target <= 0) return 1;
    return std::max(1, static_cast<int>(std::lround(static_cast<double>(nbins) / target)));
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

// The six transformed dimensions for `basis`. Slot order is
// (xTrans, yTrans, tTrans, m0, m1, m2) — the TRANSFORM's order, not the model's internal
// (t,x,y,...) data layout.
inline std::vector<DimSpec> transformedDimSpecs(MomentumBasis basis) {
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
    return specs;
}

// The six recovered physical quantities. The momentum window is widened for neutrons,
// whose spectrum extends well past the few-MeV range that covers the other species.
inline std::vector<DimSpec> physicalDimSpecs(int pdgId) {
    const bool isNeutron = (pdgId == 2112);
    const double pmax  = isNeutron ? 50. : 8.;
    const int    pnbin = 1000;

    std::vector<DimSpec> specs;
    specs.push_back({"px", "px [MeV]", " MeV", 2 * pnbin, -pmax, pmax, false, kPx});
    specs.push_back({"py", "py [MeV]", " MeV", 2 * pnbin, -pmax, pmax, false, kPy});
    specs.push_back({"pz", "pz [MeV]", " MeV", pnbin, 0., pmax, false, kPz});
    specs.push_back({"x",  "x [mm]", "mm", 800, -3904. - 2000., -3904. + 2000., false, kX});
    specs.push_back({"y",  "y [mm]", "mm", 800, -2000., 2000., false, kY});
    specs.push_back({"t",  "t [ns]", "ns", 500, 0., 5000., false, kT});
    return specs;
}

// Full per-dimension comparison list: transformed coordinates + physical quantities.
inline std::vector<DimSpec> perDimensionSpecs(MomentumBasis basis, int pdgId) {
    std::vector<DimSpec> specs = transformedDimSpecs(basis);
    const std::vector<DimSpec> phys = physicalDimSpecs(pdgId);
    specs.insert(specs.end(), phys.begin(), phys.end());
    return specs;
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
    // different coordinate systems.
    void book(const art::TFileDirectory& dir, const std::string& label,
              MomentumBasis basis, int pdgId,
              const std::string& sourceFile, const std::string& sourceTree,
              unsigned long virtualDetectorID, const InverseParams& ip,
              const std::string& moduleName)
    {
        label_      = label;
        basis_      = basis;
        pdgId_      = pdgId;
        moduleName_ = moduleName;
        specs_      = perDimensionSpecs(basis, pdgId);

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
        }

        mf::LogInfo(moduleName_) << report.str();
        return out;
    }

private:
    // "counts per <width><unit> (normalized)". %g trims trailing zeros.
    static std::string binWidthLabel(double width, const std::string& unit) {
        char buf[128];
        snprintf(buf, sizeof(buf), "counts per %g%s (normalized)", width, unit.c_str());
        return std::string(buf);
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
