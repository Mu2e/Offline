#pragma once

// Non-diffusion 1-D resampler for the total momentum |p| (pTotal, MeV/c), used as
// the V2 two-stage "stage 1". The empirical distribution IS the source data, so this
// needs no training and no saved model: at generation time the source ROOT file is
// read once, physical pTotal is extracted per accepted hit, and samples are drawn
// from the empirical distribution. Operating on PHYSICAL pTotal (not the transformed
// log(pTotal/p0)) is where the fine energy structure is defined and where sub-keV
// resolution is meaningful; the caller applies log/z-score afterward only to build
// the stage-2 conditioning input, while the RAW drawn pTotal is carried directly to
// the V2 inversion (no log/exp round-trip).
//
// Header-only; depends on ROOT (TFile/TTree, TSpline3) and CLHEP randoms.
// Yongyi Wu, Jun. 2026

#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TFile.h"
#include "TTree.h"
#include "TSpline.h"

namespace mu2e {
namespace VDResampler {

// ---------------------------------------------------------------------------
// Stage-1 generation method for pTotal.
//   DIFFUSION   : NOT handled here — a trained 1-D SBDM is used instead.
//   INVERSE_CDF : sorted order-statistics inverse-CDF (no binning; sub-keV). Default.
//   SPLINE_CDF  : monotone cubic (TSpline3) through the empirical CDF, inverted.
//   KDE         : Gaussian-kernel smoothing (blurs sharp lines; not recommended
//                 for sub-keV structure, provided as an option).
// ---------------------------------------------------------------------------
enum class Stage1Method { DIFFUSION = 0, INVERSE_CDF = 1, SPLINE_CDF = 2, KDE = 3 };

inline Stage1Method parseStage1Method(const std::string& m, const std::string& moduleName) {
    if (m == "DIFFUSION")   return Stage1Method::DIFFUSION;
    if (m == "INVERSE_CDF") return Stage1Method::INVERSE_CDF;
    if (m == "SPLINE_CDF")  return Stage1Method::SPLINE_CDF;
    if (m == "KDE")         return Stage1Method::KDE;
    throw cet::exception(moduleName)
        << "Unrecognized SBDMstage1Method value \"" << m
        << "\" (expected DIFFUSION / INVERSE_CDF / SPLINE_CDF / KDE).";
}

// ---------------------------------------------------------------------------
// forEachAcceptedHitRoot — open `file`/`tree` and invoke cb(x,y,z,t,px,py,pz) for
//   every hit passing the standard VD-resampler selection (matching VD id, pdgId
//   (0 = any), pz>0). Mirrors the read loop in VDResamplerTrainFromRoot so the
//   resampler sees the SAME accepted-hit set the trainer did. Branch schema:
//   double time,x,y,z,px,py,pz; int pdgId; ULong64_t virtualdetectorId.
// ---------------------------------------------------------------------------
inline void forEachAcceptedHitRoot(
    const std::string& file, const std::string& tree,
    unsigned long virtualDetectorID, int pdgID, const std::string& moduleName,
    const std::function<void(double,double,double,double,double,double,double)>& cb)
{
    TFile fin(file.c_str(), "READ");
    if (!fin.IsOpen())
        throw cet::exception(moduleName) << "Cannot open ROOT file: " << file;
    TTree* ttree = dynamic_cast<TTree*>(fin.Get(tree.c_str()));
    if (!ttree)
        throw cet::exception(moduleName) << "Cannot find TTree: " << tree;

    double time, x, y, z, px, py, pz;
    int stepPdgId;
    ULong64_t vdId;
    ttree->SetBranchAddress("time",              &time);
    ttree->SetBranchAddress("x",                 &x);
    ttree->SetBranchAddress("y",                 &y);
    ttree->SetBranchAddress("z",                 &z);
    ttree->SetBranchAddress("px",                &px);
    ttree->SetBranchAddress("py",                &py);
    ttree->SetBranchAddress("pz",                &pz);
    ttree->SetBranchAddress("pdgId",             &stepPdgId);
    ttree->SetBranchAddress("virtualdetectorId", &vdId);

    for (Long64_t i = 0; i < ttree->GetEntries(); ++i) {
        ttree->GetEntry(i);
        if (vdId != virtualDetectorID || (stepPdgId != pdgID && pdgID != 0) || pz <= 0)
            continue;
        cb(x, y, z, time, px, py, pz);
    }
    fin.Close();
}

// ---------------------------------------------------------------------------
// PtotResampler — holds the sorted physical pTotal samples and draws from them.
//   Built in memory each generation job (no serialization). For SPLINE_CDF/KDE the
//   same sorted array is kept and the spline / kernel sum is evaluated on demand.
// ---------------------------------------------------------------------------
class PtotResampler {
public:
    // Below this, the empirical distribution is sparsely sampled — a generic
    // small-dataset warning (statistics of any method will be noisy), independent of
    // the basis/transform resolution.
    static constexpr std::size_t kSmallSampleWarn = 1000;

    // Bound on redraws when SPLINE_CDF / KDE produce a non-positive pTotal (physical
    // magnitude must be > 0). Exceeding this points at a malformed spline / oversized
    // KDE bandwidth rather than an unlucky draw, so we throw instead of looping forever.
    static constexpr int kMaxNonPositiveRedraws = 100;

    // Populate from a source ROOT file using the standard selection. After this the
    // sorted physical pTotal vector is ready; method-specific structures are built lazily.
    void buildFromRoot(const std::string& file, const std::string& tree,
                       unsigned long virtualDetectorID, int pdgID,
                       Stage1Method method, const std::string& moduleName)
    {
        method_ = method;
        ptot_.clear();
        forEachAcceptedHitRoot(file, tree, virtualDetectorID, pdgID, moduleName,
            [&](double, double, double, double, double px, double py, double pz) {
                ptot_.push_back(std::sqrt(px * px + py * py + pz * pz));
            });
        if (ptot_.empty())
            throw cet::exception(moduleName)
                << "PtotResampler: no accepted hits in " << file << ":" << tree
                << " (VDid=" << virtualDetectorID << ", pdgId=" << pdgID << ", pz>0).";
        std::sort(ptot_.begin(), ptot_.end());

        // Report the accumulated sample size and, when small, an estimate of the
        // achievable resolution. Median adjacent-sample spacing is a rough proxy for
        // the inverse-CDF's local resolution (a draw interpolates between neighbors).
        const std::size_t N = ptot_.size();
        mf::LogInfo(moduleName)
            << "PtotResampler: accumulated " << N << " physical pTotal sample(s) from "
            << file << ":" << tree << " (range [" << ptot_.front() << ", " << ptot_.back()
            << "] MeV/c).";
        if (N < kSmallSampleWarn) {
            mf::LogWarning(moduleName)
                << "PtotResampler: only " << N << " samples (< " << kSmallSampleWarn
                << ") — a small source dataset; resampled statistics will be noisy. "
                << "Estimated median adjacent-sample spacing ~ " << medianSpacing()
                << " MeV/c (roughly the finest structure INVERSE_CDF can resolve). "
                << "Consider a larger source file.";
        }

        if (method_ == Stage1Method::SPLINE_CDF) buildSpline(moduleName);
        // KDE bandwidth: Silverman's rule on the sorted data (used only for KDE draws).
        if (method_ == Stage1Method::KDE) bandwidth_ = silvermanBandwidth();
    }

    std::size_t size() const { return ptot_.size(); }

    // Draw one physical pTotal (MeV/c) via the configured method.
    double draw(CLHEP::RandFlat& rf, CLHEP::RandGaussQ& rg) const {
        switch (method_) {
            case Stage1Method::SPLINE_CDF:  return drawSpline(rf);
            case Stage1Method::KDE:         return drawKde(rf, rg);
            case Stage1Method::INVERSE_CDF: return drawInverseCdf(rf);
            case Stage1Method::DIFFUSION:
            default:
                throw cet::exception("PtotResampler")
                    << "draw() called with DIFFUSION method (use the stage-1 diffusion model instead).";
        }
    }

private:
    Stage1Method method_ = Stage1Method::INVERSE_CDF;
    std::vector<double> ptot_;            // sorted physical pTotal (MeV/c)
    std::unique_ptr<TSpline3> cdfSpline_; // SPLINE_CDF: pTotal as a function of CDF value
    double bandwidth_ = 0.0;              // KDE bandwidth

    // Quantile (pTotal value) at CDF level c in [0,1] by continuous-rank linear
    // interpolation between order statistics. With rank q = c*(N-1), c=0 -> ptot_[0]
    // and c=1 -> ptot_[N-1]; intermediate c interpolates the two bracketing samples.
    // No grid quantization, so local resolution = inter-sample spacing.
    double quantile(double c) const {
        const std::size_t N = ptot_.size();
        if (N == 1) return ptot_[0];
        const double q = c * static_cast<double>(N - 1);
        const std::size_t i = static_cast<std::size_t>(std::floor(q));
        if (i + 1 >= N) return ptot_.back(); // only reached at c==1 exactly; no pileup for c<1
        const double frac = q - static_cast<double>(i);
        return ptot_[i] + frac * (ptot_[i + 1] - ptot_[i]);
    }

    // INVERSE_CDF: draw u~U[0,1) (RandFlat::fire is half-open, so c<1 and the c==1
    // branch in quantile() never fires -> no artificial pileup at the max).
    double drawInverseCdf(CLHEP::RandFlat& rf) const {
        return quantile(rf.fire());
    }

    // Knot count for the CDF spline. The spline is only a smoothing convenience (the
    // sub-keV-faithful path is INVERSE_CDF, which has NO binning); knots merely need to
    // be dense enough that the spline tracks the empirical CDF's shape without
    // overfitting per-sample noise. Heuristic: ~1 knot per kSamplesPerKnot samples
    // (so each spline segment averages over that many points), floored for a smooth
    // curve on small datasets and capped to bound the spline size.
    int splineKnotCount() const {
        constexpr std::size_t kSamplesPerKnot = 100;  // ~averaging window per segment
        constexpr int kMinKnots = 50;                 // floor: smoothness on small N
        constexpr int kMaxKnots = 100000;             // cap: bound TSpline3 size
        const std::size_t N = ptot_.size();
        const int byDensity = static_cast<int>(N / kSamplesPerKnot);
        int n = std::max(kMinKnots, std::min(kMaxKnots, byDensity));
        return std::min<int>(n, static_cast<int>(N)); // never more knots than samples
    }

    // SPLINE_CDF: invert a monotone-cubic fit of pTotal vs the empirical CDF level.
    // Knots at evenly spaced CDF levels; each value is the corresponding quantile.
    void buildSpline(const std::string& moduleName) {
        const int nKnots = splineKnotCount();
        mf::LogInfo(moduleName)
            << "PtotResampler: SPLINE_CDF using " << nKnots << " knots over "
            << ptot_.size() << " samples.";
        std::vector<double> cdf(nKnots), val(nKnots);
        for (int k = 0; k < nKnots; ++k) {
            const double c = (nKnots == 1) ? 0.5 : static_cast<double>(k) / (nKnots - 1);
            cdf[k] = c;
            val[k] = quantile(c);
        }
        cdfSpline_ = std::make_unique<TSpline3>("ptotCdfSpline", cdf.data(), val.data(), nKnots);
    }

    // Cubic-spline inversion can extrapolate below 0 near the CDF ends; pTotal is a
    // physical magnitude, so reject a non-positive draw and try again (bounded).
    double drawSpline(CLHEP::RandFlat& rf) const {
        for (int attempt = 0; attempt < kMaxNonPositiveRedraws; ++attempt) {
            const double p = cdfSpline_->Eval(rf.fire());
            if (p > 0.0) return p;
            mf::LogWarning("PtotResampler")
                << "SPLINE_CDF drew non-positive pTotal " << p << " MeV/c (spline extrapolation); redrawing.";
        }
        throw cet::exception("PtotResampler")
            << "SPLINE_CDF failed to draw a positive pTotal after " << kMaxNonPositiveRedraws
            << " attempts; the spline fit is likely producing negative values.";
    }

    // KDE: pick a sample at random and jitter by a Gaussian of width bandwidth_. The
    // Gaussian jitter can push a small-pTotal sample below 0; reject and redraw (bounded).
    double drawKde(CLHEP::RandFlat& rf, CLHEP::RandGaussQ& rg) const {
        for (int attempt = 0; attempt < kMaxNonPositiveRedraws; ++attempt) {
            const std::size_t i = static_cast<std::size_t>(rf.fireInt(static_cast<long>(ptot_.size())));
            const double p = ptot_[i] + bandwidth_ * rg.fire();
            if (p > 0.0) return p;
            mf::LogWarning("PtotResampler")
                << "KDE drew non-positive pTotal " << p << " MeV/c (Gaussian jitter); redrawing.";
        }
        throw cet::exception("PtotResampler")
            << "KDE failed to draw a positive pTotal after " << kMaxNonPositiveRedraws
            << " attempts; the bandwidth is likely too large for the smallest pTotal samples.";
    }

    // Median spacing between adjacent sorted samples — a rough resolution proxy.
    double medianSpacing() const {
        const std::size_t N = ptot_.size();
        if (N < 2) return 0.0;
        std::vector<double> gaps;
        gaps.reserve(N - 1);
        for (std::size_t i = 1; i < N; ++i) gaps.push_back(ptot_[i] - ptot_[i - 1]);
        std::nth_element(gaps.begin(), gaps.begin() + gaps.size() / 2, gaps.end());
        return gaps[gaps.size() / 2];
    }

    double silvermanBandwidth() const {
        const std::size_t N = ptot_.size();
        if (N < 2) return 0.0;
        double mean = 0.0;
        for (double v : ptot_) mean += v;
        mean /= static_cast<double>(N);
        double var = 0.0;
        for (double v : ptot_) var += (v - mean) * (v - mean);
        var /= static_cast<double>(N - 1);
        const double sigma = std::sqrt(var);
        return 1.06 * sigma * std::pow(static_cast<double>(N), -0.2);
    }
};

} // namespace VDResampler
} // namespace mu2e
