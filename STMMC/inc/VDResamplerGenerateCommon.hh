#pragma once

// Shared generation + inversion logic for the VDResampler generate modules
// (VDResamplerGenerateFromModel and VDResamplerGenerateMix). Both modules hold the
// same kinds of pieces — an all-at-once model OR a stage1+stage2 pair, an optional
// pTotal resampler, and sampler settings — but organize ownership differently
// (member models vs. per-ParticleModelSet). These free functions take those pieces by
// parameter so both modules share the V1/V2 + diffusion/resampler logic.
//
// Header-only. Yongyi Wu, Jun. 2026

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "cetlib_except/exception.h"

#include "Offline/MachineLearningTools/inc/ScoreBasedDiffusionModel.hh"
#include "Offline/STMMC/inc/VDResamplerTransforms.hh"
#include "Offline/STMMC/inc/VDResamplerPtotResampler.hh"

namespace mu2e {
namespace VDResampler {

// Reverse-sampler settings shared by every generateSample call in a job.
struct SamplerSettings {
    bool   useEMANetworkIfAvailable = true;
    bool   useHeun                  = true;
    bool   useSDE                   = true;
    int    diffusionSteps           = 200;
    double sdeToOdeSigmaThreshold   = -1.0;
};

// Result of one transformed generation: (t,x,y) in slots 0..2 and the momentum triple
// (m0,m1,m2) in the model's basis. For the V2 resampler path, rawPtot holds the RAW
// physical pTotal (MeV/c) and rawPtotValid is true, so the caller inverts momentum from
// it directly (no log/exp round-trip).
struct GeneratedTransformed {
    double xTrans = 0.0, yTrans = 0.0, tTrans = 0.0;
    double m0 = 0.0, m1 = 0.0, m2 = 0.0;
    MomentumBasis basis = MomentumBasis::V1_CylindricalTransformed;
    PositionBasis posBasis = PositionBasis::V1_Atanh;
    double rawPtot = 0.0;
    bool   rawPtotValid = false;
};

// Detector-center / scale parameters for the inverse transform (same as training).
struct InverseParams {
    double x0, y0, t0, tScale, p0, VDr, VDz0;
};

// ---------------------------------------------------------------------------
// Training-data statistics for one TRANSFORMED slot, in raw (pre-z-score) transformed
// units — the units the validation histograms of the model's own coordinates are binned in.
// This is ScoreBasedDiffusionModel::DimStats carried across in a form the plotting header
// can use without depending on the model class.
//
// `valid` is false for a slot whose statistics are unavailable: the model predates the
// stored normalization, or (V2 two-stage) the slot is supplied by the pTotal resampler
// rather than a diffusion model and so has no dimStats to report. Consumers fall back to a
// static axis range for such a slot.
// ---------------------------------------------------------------------------
struct TransformedDimStats {
    bool   valid = false;
    double mean = 0.0, stdev = 0.0, min = 0.0, max = 0.0;
};

// Statistics for the six transformed slots in TRANSFORM order (xTrans, yTrans, tTrans, m0,
// m1, m2) — NOT the models' own (t,x,y,...) data order. collectTransformedStats below does
// that re-ordering.
using TransformedStatsBySlot = std::array<TransformedDimStats, 6>;

// ---------------------------------------------------------------------------
// generateAllAtOnce — sample the 6-D all-at-once model; basis from its tag.
// ---------------------------------------------------------------------------
inline GeneratedTransformed generateAllAtOnce(
    ScoreBasedDiffusionModel& model, const SamplerSettings& s, const std::string& moduleName)
{
    GeneratedTransformed g;
    g.basis    = unpackMomentumBasis(model.basisTag());
    g.posBasis = unpackPositionBasis(model.basisTag());
    const SBDMGeneratedSample sample = model.generateSample(
        {}, s.useEMANetworkIfAvailable, s.useHeun, s.useSDE, s.diffusionSteps, s.sdeToOdeSigmaThreshold);
    if (sample.zscore.size() != 6u)
        throw cet::exception(moduleName)
            << "All-at-once model returned " << sample.zscore.size() << " values, expected 6.";
    g.tTrans = sample.value[0]; g.xTrans = sample.value[1]; g.yTrans = sample.value[2];
    g.m0 = sample.value[3]; g.m1 = sample.value[4]; g.m2 = sample.value[5];
    return g;
}

// ---------------------------------------------------------------------------
// generateTwoStage — sample the two-stage pair. Basis is read from stage2's tag.
//   V1: stage1 {t,x,y}; stage2 {pr,pphi,pz} | {t,x,y}.
//   V2: stage1 {pTotal} (diffusion model OR resampler); stage2 {t,x,y,ur,uphi} | {pTotal}.
//   stage1Model may be null iff stage1Method != DIFFUSION (resampler supplies pTotal).
//   resampler is used only on the V2 non-diffusion path.
// ---------------------------------------------------------------------------
inline GeneratedTransformed generateTwoStage(
    ScoreBasedDiffusionModel* stage1Model, ScoreBasedDiffusionModel& stage2Model,
    Stage1Method stage1Method, const PtotResampler& resampler,
    CLHEP::RandFlat& rf, CLHEP::RandGaussQ& rg,
    const SamplerSettings& s, double p0, const std::string& moduleName)
{
    GeneratedTransformed g;
    g.basis    = unpackMomentumBasis(stage2Model.basisTag());
    g.posBasis = unpackPositionBasis(stage2Model.basisTag());
    const bool isV2 = (g.basis != MomentumBasis::V1_CylindricalTransformed);

    if (!isV2) {
        if (!stage1Model)
            throw cet::exception(moduleName) << "V1 two-stage requires a stage-1 model.";
        const SBDMGeneratedSample s1 = stage1Model->generateSample(
            {}, s.useEMANetworkIfAvailable, s.useHeun, s.useSDE, s.diffusionSteps, s.sdeToOdeSigmaThreshold);
        if (s1.zscore.size() != 3u)
            throw cet::exception(moduleName)
                << "V1 stage-1 model returned " << s1.zscore.size() << " values, expected 3.";
        g.tTrans = s1.value[0]; g.xTrans = s1.value[1]; g.yTrans = s1.value[2];
        const std::vector<double> cond = {s1.zscore[0], s1.zscore[1], s1.zscore[2]};
        const SBDMGeneratedSample s2 = stage2Model.generateSample(
            cond, s.useEMANetworkIfAvailable, s.useHeun, s.useSDE, s.diffusionSteps, s.sdeToOdeSigmaThreshold);
        if (s2.zscore.size() != 3u)
            throw cet::exception(moduleName)
                << "V1 stage-2 model returned " << s2.zscore.size() << " values, expected 3.";
        g.m0 = s2.value[0]; g.m1 = s2.value[1]; g.m2 = s2.value[2];
        return g;
    }

    // V2: obtain log(pTotal/p0) from a diffusion stage-1 model or the resampler. The
    // stage2 condition is z-scored via stage2's stored condition normalization, so it
    // matches training regardless of the pTotal source.
    double ptot_t;
    if (stage1Method == Stage1Method::DIFFUSION) {
        if (!stage1Model)
            throw cet::exception(moduleName)
                << "V2 two-stage with DIFFUSION stage-1 requires a stage-1 model.";
        const SBDMGeneratedSample s1 = stage1Model->generateSample(
            {}, s.useEMANetworkIfAvailable, s.useHeun, s.useSDE, s.diffusionSteps, s.sdeToOdeSigmaThreshold);
        if (s1.zscore.size() != 1u)
            throw cet::exception(moduleName)
                << "V2 stage-1 model returned " << s1.zscore.size() << " values, expected 1.";
        ptot_t = s1.value[0];
    } else {
        g.rawPtot      = resampler.draw(rf, rg);
        g.rawPtotValid = true;
        ptot_t         = std::log(std::max(g.rawPtot, kRadiusSafetyEpsilon) / p0);
    }
    g.m0 = ptot_t; // log(pTotal/p0) carried in slot m0
    const std::vector<double> cond = {stage2Model.normalizeCondition(0, ptot_t)};
    const SBDMGeneratedSample s2 = stage2Model.generateSample(
        cond, s.useEMANetworkIfAvailable, s.useHeun, s.useSDE, s.diffusionSteps, s.sdeToOdeSigmaThreshold);
    if (s2.zscore.size() != 5u)
        throw cet::exception(moduleName)
            << "V2 stage-2 model returned " << s2.zscore.size() << " values, expected 5.";
    // stage2 data layout: (t,x,y,ur,uphi).
    g.tTrans = s2.value[0]; g.xTrans = s2.value[1]; g.yTrans = s2.value[2];
    g.m1 = s2.value[3]; g.m2 = s2.value[4];
    return g;
}

// ---------------------------------------------------------------------------
// invertGenerated — invert one transformed sample to detector coordinates.
//   For the V2 resampler path (rawPtotValid) the RAW pTotal drives the momentum
//   reconstruction directly; otherwise the basis-dispatching wrapper is used.
// ---------------------------------------------------------------------------
inline void invertGenerated(
    const GeneratedTransformed& g, const InverseParams& ip,
    double& x, double& y, double& z, double& t, double& px, double& py, double& pz)
{
    if (g.rawPtotValid) {
        const bool asinhSlopes = basisUsesAsinhSlopes(g.basis);
        double dx, dy, r, ur, uphi;
        invertPosition(g.xTrans, g.yTrans, ip.x0, ip.y0, ip.VDr, g.posBasis, x, y, dx, dy, r);
        z = ip.VDz0;
        t = invertTimeForBasis(g.tTrans, ip.t0, ip.tScale, g.basis);
        v2DecodeSlopes(g.m1, g.m2, asinhSlopes, ur, uphi);
        invertMomentumV2FromPtot(g.rawPtot, ur, uphi, dx, dy, r, px, py, pz);
    } else {
        invertGeneratedSample(
            g.xTrans, g.yTrans, g.tTrans, g.m0, g.m1, g.m2,
            ip.x0, ip.y0, ip.t0, ip.tScale, ip.p0, ip.VDr, ip.VDz0,
            x, y, z, t, px, py, pz, g.basis, g.posBasis);
    }
}

// ---------------------------------------------------------------------------
// collectTransformedStats — the models' per-dimension training statistics, re-ordered from
//   each model's own data layout into TRANSFORM slot order (xTrans, yTrans, tTrans, m0, m1,
//   m2) for the validation plots.
//
// The re-ordering is the whole point of this helper. A model's data vector is
// (t,x,y,...)-ordered while the transform slots are (x,y,t,...)-ordered, and which model
// supplies which slot depends on the layout:
//
//   AllAtOnce6D          dims (t,x,y,m0,m1,m2)  -> slots (1,2,0,3,4,5) from the one model
//   V1 two-stage         stage1 (t,x,y)         -> slots (1,2,0)
//                        stage2 (m0,m1,m2)      -> slots (3,4,5)
//   V2/V3 two-stage      stage2 (t,x,y,m1,m2)   -> slots (1,2,0,4,5)
//                        stage1 (log pTotal)    -> slot 3, ONLY when stage-1 is a diffusion
//                                                  model; with a pTotal resampler that slot
//                                                  has no statistics and stays invalid.
//
// A slot left invalid keeps the validation plots' static fallback axis, so a missing or
// pre-normalization model degrades to the previous behaviour rather than failing.
//
// `stage1Model` may be null (V2 resampler path, or a caller with no stage-1 at all).
inline TransformedStatsBySlot collectTransformedStats(
    const ScoreBasedDiffusionModel* allAtOnceModel,
    const ScoreBasedDiffusionModel* stage1Model,
    const ScoreBasedDiffusionModel* stage2Model,
    MomentumBasis basis)
{
    TransformedStatsBySlot out;

    // Copy one model dimension into one transform slot.
    auto take = [&out](const ScoreBasedDiffusionModel* model, int modelDim, int slot) {
        if (!model || !model->hasDataNormalization()) return;
        const ScoreBasedDiffusionModel::DimStats s = model->dimStats(modelDim);
        out[slot] = TransformedDimStats{true, s.mean, s.stdev, s.min, s.max};
    };

    if (allAtOnceModel) {
        // (t,x,y,m0,m1,m2) -> slots (2,0,1,3,4,5): model dim d holds transform slot map[d].
        static constexpr int kMap[6] = {2, 0, 1, 3, 4, 5};
        for (int d = 0; d < 6; ++d) take(allAtOnceModel, d, kMap[d]);
        return out;
    }

    if (basis == MomentumBasis::V1_CylindricalTransformed) {
        // stage1 (t,x,y) -> slots (2,0,1); stage2 (m0,m1,m2) -> slots (3,4,5).
        static constexpr int kStage1Map[3] = {2, 0, 1};
        for (int d = 0; d < 3; ++d) take(stage1Model, d, kStage1Map[d]);
        for (int d = 0; d < 3; ++d) take(stage2Model, d, 3 + d);
        return out;
    }

    // V2/V3: stage2 is (t,x,y,ur,uphi) -> slots (2,0,1,4,5).
    static constexpr int kStage2Map[5] = {2, 0, 1, 4, 5};
    for (int d = 0; d < 5; ++d) take(stage2Model, d, kStage2Map[d]);
    // stage-1 is the 1-D log(pTotal/p0) model when diffusion supplies it; when the pTotal
    // resampler does, stage1Model is null and slot 3 stays invalid by construction.
    take(stage1Model, 0, 3);
    return out;
}

} // namespace VDResampler
} // namespace mu2e
