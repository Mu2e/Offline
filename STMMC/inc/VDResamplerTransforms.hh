#pragma once

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace mu2e {
  namespace VDResampler {

    // nominal center of the detector
    constexpr double kX0 = -3904.0;
    constexpr double kY0 = 0.0;
    // sensitive time scale
    constexpr double kT0 = 1700.0; // ns
    constexpr double kTScale = 1.0;
    // tunable momentum scale
    constexpr double kP0 = 1.0; // MeV/c

    // safety constants for numerical stability in the forward and inverse transforms.
    // NOTE kRadiusSafetyEpsilon has several unrelated uses (an r>eps test in mm before
    // dividing by r to build the local polar frame, and a floor on pTot in MeV before
    // log), so it is deliberately NOT reused as the radial clamp -- see
    // kRhoClampEpsilon below, which is dimensionless and guards u = r/VDr.
    constexpr double kRadiusSafetyEpsilon = 1e-6;
    constexpr double kMinSafeTime = 0.1;

    // Radial clamp: rho = r/VDr is capped at 1-kRhoClampEpsilon before the radial map,
    // so an event exactly on (or numerically past) the rim cannot produce u = inf. This
    // is a pure numerical guard. It is NOT a fix for the radial mis-modelling described
    // under PositionBasis: no generated event ever reached the clamp (0/100000
    // measured), so that was not rim saturation, and tightening or loosening eps does
    // not address it. Raising it only truncates real area near the rim (the truncated
    // fraction is ~2 eps for a uniform disc), so it is kept at the historical value.
    //
    // Both PositionBasis maps are atanh-based and so diverge at rho = 1; the clamp is
    // what keeps u finite there.
    constexpr double kRhoClampEpsilon = 1e-6;
    // Independent floor for pz in divisions (slope basis). Distinct from
    // kRadiusSafetyEpsilon because it guards a different quantity (longitudinal
    // momentum, not transverse radius). If it is ever used, the input pz was
    // pathologically small or non-positive (pz>0 is expected from the selection).
    constexpr double kPzSafetyEpsilon = 1e-9;
    // asinh slope scales for the V2 asinh variant (slopes are near 0 for a forward
    // beam, so asinh ~ identity there; this only compresses large-|u| tails).
    // when additional kSlopeScales < 1 are used, the regions near u=0 are stretched 
    // and the tails are more compressed.
    // Separate scales for the radial (ur=pr/pz) and azimuthal (uphi=pphi/pz) slopes
    // so each can be tuned to its own spread.
    constexpr double kUrSlopeScale   = 0.05;
    constexpr double kUphiSlopeScale = 0.05;

    // Tail-taming of the log-time coordinate for the V3 basis. After the shared
    // forwardTime (ln(t/t0)/tScale), V3 applies asinh((base - kTBulkCenter)/kTTailScale)
    // to compress a heavy but physical right tail (~1e-6 of events reach z-score ~100
    // after the log alone). kTBulkCenter is the empirical bulk center of ln(t/t0)/tScale
    // (measured ~-2.25 across particle species) so the bulk lands at asinh(0)=0 where the
    // map is locally linear/symmetric and only the far tail is compressed. asinh is smooth,
    // monotone, and exactly invertible with no hard cutoff (unlike asinh(pz/p0), which was
    // rejected for its wall at 0). kTTailScale sets the asinh width; smaller => tails
    // compressed harder, bulk near center stretched more.
    constexpr double kTBulkCenter = -2.25;
    constexpr double kTTailScale  = 1.0;

    // ------------------------------------------------------------------------
    // PzFallbackStats — accumulates pz-fallback occurrences across a transform
    //   loop so the CALLER can emit a SINGLE summary warning at the end (instead
    //   of one line per hit). This struct performs NO I/O so the header stays
    //   free of messagefacility/iostream: it only records. The caller (a module
    //   with mf::LogWarning available) inspects count/firstValues after its loop.
    // ------------------------------------------------------------------------
    struct PzFallbackStats {
      static constexpr std::size_t kMaxSamples = 50;
      std::size_t count = 0;                 // total fallbacks
      std::vector<double> firstValues;       // first kMaxSamples offending pz values

      void record(double pz) {
        ++count;
        if (firstValues.size() < kMaxSamples) firstValues.push_back(pz);
      }
    };

    // ------------------------------------------------------------------------
    // MomentumBasis — selects the momentum transform math used by the ALL-AT-ONCE
    //   6-vector (t,x,y, m0,m1,m2). It defines what the three momentum slots mean
    //   and therefore which forward/inverse math applies. (For two-stage models the
    //   per-stage layout is given by ModelLayout below; the V2 stages reuse the V2
    //   slope/pTotal math.)
    //   V1_CylindricalTransformed : original basis. A per-component transform of
    //       the LOCAL CYLINDRICAL momentum (pr, pphi, pz):
    //       (asinh(pr/p0), asinh(pphi/p0), log(pz/p0)). NOTE this is NOT an exact
    //       direction/magnitude factorization — it transforms the cylindrical
    //       components independently. Default, for backward compatibility.
    //       Momentum slot order: m0=asinh(pr/p0), m1=asinh(pphi/p0), m2=log(pz/p0).
    //   V2_PtotSlopes : pTotal carries the energy structure; ur,uphi are transverse
    //       slopes (theta=0 regular, unbounded). Inversion is exact and enforces
    //       |p|=pTotal, pz>0. Momentum slot order: m0=log(pTotal/p0), m1=ur=pr/pz,
    //       m2=uphi=pphi/pz  (pTotal FIRST).
    //   V2_PtotSlopesAsinh : as V2_PtotSlopes but slopes (m1,m2) wrapped in
    //       asinh(ur/kUrSlopeScale), asinh(uphi/kUphiSlopeScale) to tame heavy wide-angle tails.
    //   V3_PtotSlopesAsinhTimeAsinh : identical momentum treatment to V2_PtotSlopesAsinh,
    //       plus asinh tail-taming of the log-time coordinate (see kTBulkCenter/kTTailScale).
    //       Use V2_PtotSlopesAsinh when the extra time transform is NOT wanted.
    // In the all-at-once vector t,x,y occupy slots 0,1,2 and momentum occupies 3,4,5.
    // ------------------------------------------------------------------------
    enum class MomentumBasis {
      V1_CylindricalTransformed   = 0,
      V2_PtotSlopes               = 1,
      V2_PtotSlopesAsinh          = 2,
      V3_PtotSlopesAsinhTimeAsinh = 3
    };

    // ------------------------------------------------------------------------
    // PositionBasis - the radial boundary-removing map u(rho), rho = r/VDr in [0,1),
    //   selected PER SPECIES and orthogonal to MomentumBasis. Every variant keeps the
    //   SAME polar encoding (xTrans,yTrans) = u*(cos theta, sin theta) with
    //   theta = atan2(dy,dx), so the angular treatment is untouched; only u(rho)
    //   differs. All are monotone, send rho -> [0,inf), and invert in closed form.
    //
    // Three effects decide whether a map suits a species. They are easy to conflate, so
    // state them explicitly:
    //
    //   EFFECT 1 - RESOLUTION, drho/du. How far a fixed error in the network's output u
    //     displaces the physical radius. Where drho/du is large, a given u-error costs
    //     more in rho. This is about precision, and it is the one that matters when a
    //     species' structure is concentrated somewhere specific.
    //
    //   EFFECT 2 - DENSITY, p(u) = p(rho) drho/du. What the training distribution
    //     actually looks like in the coordinate the network fits. A density that is flat
    //     and compactly supported is easy; one with a thin tail is not, because the tail
    //     is sparsely sampled and a score model generically over-populates it. This is
    //     about how many samples land where.
    //
    //     For a uniform disc p(rho) ~ rho, so BOTH maps give p(u) -> 0 at the origin:
    //     that is inherited from the geometry, not a property of either map, and cannot
    //     be fixed by choosing between them. The maps differ in their TAILS.
    //
    //   EFFECT 3 - CORE WIDTH AFTER NORMALIZATION. The model z-scores each coordinate,
    //     so what it fits is u divided by the stdev of u over the whole training set. A
    //     map whose u grows steeply toward the rim inflates that stdev, and dividing by
    //     it squeezes the populated core into a narrow band -- which shows up directly
    //     as a generated peak that is too narrow. This is what ruled out the
    //     rho/(1-rho) style maps, whose u reaches 99 at rho=0.99. It is invisible in u
    //     itself and only appears after normalization, so it is easily missed when maps
    //     are compared on their drho/du alone.
    //
    // A caution from an earlier investigation. A radial mis-modelling seen under
    //   V1_Atanh (generated/mother ratio ~0.7 at u~0, peaking ~1.38 at u~1.7, back to
    //   ~0.7 by u~2.8 -- a deficit at BOTH ends with an excess in the middle) looked
    //   like a coordinate problem but was not: that shape is a distribution contracted
    //   toward its own mode, the signature of a weak score field. Those runs were later
    //   found to have stopped with a training loss still above 0.99, i.e. undertrained.
    //   So before reaching for this enum to explain a radial discrepancy, check the
    //   training loss and the planner's stop reason first.
    //
    //   V1_Atanh     : u = atanh(rho); rho = tanh(u); drho/du = 1-rho^2.
    //       Effect 1: drho/du = 1 at rho=0, 0.75 at rho=0.5, 0.19 at rho=0.9. Coarsest
    //         at the centre, finest at the rim. Since u ~ rho + O(rho^3) near the
    //         origin, it gives no extra ABSOLUTE resolution where the events actually
    //         are for a centrally-peaked species.
    //       Effect 2: for a uniform disc p(u) = tanh(u) sech^2(u), peaking near u~0.66
    //         with an e^{-2u} tail -- thin, and hard to terminate correctly.
    //   V2_AtanhSqrt : u = atanh(sqrt(rho)); rho = tanh^2(u);
    //         drho/du = 2 sqrt(rho) (1-rho).
    //       Effect 1: the sqrt sends drho/du -> 0 as rho -> 0, so du/drho DIVERGES at
    //         the centre. Where atanh spends a fixed amount of u range per unit rho near
    //         the origin, this spends an unbounded amount: the inner region is STRETCHED
    //         rather than compressed. Away from the centre it converges back toward
    //         atanh -- 0.71 at rho=0.5, 0.19 at rho=0.9 -- so the rim keeps the same
    //         tanh saturation and the same class of tail.
    //       Effect 3: this is the reason to prefer it. Output range stays modest
    //         (u = 1.47 at rho=0.9, 2.65 at rho=0.99), so the stdev is not inflated by a
    //         steep rim and the core survives the z-score at a usable width.
    // ------------------------------------------------------------------------
    enum class PositionBasis {
      V1_Atanh     = 0,
      V2_AtanhSqrt = 1
    };

    // ------------------------------------------------------------------------
    // PeakTag — one monoenergetic line, named by its position on the RAW physical pTotal
    //   axis. A hit is "in" the line when |pTot - center| <= halfWidth, both in MeV/c.
    //
    // WHY A DISCRETE TAG. The stage-2 condition is m0 = log(pTotal/p0), z-scored across the
    // whole spectrum. For a line at 1.809 MeV a +/-1 keV window is Delta log p ~ 5.5e-4,
    // which after the z-score is well under 1e-3 sigma. Making p(t | pTot) change character
    // across a gap that narrow would demand an effective slope of order 1e3 in the condition;
    // a Fourier condition embedding could only reach that with a base frequency ~1e3, which
    // would alias everywhere else in the spectrum. The conditional distribution really does
    // differ — a prompt line photon has a different arrival-time shape from the degraded
    // continuum at a neighbouring energy — so the clean fix is to hand the network the class
    // membership directly instead of asking it to infer it from a near-discontinuity.
    //
    // The tag is a DETERMINISTIC FUNCTION of pTotal, which is what makes it safe: pTotal is
    // already drawn (by stage 1 or the resampler) before stage 2 runs, so the tag is available
    // at generation time with no extra sampling, adds no new marginal to model, and cannot
    // desynchronize from the condition it was derived from.
    //
    // CAVEAT worth knowing when reading generated output. The tag only fires correctly if the
    // pTotal feeding it reproduces the line. A DIFFUSION stage-1 over log(pTotal/p0) smears a
    // 1 keV line, so the tagged fraction will not match the truth. The empirical pTotal
    // resampler (INVERSE_CDF / SPLINE_CDF) preserves the line exactly and is the configuration
    // this feature is meant for.
    // ------------------------------------------------------------------------
    struct PeakTag {
      double center    = 0.0; // MeV/c, RAW pTotal (not transformed)
      double halfWidth = 0.0; // MeV/c, inclusive half-window
    };

    // Index of the class label within the V2 stage-2 condition vector, which is
    // {log(pTotal/p0)} or {log(pTotal/p0), label}. The label goes LAST so the existing
    // coordinate keeps index 0 and every normalizeCondition(0, ...) call means what it always
    // meant. Defined here so the train and generate sides index the same slot by name rather
    // than by a repeated literal.
    constexpr int kPeakTagCondIdx = 1;

    // Class label for one hit: 0 = in no configured line (continuum), k+1 = inside tags[k].
    // Zero-based tag index k maps to label k+1 so that 0 always means "none", which keeps the
    // label meaningful when no tags are configured at all (every hit is 0).
    //
    // Windows are expected to be disjoint; if they do overlap the FIRST match wins, so the
    // label stays single-valued and the ordering of the configured list decides. Callers that
    // care validate disjointness up front (see assemblePeakTags).
    inline int peakTagFor(double pTot, const std::vector<PeakTag>& tags) {
      for (std::size_t k = 0; k < tags.size(); ++k)
        if (std::abs(pTot - tags[k].center) <= tags[k].halfWidth)
          return static_cast<int>(k) + 1;
      return 0;
    }

    // ------------------------------------------------------------------------
    // ModelLayout — the role/vector layout of a single saved model, so the
    //   generate side knows how to interpret and invert it from the file alone.
    //   AllAtOnce6D          : one model over (t,x,y, m0,m1,m2) per MomentumBasis.
    //   TwoStageStage1Ptot1D : V2 two-stage stage 1, a 1-D model over (log(pTotal/p0)).
    //   TwoStageStage2_5D    : V2 two-stage stage 2, a 5-D model over
    //                          (t,x,y,ur,uphi) CONDITIONED on log(pTotal/p0).
    // ------------------------------------------------------------------------
    enum class ModelLayout {
      AllAtOnce6D          = 0,
      TwoStageStage1Ptot1D = 1,
      TwoStageStage2_5D    = 2
    };

    // ------------------------------------------------------------------------
    // Opaque basis-tag packing. The SBDM stores a single int32 it never interprets;
    // VDResampler packs (ModelLayout, PositionBasis, MomentumBasis) into it and
    // unpacks on load. Encoding: peakTagged*1000 + layout*100 + position*10 + momentum,
    // one decade per field, which caps ModelLayout, PositionBasis and MomentumBasis at
    // 10 entries each. Widen the strides if that is ever reached.
    //
    // This stays backward compatible with the pre-PositionBasis encoding
    // (layout*100 + momentum): an old tag has momentum < 10 in the units place and
    // nothing in the tens, so it decodes to the same layout and momentum with
    // position = 0 = V1_Atanh, which IS the map those models were trained with. A
    // pre-tag (v<=6) model still loads as tag 0 = (AllAtOnce6D, V1_Atanh,
    // V1_CylindricalTransformed), preserving the original behaviour.
    //
    // The thousands digit is the peak-tag flag: 1 when the model carries the extra
    // categorical condition dim described by PeakTag, 0 otherwise. Every tag written
    // before that feature has nothing in the thousands place and so decodes to 0 = not
    // tagged, again matching how those models were trained. The flag is what lets the
    // generate side decide, FROM THE CHECKPOINT ALONE, whether stage 2 wants a one- or
    // two-element condition vector, so a tagged model and an untagged one cannot be fed
    // the wrong shape.
    // ------------------------------------------------------------------------
    constexpr int kBasisTagPeakStride     = 1000;
    constexpr int kBasisTagLayoutStride   = 100;
    constexpr int kBasisTagPositionStride = 10;

    inline int packBasisTag(ModelLayout layout, PositionBasis position, MomentumBasis basis,
                            bool peakTagged = false) {
      return (peakTagged ? kBasisTagPeakStride : 0)
           + static_cast<int>(layout)   * kBasisTagLayoutStride
           + static_cast<int>(position) * kBasisTagPositionStride
           + static_cast<int>(basis);
    }
    inline ModelLayout unpackModelLayout(int tag) {
      return static_cast<ModelLayout>((tag % kBasisTagPeakStride) / kBasisTagLayoutStride);
    }
    inline PositionBasis unpackPositionBasis(int tag) {
      return static_cast<PositionBasis>((tag % kBasisTagLayoutStride) / kBasisTagPositionStride);
    }
    inline MomentumBasis unpackMomentumBasis(int tag) {
      return static_cast<MomentumBasis>(tag % kBasisTagPositionStride);
    }
    // True when the model carries the peak-tag condition dim (see PeakTag).
    inline bool unpackPeakTagged(int tag) {
      return (tag / kBasisTagPeakStride) != 0;
    }

    // Human-readable enum names and a full basisTag decode, shared by the train and
    // generate modules so an opaque tag is never printed raw. Pure string helpers —
    // they add no I/O dependency, keeping this header messagefacility-free.
    inline const char* modelLayoutName(ModelLayout l) {
      switch (l) {
        case ModelLayout::AllAtOnce6D:          return "AllAtOnce6D";
        case ModelLayout::TwoStageStage1Ptot1D: return "TwoStageStage1Ptot1D";
        case ModelLayout::TwoStageStage2_5D:    return "TwoStageStage2_5D";
      }
      return "unknown";
    }
    inline const char* positionBasisName(PositionBasis p) {
      switch (p) {
        case PositionBasis::V1_Atanh:     return "V1_Atanh";
        case PositionBasis::V2_AtanhSqrt: return "V2_AtanhSqrt";
      }
      return "unknown";
    }
    inline const char* momentumBasisName(MomentumBasis b) {
      switch (b) {
        case MomentumBasis::V1_CylindricalTransformed: return "V1_CylindricalTransformed";
        case MomentumBasis::V2_PtotSlopes:             return "V2_PtotSlopes";
        case MomentumBasis::V2_PtotSlopesAsinh:        return "V2_PtotSlopesAsinh";
        case MomentumBasis::V3_PtotSlopesAsinhTimeAsinh: return "V3_PtotSlopesAsinhTimeAsinh";
      }
      return "unknown";
    }
    // Decode an opaque basisTag (peakTagged*1000 + layout*100 + position*10 + momentum)
    // into "layout=<...>, basis=<...> (tag <n>)" for logs and error messages.
    inline std::string basisTagToString(int tag) {
      std::ostringstream os;
      os << "layout=" << modelLayoutName(unpackModelLayout(tag))
         << ", position=" << positionBasisName(unpackPositionBasis(tag))
         << ", basis=" << momentumBasisName(unpackMomentumBasis(tag))
         << ", peakTagged=" << (unpackPeakTagged(tag) ? "yes" : "no")
         << " (tag " << tag << ")";
      return os.str();
    }

    // ------------------------------------------------------------------------
    // Shared geometry helpers (basis-independent)
    // ------------------------------------------------------------------------

    // Extrapolate (x,y) along the momentum to the nominal VDz0 plane, then shift to
    // detector-centered coordinates. Outputs dx, dy and the in-plane radius r.
    inline void extrapolateAndCenter(
      double x, double y, double z, double px, double py, double pz,
      double x0, double y0, double VDz0,
      double& dx, double& dy, double& r)
    {
      const double extrapolationFactor = (VDz0 - z) / pz;
      const double xExtrapolated = x + extrapolationFactor * px;
      const double yExtrapolated = y + extrapolationFactor * py;
      dx = xExtrapolated - x0;
      dy = yExtrapolated - y0;
      r  = std::sqrt(dx * dx + dy * dy);
    }

    // The radial map u(rho) and its inverse, selected by PositionBasis. Split out from
    // forwardPosition/invertPosition so the two directions cannot drift apart: they are
    // exact inverses for every enum value, which invertPosition relies on. Also used by
    // the validation plots to build the transformed radial coordinate.
    inline double radialForward(double rho, PositionBasis basis) {
      switch (basis) {
        case PositionBasis::V2_AtanhSqrt: {
          const double s = std::sqrt(rho);
          return 0.5 * std::log((1.0 + s) / (1.0 - s)); // atanh(sqrt(rho))
        }
        case PositionBasis::V1_Atanh:
        default:
          return 0.5 * std::log((1.0 + rho) / (1.0 - rho)); // atanh(rho)
      }
    }
    inline double radialInverse(double u, PositionBasis basis) {
      switch (basis) {
        case PositionBasis::V2_AtanhSqrt: {
          const double th = std::tanh(u);
          return th * th;                                  // tanh^2(u)
        }
        case PositionBasis::V1_Atanh:
        default:
          return std::tanh(u);
      }
    }

    // Forward position transform: (dx,dy,r) -> (xTrans,yTrans) = u(rho)*(cos,sin theta)
    // with rho = r/VDr. The angular part is identical for every PositionBasis; only the
    // radial map differs (see PositionBasis for which to use and why).
    inline void forwardPosition(double dx, double dy, double r, double VDr,
                                PositionBasis basis,
                                double& xTrans, double& yTrans)
    {
      double rho = r / VDr;
      // Numerical guard only, so rho=1 cannot produce u=inf; see kRhoClampEpsilon.
      rho = std::min(rho, 1.0 - kRhoClampEpsilon);
      const double u = radialForward(rho, basis);
      const double theta = std::atan2(dy, dx);
      xTrans = u * std::cos(theta);
      yTrans = u * std::sin(theta);
    }

    // Inverse position transform: (xTrans,yTrans) -> detector-space (x,y) and the local
    // frame (dx,dy,r) needed to rotate momentum back. MUST be given the same
    // PositionBasis the forward transform used; radialInverse is its exact inverse.
    //
    // Deliberately NOT clamped to the forward map's kRhoClampEpsilon range: every
    // radialInverse already returns rho<1 by construction, and a generated u above the
    // forward u_max is better left to land smoothly within the last fraction of a mm
    // than snapped onto a hard edge, which would build a rim spike. The generated-vs-
    // source radial ratio near rho=1 is the diagnostic for whether that tail matters.
    inline void invertPosition(double xTrans, double yTrans, double x0, double y0,
                               double VDr, PositionBasis basis,
                               double& x, double& y, double& dx, double& dy, double& r)
    {
      const double u = std::sqrt(xTrans * xTrans + yTrans * yTrans);
      const double theta = std::atan2(yTrans, xTrans);
      const double rho = radialInverse(u, basis);
      r  = rho * VDr;
      dx = r * std::cos(theta);
      dy = r * std::sin(theta);
      x  = dx + x0;
      y  = dy + y0;
    }

    // Project Cartesian (px,py) onto the LOCAL cylindrical frame defined by (dx,dy).
    // Returns radial (pr) and azimuthal (pphi) momentum. Universal across bases.
    inline void cartesianToLocalPolar(double px, double py, double dx, double dy, double r,
                                      double& pr, double& pphi)
    {
      if (r > kRadiusSafetyEpsilon) { // else r~0: pr~px, pphi~py
        const double rx = dx / r, ry = dy / r;     // radial unit vector
        const double phix = -ry,  phiy = rx;       // azimuthal unit vector
        pr   = px * rx + py * ry;
        pphi = px * phix + py * phiy;
      } else {
        pr   = px;
        pphi = py;
      }
    }

    // Rotate local cylindrical (pr,pphi) back to Cartesian (px,py) using (dx,dy).
    // Inverse of cartesianToLocalPolar. Universal across bases.
    inline void localPolarToCartesian(double pr, double pphi, double dx, double dy, double r,
                                      double& px, double& py)
    {
      if (r > kRadiusSafetyEpsilon) {
        const double rx = dx / r, ry = dy / r;
        const double phix = -ry,  phiy = rx;
        px = pr * rx + pphi * phix;
        py = pr * ry + pphi * phiy;
      } else {
        px = pr;
        py = pphi;
      }
    }

    // Forward time transform (shared).
    inline double forwardTime(double t, double t0, double tScale) {
      const double tSafe = (t > kMinSafeTime) ? t : kMinSafeTime; // avoid log(0)
      return std::log(tSafe / t0) / tScale;
    }

    // Inverse time transform (shared).
    inline double invertTime(double tTrans, double t0, double tScale) {
      return t0 * std::exp(tTrans * tScale);
    }

    // V3 time transform: asinh tail-taming layered on top of the shared log-time.
    // Centered at kTBulkCenter so the bulk maps near asinh(0)=0. Exactly invertible.
    inline double forwardTimeAsinh(double t, double t0, double tScale) {
      const double base = forwardTime(t, t0, tScale); // ln(tSafe/t0)/tScale
      return std::asinh((base - kTBulkCenter) / kTTailScale);
    }
    inline double invertTimeAsinh(double tTrans, double t0, double tScale) {
      const double base = kTBulkCenter + kTTailScale * std::sinh(tTrans);
      return t0 * std::exp(base * tScale);
    }

    // Basis feature predicates — single source of truth for "which knobs does this
    // basis turn on", so callers never enumerate basis values themselves. A basis is
    // added here once and every dispatch site follows.
    inline bool basisUsesAsinhSlopes(MomentumBasis b) {
      return b == MomentumBasis::V2_PtotSlopesAsinh
          || b == MomentumBasis::V3_PtotSlopesAsinhTimeAsinh;
    }
    inline bool basisUsesAsinhTime(MomentumBasis b) {
      return b == MomentumBasis::V3_PtotSlopesAsinhTimeAsinh;
    }

    // Basis-aware time transforms: dispatch to the asinh variant for bases that
    // request it, else the plain log-time. Use these wherever a bare forwardTime/
    // invertTime would otherwise be called against a known basis.
    inline double forwardTimeForBasis(double t, double t0, double tScale, MomentumBasis b) {
      return basisUsesAsinhTime(b) ? forwardTimeAsinh(t, t0, tScale)
                                   : forwardTime(t, t0, tScale);
    }
    inline double invertTimeForBasis(double tTrans, double t0, double tScale, MomentumBasis b) {
      return basisUsesAsinhTime(b) ? invertTimeAsinh(tTrans, t0, tScale)
                                   : invertTime(tTrans, t0, tScale);
    }

    // Shared momentum reconstruction for V2 given a RAW physical pTotal (MeV/c) and
    // the (already de-asinh'd) slopes ur,uphi. Lets the resampler path pass its raw
    // drawn pTotal directly (no log/exp round-trip); also used by invertGeneratedSampleV2.
    inline void invertMomentumV2FromPtot(
      double pTot, double ur, double uphi, double dx, double dy, double r,
      double& px, double& py, double& pz)
    {
      pz = pTot / std::sqrt(1.0 + ur * ur + uphi * uphi); // >0; denom >= 1, no guard needed
      const double pr   = ur   * pz;
      const double pphi = uphi * pz;
      localPolarToCartesian(pr, pphi, dx, dy, r, px, py);
    }

    // ========================================================================
    // V1_CylindricalTransformed  (original basis; default)
    //   momentum slots: m0=asinh(pr/p0), m1=asinh(pphi/p0), m2=log(pz/p0)
    // ========================================================================
    inline void forwardTransformSampleV1(
      const double x, const double y, const double z, const double t,
      const double px, const double py, const double pz,
      const double x0, const double y0, const double t0,
      const double tScale, const double p0, const double VDr, const double VDz0,
      double& xTrans, double& yTrans, double& tTrans,
      double& prTrans, double& pphiTrans, double& pzTrans,
      const PositionBasis posBasis = PositionBasis::V1_Atanh)
    {
      double dx, dy, r;
      extrapolateAndCenter(x, y, z, px, py, pz, x0, y0, VDz0, dx, dy, r);
      forwardPosition(dx, dy, r, VDr, posBasis, xTrans, yTrans);

      double pr, pphi;
      cartesianToLocalPolar(px, py, dx, dy, r, pr, pphi);
      prTrans   = std::asinh(pr / p0);
      pphiTrans = std::asinh(pphi / p0);
      pzTrans   = std::log(pz / p0); // tried asinh(pz/p0) but hard cutoff at 0 was unfriendly for DM

      tTrans = forwardTime(t, t0, tScale);
    }

    inline void invertGeneratedSampleV1(
      const double xTrans, const double yTrans, const double tTrans,
      const double prTrans, const double pphiTrans, const double pzTrans,
      const double x0, const double y0, const double t0,
      const double tScale, const double p0, const double VDr, const double VDz0,
      double& x, double& y, double& z, double& t,
      double& px, double& py, double& pz,
      const PositionBasis posBasis = PositionBasis::V1_Atanh)
    {
      double dx, dy, r;
      invertPosition(xTrans, yTrans, x0, y0, VDr, posBasis, x, y, dx, dy, r);
      z = VDz0;
      t = invertTime(tTrans, t0, tScale);

      const double pr   = p0 * std::sinh(prTrans);
      const double pphi = p0 * std::sinh(pphiTrans);
      pz = p0 * std::exp(pzTrans);
      localPolarToCartesian(pr, pphi, dx, dy, r, px, py);
    }

    // ========================================================================
    // V2_PtotSlopes / V2_PtotSlopesAsinh  (all-at-once 6-vector)
    //   momentum slots: m0=log(pTotal/p0), m1=ur=pr/pz, m2=uphi=pphi/pz
    //   asinh variant wraps the slopes: m1=asinh(ur/kUrSlopeScale), m2=asinh(uphi/kUphiSlopeScale)
    //   Inversion: pz = pTotal/sqrt(1+ur^2+uphi^2) (>0), pr=ur*pz, pphi=uphi*pz.
    //   pzStats (optional): records pz-fallback occurrences for a single summary
    //   warning by the caller (see PzFallbackStats). The inverse needs no such
    //   guard — its denominator sqrt(1+ur^2+uphi^2) >= 1.
    // ========================================================================
    inline void forwardTransformSampleV2(
      const double x, const double y, const double z, const double t,
      const double px, const double py, const double pz,
      const double x0, const double y0, const double t0,
      const double tScale, const double p0, const double VDr, const double VDz0,
      double& xTrans, double& yTrans, double& tTrans,
      double& pTotTrans, double& urTrans, double& uphiTrans,
      const bool asinhSlopes, const bool asinhTime = false,
      PzFallbackStats* pzStats = nullptr,
      const PositionBasis posBasis = PositionBasis::V1_Atanh)
    {
      double dx, dy, r;
      extrapolateAndCenter(x, y, z, px, py, pz, x0, y0, VDz0, dx, dy, r);
      forwardPosition(dx, dy, r, VDr, posBasis, xTrans, yTrans);

      double pr, pphi;
      cartesianToLocalPolar(px, py, dx, dy, r, pr, pphi);

      double pzSafe = pz;
      if (std::abs(pz) < kPzSafetyEpsilon) {
        if (pzStats) pzStats->record(pz);
        pzSafe = kPzSafetyEpsilon;
      }
      double ur   = pr   / pzSafe;
      double uphi = pphi / pzSafe;
      if (asinhSlopes) {
        ur   = std::asinh(ur   / kUrSlopeScale);
        uphi = std::asinh(uphi / kUphiSlopeScale);
      }

      const double pTot = std::sqrt(px * px + py * py + pz * pz);
      pTotTrans = std::log(std::max(pTot, kRadiusSafetyEpsilon) / p0);
      urTrans   = ur;
      uphiTrans = uphi;

      tTrans = asinhTime ? forwardTimeAsinh(t, t0, tScale) : forwardTime(t, t0, tScale);
    }

    inline void invertGeneratedSampleV2(
      const double xTrans, const double yTrans, const double tTrans,
      const double pTotTrans, const double urTrans, const double uphiTrans,
      const double x0, const double y0, const double t0,
      const double tScale, const double p0, const double VDr, const double VDz0,
      double& x, double& y, double& z, double& t,
      double& px, double& py, double& pz,
      const bool asinhSlopes, const bool asinhTime = false,
      const PositionBasis posBasis = PositionBasis::V1_Atanh)
    {
      double dx, dy, r;
      invertPosition(xTrans, yTrans, x0, y0, VDr, posBasis, x, y, dx, dy, r);
      z = VDz0;
      t = asinhTime ? invertTimeAsinh(tTrans, t0, tScale) : invertTime(tTrans, t0, tScale);

      double ur = urTrans, uphi = uphiTrans;
      if (asinhSlopes) {
        ur   = kUrSlopeScale   * std::sinh(ur);
        uphi = kUphiSlopeScale * std::sinh(uphi);
      }
      const double pTot = p0 * std::exp(pTotTrans);
      invertMomentumV2FromPtot(pTot, ur, uphi, dx, dy, r, px, py, pz);
    }

    // ------------------------------------------------------------------------
    // V2 slope helpers (de-asinh / asinh of the slope pair) — used by the two-stage
    // generate path which assembles (pTotal, ur, uphi) from stage1 + stage2 and the
    // raw resampled pTotal, then calls invertMomentumV2FromPtot directly.
    // ------------------------------------------------------------------------
    inline void v2DecodeSlopes(double m1, double m2, bool asinhSlopes,
                               double& ur, double& uphi)
    {
      if (asinhSlopes) { ur = kUrSlopeScale * std::sinh(m1); uphi = kUphiSlopeScale * std::sinh(m2); }
      else             { ur = m1;                            uphi = m2;                             }
    }

    // ========================================================================
    // Dispatching wrappers for the ALL-AT-ONCE 6-vector. Default basis =
    // V1_CylindricalTransformed so existing callers (which omit the basis argument)
    // are unchanged / backward compatible. Momentum out-params are basis-neutral
    // here (m0,m1,m2). (Two-stage models assemble the V2 vector explicitly and use
    // invertMomentumV2FromPtot / v2DecodeSlopes rather than these wrappers.)
    // ========================================================================
    inline void forwardTransformSample(
      const double x, const double y, const double z, const double t,
      const double px, const double py, const double pz,
      const double x0, const double y0, const double t0,
      const double tScale, const double p0, const double VDr, const double VDz0,
      double& xTrans, double& yTrans, double& tTrans,
      double& m0, double& m1, double& m2,
      const MomentumBasis basis = MomentumBasis::V1_CylindricalTransformed,
      PzFallbackStats* pzStats = nullptr,
      const PositionBasis posBasis = PositionBasis::V1_Atanh)
    {
      switch (basis) {
        case MomentumBasis::V2_PtotSlopes:
          forwardTransformSampleV2(x, y, z, t, px, py, pz, x0, y0, t0, tScale, p0, VDr, VDz0,
                                   xTrans, yTrans, tTrans, m0, m1, m2,
                                   /*asinhSlopes=*/false, /*asinhTime=*/false, pzStats, posBasis);
          break;
        case MomentumBasis::V2_PtotSlopesAsinh:
          forwardTransformSampleV2(x, y, z, t, px, py, pz, x0, y0, t0, tScale, p0, VDr, VDz0,
                                   xTrans, yTrans, tTrans, m0, m1, m2,
                                   /*asinhSlopes=*/true, /*asinhTime=*/false, pzStats, posBasis);
          break;
        case MomentumBasis::V3_PtotSlopesAsinhTimeAsinh:
          forwardTransformSampleV2(x, y, z, t, px, py, pz, x0, y0, t0, tScale, p0, VDr, VDz0,
                                   xTrans, yTrans, tTrans, m0, m1, m2,
                                   /*asinhSlopes=*/true, /*asinhTime=*/true, pzStats, posBasis);
          break;
        case MomentumBasis::V1_CylindricalTransformed:
        default:
          forwardTransformSampleV1(x, y, z, t, px, py, pz, x0, y0, t0, tScale, p0, VDr, VDz0,
                                   xTrans, yTrans, tTrans, m0, m1, m2, posBasis);
          break;
      }
    }

    inline void invertGeneratedSample(
      const double xTrans, const double yTrans, const double tTrans,
      const double m0, const double m1, const double m2,
      const double x0, const double y0, const double t0,
      const double tScale, const double p0, const double VDr, const double VDz0,
      double& x, double& y, double& z, double& t,
      double& px, double& py, double& pz,
      const MomentumBasis basis = MomentumBasis::V1_CylindricalTransformed,
      const PositionBasis posBasis = PositionBasis::V1_Atanh)
    {
      switch (basis) {
        case MomentumBasis::V2_PtotSlopes:
          invertGeneratedSampleV2(xTrans, yTrans, tTrans, m0, m1, m2, x0, y0, t0, tScale, p0, VDr, VDz0,
                                  x, y, z, t, px, py, pz, /*asinhSlopes=*/false, /*asinhTime=*/false,
                                  posBasis);
          break;
        case MomentumBasis::V2_PtotSlopesAsinh:
          invertGeneratedSampleV2(xTrans, yTrans, tTrans, m0, m1, m2, x0, y0, t0, tScale, p0, VDr, VDz0,
                                  x, y, z, t, px, py, pz, /*asinhSlopes=*/true, /*asinhTime=*/false,
                                  posBasis);
          break;
        case MomentumBasis::V3_PtotSlopesAsinhTimeAsinh:
          invertGeneratedSampleV2(xTrans, yTrans, tTrans, m0, m1, m2, x0, y0, t0, tScale, p0, VDr, VDz0,
                                  x, y, z, t, px, py, pz, /*asinhSlopes=*/true, /*asinhTime=*/true,
                                  posBasis);
          break;
        case MomentumBasis::V1_CylindricalTransformed:
        default:
          invertGeneratedSampleV1(xTrans, yTrans, tTrans, m0, m1, m2, x0, y0, t0, tScale, p0, VDr, VDz0,
                                  x, y, z, t, px, py, pz, posBasis);
          break;
      }
    }

  } // namespace VDResampler
} // namespace mu2e
