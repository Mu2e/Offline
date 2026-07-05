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

    // safety constants for numerical stability in the forward and inverse transforms
    constexpr double kRadiusSafetyEpsilon = 1e-6;
    constexpr double kMinSafeTime = 0.1;
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
    // In the all-at-once vector t,x,y occupy slots 0,1,2 and momentum occupies 3,4,5.
    // ------------------------------------------------------------------------
    enum class MomentumBasis {
      V1_CylindricalTransformed = 0,
      V2_PtotSlopes             = 1,
      V2_PtotSlopesAsinh        = 2
    };

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
    // Opaque basis-tag packing. The SBDM stores a single int32 it never
    // interprets; VDResampler packs (ModelLayout, MomentumBasis) into it and
    // unpacks on load. Encoding: layout*100 + basis. A pre-tag (v<=6) model loads
    // as tag 0 = (AllAtOnce6D, V1_CylindricalTransformed), preserving old behavior.
    // ------------------------------------------------------------------------
    constexpr int kBasisTagLayoutStride = 100;

    inline int packBasisTag(ModelLayout layout, MomentumBasis basis) {
      return static_cast<int>(layout) * kBasisTagLayoutStride + static_cast<int>(basis);
    }
    inline ModelLayout unpackModelLayout(int tag) {
      return static_cast<ModelLayout>(tag / kBasisTagLayoutStride);
    }
    inline MomentumBasis unpackMomentumBasis(int tag) {
      return static_cast<MomentumBasis>(tag % kBasisTagLayoutStride);
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
    inline const char* momentumBasisName(MomentumBasis b) {
      switch (b) {
        case MomentumBasis::V1_CylindricalTransformed: return "V1_CylindricalTransformed";
        case MomentumBasis::V2_PtotSlopes:             return "V2_PtotSlopes";
        case MomentumBasis::V2_PtotSlopesAsinh:        return "V2_PtotSlopesAsinh";
      }
      return "unknown";
    }
    // Decode an opaque basisTag (layout*100 + basis) into
    // "layout=<...>, basis=<...> (tag <n>)" for logs and error messages.
    inline std::string basisTagToString(int tag) {
      std::ostringstream os;
      os << "layout=" << modelLayoutName(unpackModelLayout(tag))
         << ", basis=" << momentumBasisName(unpackMomentumBasis(tag))
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

    // Forward position transform: (dx,dy,r) -> (xTrans,yTrans) via u=atanh(r/VDr)
    // mapped back through the angle. Shared by all bases.
    inline void forwardPosition(double dx, double dy, double r, double VDr,
                                double& xTrans, double& yTrans)
    {
      double rho = r / VDr;
      rho = std::min(rho, 1.0 - kRadiusSafetyEpsilon); // avoid rho >= 1
      const double u = 0.5 * std::log((1.0 + rho) / (1.0 - rho)); // atanh(r/VDr)
      const double theta = std::atan2(dy, dx);
      xTrans = u * std::cos(theta);
      yTrans = u * std::sin(theta);
    }

    // Inverse position transform: (xTrans,yTrans) -> detector-space (x,y) and the
    // local frame (dx,dy,r) needed to rotate momentum back. Shared by all bases.
    inline void invertPosition(double xTrans, double yTrans, double x0, double y0,
                               double VDr,
                               double& x, double& y, double& dx, double& dy, double& r)
    {
      const double u = std::sqrt(xTrans * xTrans + yTrans * yTrans);
      const double theta = std::atan2(yTrans, xTrans);
      const double rho = std::tanh(u);
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
      double& prTrans, double& pphiTrans, double& pzTrans)
    {
      double dx, dy, r;
      extrapolateAndCenter(x, y, z, px, py, pz, x0, y0, VDz0, dx, dy, r);
      forwardPosition(dx, dy, r, VDr, xTrans, yTrans);

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
      double& px, double& py, double& pz)
    {
      double dx, dy, r;
      invertPosition(xTrans, yTrans, x0, y0, VDr, x, y, dx, dy, r);
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
      const bool asinhSlopes, PzFallbackStats* pzStats = nullptr)
    {
      double dx, dy, r;
      extrapolateAndCenter(x, y, z, px, py, pz, x0, y0, VDz0, dx, dy, r);
      forwardPosition(dx, dy, r, VDr, xTrans, yTrans);

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

      tTrans = forwardTime(t, t0, tScale);
    }

    inline void invertGeneratedSampleV2(
      const double xTrans, const double yTrans, const double tTrans,
      const double pTotTrans, const double urTrans, const double uphiTrans,
      const double x0, const double y0, const double t0,
      const double tScale, const double p0, const double VDr, const double VDz0,
      double& x, double& y, double& z, double& t,
      double& px, double& py, double& pz,
      const bool asinhSlopes)
    {
      double dx, dy, r;
      invertPosition(xTrans, yTrans, x0, y0, VDr, x, y, dx, dy, r);
      z = VDz0;
      t = invertTime(tTrans, t0, tScale);

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
      PzFallbackStats* pzStats = nullptr)
    {
      switch (basis) {
        case MomentumBasis::V2_PtotSlopes:
          forwardTransformSampleV2(x, y, z, t, px, py, pz, x0, y0, t0, tScale, p0, VDr, VDz0,
                                   xTrans, yTrans, tTrans, m0, m1, m2, /*asinhSlopes=*/false, pzStats);
          break;
        case MomentumBasis::V2_PtotSlopesAsinh:
          forwardTransformSampleV2(x, y, z, t, px, py, pz, x0, y0, t0, tScale, p0, VDr, VDz0,
                                   xTrans, yTrans, tTrans, m0, m1, m2, /*asinhSlopes=*/true, pzStats);
          break;
        case MomentumBasis::V1_CylindricalTransformed:
        default:
          forwardTransformSampleV1(x, y, z, t, px, py, pz, x0, y0, t0, tScale, p0, VDr, VDz0,
                                   xTrans, yTrans, tTrans, m0, m1, m2);
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
      const MomentumBasis basis = MomentumBasis::V1_CylindricalTransformed)
    {
      switch (basis) {
        case MomentumBasis::V2_PtotSlopes:
          invertGeneratedSampleV2(xTrans, yTrans, tTrans, m0, m1, m2, x0, y0, t0, tScale, p0, VDr, VDz0,
                                  x, y, z, t, px, py, pz, /*asinhSlopes=*/false);
          break;
        case MomentumBasis::V2_PtotSlopesAsinh:
          invertGeneratedSampleV2(xTrans, yTrans, tTrans, m0, m1, m2, x0, y0, t0, tScale, p0, VDr, VDz0,
                                  x, y, z, t, px, py, pz, /*asinhSlopes=*/true);
          break;
        case MomentumBasis::V1_CylindricalTransformed:
        default:
          invertGeneratedSampleV1(xTrans, yTrans, tTrans, m0, m1, m2, x0, y0, t0, tScale, p0, VDr, VDz0,
                                  x, y, z, t, px, py, pz);
          break;
      }
    }

  } // namespace VDResampler
} // namespace mu2e
