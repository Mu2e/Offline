#pragma once

#include <algorithm>
#include <cmath>

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

    // Transform detector-space quantities to model-space variables used for training.
    inline void forwardTransformSample(
      const double x,
      const double y,
      const double z,
      const double t,
      const double px,
      const double py,
      const double pz,
      const double x0,
      const double y0,
      const double t0,
      const double tScale,
      const double p0,
      const double VDr,
      const double VDz0,
      double& xTrans,
      double& yTrans,
      double& tTrans,
      double& prTrans,
      double& pphiTrans,
      double& pzTrans
    ) {
      // as z maybe slightly different from the nominal VDz0 due to the step size, we will extrapolate the (x, y) coordinates
      // to the nominal VDz0 for all hits to compute the training parameters to be fed into the SBDM.
      const double extrapolationFactor = (VDz0 - z) / pz;
      const double xExtrapolated = x + extrapolationFactor * px;
      const double yExtrapolated = y + extrapolationFactor * py;

      // Now we have the extrapolated (x, y) at the nominal VDz0, we can compute the training parameters for the SBDM.
      // We convert (t, x_extrapolated, y_extrapolated, px, py, pz) to (t', x_extrapolated, y_extrapolated, p_r', p_phi', p_z')

      // Shift to detector-centered coordinates and map to transformed position variables.
      const double dx = xExtrapolated - x0;
      const double dy = yExtrapolated - y0;
      // polar coordinates
      const double r = std::sqrt(dx * dx + dy * dy);
      double rho = r / VDr;
      // numerical safety (avoid rho >= 1)
      rho = std::min(rho, 1.0 - kRadiusSafetyEpsilon);
      // boundary-removing transform
      // u = atanh(r/R)
      const double u = 0.5 * std::log((1.0 + rho) / (1.0 - rho));
      // angle
      const double theta = std::atan2(dy, dx);
      // map back to Cartesian-like coordinates
      xTrans = u * std::cos(theta);
      yTrans = u * std::sin(theta);

      // compute momentum components in the local polar coordinate system (r, phi, z)
      double pr = 0.0;
      double pphi = 0.0;
      if (r > kRadiusSafetyEpsilon) { // avoid division by zero, if r is very small, we can approximate pr ~ px and pphi ~ py
        const double rx = dx / r; // unit vector in the radial direction
        const double ry = dy / r;
        const double phix = -ry; // unit vector in the angular direction
        const double phiy = rx;
        pr = px * rx + py * ry; // radial momentum component
        pphi = px * phix + py * phiy; // angular momentum component
      } else {
        pr = px;
        pphi = py;
      }
      // momentum scaling
      prTrans = std::asinh(pr / p0); // tunable scale where I want best resolution
      pphiTrans = std::asinh(pphi / p0);
      pzTrans = std::asinh(pz / p0);

      // time transform
      const double tSafe = (t > kMinSafeTime) ? t : kMinSafeTime; // avoid log(0)
      tTrans = std::log(tSafe / t0) / tScale;
    }

    // Invert transformed sample coordinates/momenta back to detector-space quantities.
    inline void invertGeneratedSample(
      const double xTrans,
      const double yTrans,
      const double tTrans,
      const double prTrans,
      const double pphiTrans,
      const double pzTrans,
      const double x0,
      const double y0,
      const double t0,
      const double tScale,
      const double p0,
      const double VDr,
      const double VDz0,
      double& x,
      double& y,
      double& z,
      double& t,
      double& px,
      double& py,
      double& pz
    ) {
      const double u = std::sqrt(xTrans * xTrans + yTrans * yTrans);
      const double theta = std::atan2(yTrans, xTrans);
      const double rho = std::tanh(u);
      const double r = rho * VDr;
      const double dx = r * std::cos(theta);
      const double dy = r * std::sin(theta);
      x = dx + x0;
      y = dy + y0;
      z = VDz0;

      t = t0 * std::exp(tTrans * tScale);

      const double pr = p0 * std::sinh(prTrans);
      const double pphi = p0 * std::sinh(pphiTrans);
      pz = p0 * std::sinh(pzTrans);

      if (r > kRadiusSafetyEpsilon) {
        const double rx = dx / r;
        const double ry = dy / r;
        const double phix = -ry;
        const double phiy = rx;
        px = pr * rx + pphi * phix;
        py = pr * ry + pphi * phiy;
      } else {
        px = pr;
        py = pphi;
      }
    }

  } // namespace VDResampler
} // namespace mu2e
