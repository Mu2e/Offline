// predicate to extrapolate to passive calorimeter front panel materials (foam, carbon)
// Sophie Middleton (2026)
#ifndef Mu2eKinKal_ExtrapolateCaloMaterial_hh
#define Mu2eKinKal_ExtrapolateCaloMaterial_hh
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/KinKalGeom/inc/Calo.hh"
#include "cetlib_except/exception.h"
#include <limits>
#include "TH2F.h"
namespace mu2e {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KinKal::Intersection;
  using KKGeom::Calo;
  using KinKal::Annulus;
  using AnnPtr = std::shared_ptr<KinKal::Annulus>;

  class ExtrapolateCaloMaterial {
    public:
      // Surface types for calo passive materials
      enum SurfaceType { FrontPanelFoam = 0, FrontPanelCarbon = 1, Unknown = 2 };

      ExtrapolateCaloMaterial() : maxDt_(-1.0), dptol_(1e10), intertol_(1e10),
      fp_z_(std::numeric_limits<double>::max()), disk_(0),
      inter_(), surftype_(Unknown), intersection_found_(false), fpann_(nullptr), debug_(0) {}

      // Constructor that initializes from calorimeter geometry
      // depth: which disk (0 or 1)
      ExtrapolateCaloMaterial(double maxdt, double dptol, double intertol, Calo const& calo,
                             int depth, int debug=0);

      // interface for extrapolation
      double maxDt() const { return maxDt_; }
      double dpTolerance() const { return dptol_; }
      double interTolerance() const { return intertol_; }
      auto const& frontPanelAnnulus() const { return fpann_; }
      auto const& intersection() const { return inter_; }
      auto surfaceType() const { return surftype_; }
      int debug() const { return debug_; }

      // extrapolation predicate: returns true if track should continue through material region
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const;
      // reset between tracks
      void reset() const { inter_ = Intersection(); surftype_ = Unknown; intersection_found_ = false; }

      // Validation plot methods - histograms created externally via TFileService
      void setValidationHistograms(TH2F* h_eff, TH2F* h_pos, TH1F* h_deg = nullptr, TH1F* h_scat100 = nullptr, TH1F* h_scat80 = nullptr) {
        h_intersection_efficiency_ = h_eff;
        h_frontpanel_hits_ = h_pos;
        h_momentum_degradation_ = h_deg;
        h_scatter_100mev_ = h_scat100;
        h_scatter_80mev_ = h_scat80;
      }
      void fillValidationPlots(double momentum, int disk, double pos_x, double pos_y, double dmom = 0.0, double scatter_angle = 0.0) const;

      // thickness of passive materials (should come from geometry TODO)
      static constexpr double foamThickness_ = 21.75;  // mm, front panel foam
      static constexpr double carbonThickness_ = 3.0;  // mm, front panel carbon

    private:
      double maxDt_; // maximum extrapolation time
      double dptol_; // fractional momentum tolerance
      double intertol_; // intersection tolerance (mm)
      double fp_z_; // z position of front panel
      int disk_; // disk ID (0 or 1)
      mutable Intersection inter_; // cache of most recent intersection
      mutable SurfaceType surftype_; // type of surface that was intersected
      mutable bool intersection_found_; // flag to prevent finding same intersection twice
      AnnPtr fpann_; // annulus surface for front panel
      int debug_; // debug level

      // Validation histograms (managed externally by TFileService)
      mutable TH2F* h_intersection_efficiency_ = nullptr; // Plot 1: efficiency vs momentum and disk
      mutable TH2F* h_frontpanel_hits_ = nullptr;         // Plot 2: hit distribution on front panel (phi vs r)
      mutable TH1F* h_momentum_degradation_ = nullptr;    // Plot 3: momentum degradation distribution
      mutable TH1F* h_scatter_100mev_ = nullptr;          // Plot 4: scattering angle at 100 MeV
      mutable TH1F* h_scatter_80mev_ = nullptr;           // Plot 4: scattering angle at 80 MeV
  };

  template <class KTRAJ> bool ExtrapolateCaloMaterial::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& fittraj, TimeDir tdir) const {
    // Calorimeter has only ONE surface (unlike multi-foil predicates like ST).
    // Once we've found the intersection, we must stop - don't search again!
    if(debug_ == -300) {
      std::cout << "    [needsExtrapolation] CHECK: intersection_found_ = " << intersection_found_ << std::endl;
    }
    if(intersection_found_) {
      if(debug_ == -300) std::cout << "    [needsExtrapolation] Already found intersection, returning false" << std::endl;
      return false;  // Already found the one intersection - stop
    }

    // Check if we need to continue extrapolating to find calorimeter front panel intersection
    // Use FULL trajectory range (particles already past calorimeter position after tracker fit)
    auto const& trange = fittraj.range();
    auto stime = trange.begin();
    auto etime = trange.end();

    auto spos = fittraj.position3(stime);
    auto epos = fittraj.position3(etime);

    // Get velocity at start
    auto vel_unit = fittraj.direction(stime);
    double zvel = vel_unit.Z();
    double zvel_signed = zvel * (tdir == TimeDir::forwards ? 1.0 : -1.0);

    if(debug_ == -300) {
      std::cout << "\n[needsExtrapolation] Checking front panel crossing (full trajectory):" << std::endl;
      std::cout << "  Trajectory Z range: [" << spos.Z() << ", " << epos.Z() << "] mm" << std::endl;
      std::cout << "  Z velocity: " << zvel << " (signed: " << zvel_signed << ")" << std::endl;
      std::cout << "  Front panel Z: " << fp_z_ << " mm" << std::endl;
      std::cout << "  Time range: [" << trange.begin() << ", " << trange.end() << "] ns" << std::endl;
    }

    // Quick Z range check - is the front panel even in our trajectory volume?
    double zmin = std::min(spos.Z(), epos.Z());
    double zmax = std::max(spos.Z(), epos.Z());

    if(debug_ == -300) {
      std::cout << "    [needsExtrapolation] Z-range check: zmin=" << zmin << " zmax=" << zmax << " fp_z_=" << fp_z_ << std::endl;
      std::cout << "    [needsExtrapolation] Panel in range? " << (fp_z_ >= zmin && fp_z_ <= zmax) << std::endl;
    }

    if(fp_z_ < zmin || fp_z_ > zmax) {
      // Front panel not in trajectory range yet
      if(debug_ == -300) {
        std::cout << "  Front panel Z " << fp_z_ << " outside trajectory Z range [" << zmin << ", " << zmax << "]" << std::endl;
      }
      return true;  // keep going
    }

    // Front panel is within trajectory Z range - test for intersection with annulus
    if(debug_ == -300) {
      std::cout << "  Testing intersection with annulus (full trajectory):" << std::endl;
    }
    Intersection newinter = KinKal::intersect(fittraj, *fpann_, trange, intertol_);

    if(debug_ > 0 || debug_ == -300) {
      std::cout << "[ExtrapolateCaloMaterial::needsExtrapolation] Intersection test result:" << std::endl;
      std::cout << "  good() = " << newinter.good() << std::endl;
      if(newinter.good()) {
        std::cout << "  *** INTERSECTION FOUND ***" << std::endl;
        std::cout << "  Position: (" << newinter.pos_.X() << ", " \
                  << newinter.pos_.Y() << ", " << newinter.pos_.Z() << ")" << std::endl;
      } else {
        std::cout << "  No intersection detected" << std::endl;
      }
    }

    if(newinter.good()){
      inter_ = newinter;
      intersection_found_ = true;  // Mark that we've found THE intersection

      // Fill validation plots
      double momentum = fittraj.momentum(newinter.time_);
      auto pos = newinter.pos_;
      fillValidationPlots(momentum, disk_, pos.X(), pos.Y());

      if(debug_ == -300) {
        std::cout << "  Good intersection with front panel FOUND" << std::endl;
      }
      return false;  // stop, we found material crossing
    } else {
      if(debug_ == -300) {
        std::cout << "  No intersection found with annulus in Z-range" << std::endl;
      }
      return false;  // stop - panel in range but no intersection means we're past it
    }
  }
}
#endif
