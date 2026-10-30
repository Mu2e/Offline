// Targets and predicate for the KKExtrap sweep: all requested surfaces are searched together from the last processed crossing.
#ifndef Mu2eKinKal_ExtrapolateSweep_hh
#define Mu2eKinKal_ExtrapolateSweep_hh
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "KinKal/Geometry/Annulus.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Plane.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>

namespace mu2e {
  // one surface the sweep can cross
  struct SweepTarget {
    enum Effect { record=0, cylshell, ipashell, stfoil, planeshell, crvsector }; // what a crossing adds to the track
    enum Reach { radial=0, plane, zplane, zrange, inside }; // when a driving target can still be reached
    static constexpr double inf = std::numeric_limits<double>::max();
    Effect effect_ = record;
    Reach reach_ = radial;
    SurfaceId sid_;
    KinKal::Surface const* surf_ = nullptr; // surface to intersect
    KinKal::Plane const* pl_ = nullptr; // the same surface as a plane, for plane reach
    std::shared_ptr<KinKal::Cylinder> cyl_; // typed surfaces for the material xings
    std::shared_ptr<KinKal::Annulus> ann_;
    std::shared_ptr<KinKal::Plane> plane_;
    std::shared_ptr<KinKal::Rectangle> rect_;
    MatEnv::DetMaterial const* mat_ = nullptr;
    double thick_ = 0.0; // material thickness (mm)
    unsigned maxcross_ = 1; // crossings kept per sweep, 0 for no limit
    int group_ = -1; // targets in a group keep only the first crossing of any of them
    bool drives_ = true; // keep extending while this target is reachable
    bool terminal_ = false; // a crossing ends the sweep
    bool inbounds_ = true; // require the crossing inside the surface bounds
    bool active_ = true; // initial state; the upstream reflection logic switches targets on and off
    double mincos_ = 0.0; // minimum |cos| between the surface normal and the direction at the crossing
    double rmax_ = inf; // unreachable beyond this transverse radius when moving outward (mm)
    double zval_ = 0.0; // plane z for zplane reach (mm)
    double rreg_ = inf, zreglo_ = -inf, zreghi_ = inf; // region for inside and zrange reach (mm)
  };
  using SweepTargetCollection = std::vector<SweepTarget>;

  class ExtrapolateSweep {
    public:
      struct Crossing { size_t itarget_ = 0; KinKal::Intersection inter_; };
      ExtrapolateSweep(SweepTargetCollection const& targets, KinKal::TimeDir tdir, double maxdt, double maxdtstep, double dptol,
          double intertol, double minv, int debug) : targets_(targets), tdir_(tdir), maxDt_(maxdt), maxDtStep_(maxdtstep),
      dptol_(dptol), intertol_(intertol), minv_(minv), debug_(debug), active_(targets.size()), count_(targets.size(),0),
      tlast_(targets.size(),0.0), coslast_(targets.size(),1.0), episodes_(targets.size()), planes_(targets.size(),nullptr), cylinders_(targets.size(),nullptr) {
        for(size_t itg = 0; itg < targets_.size(); ++itg){
          active_[itg] = targets_[itg].active_;
          planes_[itg] = dynamic_cast<KinKal::Plane const*>(targets_[itg].surf_);
          cylinders_[itg] = dynamic_cast<KinKal::Cylinder const*>(targets_[itg].surf_);
        }
      }
      // interface used by KinKal::Track::extrapolate
      double maxDt() const { return maxDt_; }
      double maxDtStep() const { return maxDtStep_; }
      double dpTolerance() const { return dptol_; }
      int debug() const { return debug_; }
      // false (stop extending) once a crossing is cached, nothing driving is reachable, a watched reflection happened or the search failed
      template <class KTRAJ> bool needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, KinKal::TimeDir tdir) const;
      bool failed() const { return !error_.empty(); }
      std::string const& error() const { return error_; } // why the search failed; KinKal's walk never sees the exception
      // sweep control
      void start(double tcursor) { searched_ = tcursor; found_ = reflected_ = false; }
      void limit(double tlimit, bool recordonly) { limited_ = true; tlimit_ = tlimit; recordonly_ = recordonly; }
      void unlimit() { limited_ = recordonly_ = false; }
      void watchReflection(double vbsign) { vbsign_ = vbsign; reflected_ = false; }
      void setActive(size_t itg, bool active) { active_[itg] = active; }
      bool hasCrossing() const { return found_; }
      bool reflected() const { return reflected_; }
      Crossing const& crossing() const { return crossing_; }
      SweepTargetCollection const& targets() const { return targets_; }
      // mark the cached crossing processed and move the search start to it
      void accept();
    private:
      template <class KTRAJ> bool extend(KinKal::ParticleTrajectory<KTRAJ> const& ptraj) const;
      bool usable(size_t itg) const;
      bool withinReach(size_t itg, KinKal::VEC3 const& pos, double reach) const;
      template <class KTRAJ> bool search(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, double tend) const;
      template <class KTRAJ> bool firstCrossing(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, KTRAJ const& piece, size_t itg, KinKal::TimeRange win, KinKal::Intersection& inter) const;
      template <class KTRAJ, class VISIT> bool visitPlaneExtrema(KTRAJ const& piece, KinKal::Plane const& plane, KinKal::TimeRange const& range, VISIT&& visit) const;
      template <class KTRAJ> KinKal::Intersection planeCrossing(KTRAJ const& piece, KinKal::Plane const& plane, KinKal::TimeRange const& query) const;
      template <class KTRAJ> void ownBoundary(KTRAJ const& piece, SweepTarget const& tg, KinKal::Intersection& inter) const;
      template <class KTRAJ> bool seam(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, size_t itg, KinKal::Intersection const& inter) const;
      template <class KTRAJ> bool departed(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, size_t itg, size_t previous, double accepted, double proposed) const;
      template <class KTRAJ> bool reachable(KinKal::ParticleTrajectory<KTRAJ> const& ptraj) const;
      SweepTargetCollection const& targets_;
      KinKal::TimeDir tdir_;
      double maxDt_, maxDtStep_, dptol_, intertol_, minv_;
      int debug_;
      std::vector<bool> active_; // per target
      std::vector<unsigned> count_; // crossings processed per target
      std::vector<double> tlast_; // time of the last processed crossing per target
      std::vector<double> coslast_; // |cos| between the surface normal and the direction at that crossing
      std::vector<int> doneGroups_; // groups whose first crossing is processed
      struct SeamEpisode { size_t itarget_; double time_; };
      std::vector<std::vector<SeamEpisode>> episodes_; // accepted seam episodes per eligible target, including self
      std::vector<KinKal::Plane const*> planes_; // per target, the surface as a plane or cylinder when it is one
      std::vector<KinKal::Cylinder const*> cylinders_;
      bool limited_ = false, recordonly_ = false;
      double tlimit_ = 0.0;
      double vbsign_ = 0.0; // sign of the B-parallel velocity to watch for a reflection; 0 to not watch
      mutable double searched_ = 0.0; // the trajectory is searched up to this time
      mutable bool found_ = false, reflected_ = false;
      mutable Crossing crossing_;
      mutable std::string error_;
      static constexpr unsigned maxRootTrials_ = 256; // per target, piece and window
  };

  inline bool ExtrapolateSweep::usable(size_t itg) const {
    auto const& tg = targets_[itg];
    if(!active_[itg]) return false;
    if(recordonly_ && tg.effect_ != SweepTarget::record) return false;
    if(tg.maxcross_ > 0 && count_[itg] >= tg.maxcross_) return false;
    if(tg.group_ >= 0 && std::find(doneGroups_.begin(),doneGroups_.end(),tg.group_) != doneGroups_.end()) return false;
    return true;
  }

  // A trajectory moves at most speed x time, so a surface farther than that from the window start cannot be crossed in it.
  inline bool ExtrapolateSweep::withinReach(size_t itg, KinKal::VEC3 const& pos, double reach) const {
    if(planes_[itg]) return std::fabs(planes_[itg]->distance(pos)) <= reach;
    if(cylinders_[itg]){
      auto const& cyl = *cylinders_[itg];
      KinKal::VEC3 rel = pos - cyl.center();
      double along = rel.Dot(cyl.axis());
      double rho = (rel - along*cyl.axis()).R();
      if(std::fabs(rho - cyl.radius()) > reach) return false;
      return !targets_[itg].inbounds_ || std::fabs(along) <= cyl.halfLength() + reach;
    }
    return true;
  }

  inline void ExtrapolateSweep::accept() {
    if(!found_) return;
    auto itg = crossing_.itarget_;
    auto const& tg = targets_[itg];
    double time = crossing_.inter_.time_;
    ++count_[itg];
    tlast_[itg] = time;
    coslast_[itg] = std::fabs(crossing_.inter_.norm_.Dot(crossing_.inter_.pdir_));
    if(tg.group_ >= 0) doneGroups_.push_back(tg.group_);
    if((tg.effect_ == SweepTarget::planeshell || tg.effect_ == SweepTarget::crvsector || tg.effect_ == SweepTarget::stfoil) && tg.pl_){
      for(size_t other = 0; other < targets_.size(); ++other){
        auto const& neighbor = targets_[other];
        if(neighbor.effect_ != tg.effect_ || !neighbor.pl_ || neighbor.mat_ != tg.mat_ || neighbor.thick_ != tg.thick_) continue;
        if(std::fabs(neighbor.pl_->normal().Dot(tg.pl_->normal())) < 1.0 - 1.0e-12) continue;
        if(std::fabs((neighbor.pl_->center()-tg.pl_->center()).Dot(tg.pl_->normal())) > intertol_) continue;
        if(!tg.surf_->inBounds(crossing_.inter_.pos_,intertol_) || !neighbor.surf_->inBounds(crossing_.inter_.pos_,intertol_)) continue;
        auto& episodes = episodes_[other];
        auto previous = std::find_if(episodes.begin(),episodes.end(),[itg](SeamEpisode const& episode){ return episode.itarget_ == itg; });
        // A later acceptance of the same source supersedes its earlier episode without erasing other planes.
        if(previous == episodes.end()) episodes.push_back({itg,time});
        else previous->time_ = time;
      }
    }
    searched_ = time;
    found_ = false;
  }

  template <class KTRAJ> bool ExtrapolateSweep::needsExtrapolation(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, KinKal::TimeDir) const {
    if(failed()) return false;
    try {
      return extend(ptraj);
    } catch(std::bad_alloc const&) {
      throw;
    } catch(std::exception const& err) {
      // a problem with this track must not reach KinKal's field walk, which would mark the whole fit unusable
      error_ = err.what();
      found_ = false;
      return false;
    }
  }

  template <class KTRAJ> bool ExtrapolateSweep::extend(KinKal::ParticleTrajectory<KTRAJ> const& ptraj) const {
    found_ = false;
    double sgn = tdir_ == KinKal::TimeDir::forwards ? 1.0 : -1.0;
    double tend = tdir_ == KinKal::TimeDir::forwards ? ptraj.range().end() : ptraj.range().begin();
    if(limited_ && sgn*(tend - tlimit_) > 0.0) tend = tlimit_;
    if(search(ptraj,tend)) return false;
    if(limited_) return false;
    if(vbsign_ != 0.0){
      auto const& ktraj = tdir_ == KinKal::TimeDir::forwards ? ptraj.back() : ptraj.front();
      double vb = ktraj.velocity(tend).Dot(ktraj.bnom().Unit())*sgn;
      if(vb*vbsign_ < 0.0){
        if(debug_ > 1) std::cout << "Sweep: reflection at time " << tend << " z " << ktraj.position3(tend).Z() << std::endl;
        reflected_ = true;
        return false;
      }
    }
    return reachable(ptraj);
  }

  template <class KTRAJ, class VISIT> bool ExtrapolateSweep::visitPlaneExtrema(KTRAJ const& piece,
      KinKal::Plane const& plane, KinKal::TimeRange const& range, VISIT&& visit) const {
    double omega = 0.0;
    if constexpr(requires { piece.Omega(); }) omega = std::fabs(piece.Omega());
    else if constexpr(requires { piece.omega(); }) omega = std::fabs(piece.omega());
    if(!std::isfinite(omega)) throw std::runtime_error("KKExtrap: nonfinite plane angular frequency");
    if(omega == 0.0) return false;
    double sign = tdir_ == KinKal::TimeDir::forwards ? 1.0 : -1.0;
    double pi = std::acos(-1.0), base = sign > 0.0 ? range.begin() : range.end();
    double stop = sign > 0.0 ? range.end() : range.begin();
    double half = base+pi/omega, quarter = base+0.5*pi/omega;
    if(!std::isfinite(half) || !std::isfinite(quarter))
      throw std::runtime_error("KKExtrap: nonfinite plane velocity sample time");
    double v0 = plane.normal().Dot(piece.velocity(base));
    double vhalf = plane.normal().Dot(piece.velocity(half));
    double constant = 0.5*(v0+vhalf), cosine = 0.5*(v0-vhalf);
    double sine = plane.normal().Dot(piece.velocity(quarter))-constant;
    double amplitude = std::hypot(cosine,sine);
    if(!std::isfinite(constant) || !std::isfinite(amplitude))
      throw std::runtime_error("KKExtrap: nonfinite plane velocity coefficients");
    if(amplitude == 0.0 || std::fabs(constant) > amplitude) return false;
    double phase = std::atan2(sine,cosine), angle = std::acos(std::clamp(-constant/amplitude,-1.0,1.0));
    std::array<double,2> next;
    unsigned index = 0;
    double period = 2.0*pi/omega;
    for(double solution : {phase-angle,phase+angle}){
      double normalized = std::fmod(sign*solution,2.0*pi);
      if(normalized <= 0.0) normalized += 2.0*pi;
      next[index++] = base+sign*normalized/omega;
    }
    unsigned trials = 0;
    while(true){
      index = sign*next[0] < sign*next[1] ? 0 : 1;
      double time = next[index];
      if(!std::isfinite(time) || sign*(stop-time) <= 0.0) return false;
      if(sign*(time-base) <= 0.0) throw std::runtime_error("KKExtrap: plane extremum did not advance");
      if(++trials > maxRootTrials_) throw std::runtime_error("KKExtrap: plane search exceeded extrema budget");
      if(visit(time)) return true;
      next[index] = time+sign*period;
      if(!std::isfinite(next[index]) || sign*(next[index]-time) <= 0.0)
        throw std::runtime_error("KKExtrap: plane extremum did not advance");
    }
  }

  template <class KTRAJ> KinKal::Intersection ExtrapolateSweep::planeCrossing(KTRAJ const& piece,
      KinKal::Plane const& plane, KinKal::TimeRange const& query) const {
    auto distance = [&](double time){
      double value = plane.distance(piece.position3(time));
      if(!std::isfinite(value)) throw std::runtime_error("KKExtrap: nonfinite plane distance");
      return value;
    };
    auto root = [&](double time){
      KinKal::Intersection inter;
      inter.time_ = time;
      inter.pos_ = piece.position3(time);
      inter.pdir_ = piece.direction(time);
      inter.norm_ = plane.normal();
      inter.onsurface_ = plane.onSurface(inter.pos_,intertol_);
      inter.inbounds_ = plane.inBounds(inter.pos_,intertol_);
      inter.inrange_ = query.inRange(time) && piece.range().inRange(time);
      return inter;
    };
    double cursor = tdir_ == KinKal::TimeDir::forwards ? query.begin() : query.end();
    auto endpointContact = [&](double time){
      auto inter = root(time);
      double normalSpeed = std::fabs(plane.normal().Dot(piece.velocity(time)));
      return inter.onsurface_ && inter.inrange_ && normalSpeed <= 1.0e-12*std::max(piece.speed(time),1.0);
    };
    if(endpointContact(cursor)) return root(cursor);
    KinKal::Intersection result;
    // Plane distance is monotonic between consecutive extrema of the supported line and helix pieces.
    auto bracket = [&](double stop){
      double a = cursor, b = stop;
      cursor = stop;
      double fa = distance(a), fb = distance(b);
      if(fa == 0.0 && query.inRange(a)){ result = root(a); return true; }
      if(fb == 0.0 && query.inRange(b)){ result = root(b); return true; }
      if(fa == 0.0 || fb == 0.0 || std::signbit(fa) == std::signbit(fb)) return false;
      for(unsigned iteration = 0; iteration < maxRootTrials_; ++iteration){
        double mid = a+0.5*(b-a);
        if(mid == a || mid == b){
          for(double time : {a,b}){
            auto inter = root(time);
            if(inter.onsurface_ && inter.inrange_){ result = inter; return true; }
          }
          throw std::runtime_error("KKExtrap: plane bracket has no representable root");
        }
        double fm = distance(mid), speed = piece.speed(mid);
        if(!std::isfinite(speed) || !(speed > 0.0)) throw std::runtime_error("KKExtrap: invalid plane-root speed");
        if(fm == 0.0 || (std::fabs(b-a)*speed <= 0.125*intertol_ && std::fabs(fm) < intertol_)){
          result = root(mid);
          return true;
        }
        if(std::signbit(fa) == std::signbit(fm)){ a = mid; fa = fm; }
        else { b = mid; fb = fm; }
      }
      throw std::runtime_error("KKExtrap: plane bracket exceeded bisection budget");
    };
    auto segment = [&](double time){
      if(bracket(time)) return true;
      // A tolerance-level extremum is a recordable contact; material incidence cuts still apply.
      auto contact = root(time);
      if(contact.onsurface_ && contact.inrange_){ result = contact; return true; }
      return false;
    };
    if(visitPlaneExtrema(piece,plane,query,segment)) return result;
    double stop = tdir_ == KinKal::TimeDir::forwards ? query.end() : query.begin();
    if(!bracket(stop) && endpointContact(stop)) return root(stop);
    return result;
  }

  template <class KTRAJ> bool ExtrapolateSweep::departed(KinKal::ParticleTrajectory<KTRAJ> const& ptraj,
      size_t itg, size_t previous, double accepted, double proposed) const {
    auto const& plane = *targets_[itg].pl_;
    auto const& prev = *targets_[previous].pl_;
    double sgn = tdir_ == KinKal::TimeDir::forwards ? 1.0 : -1.0;
    if(sgn*(proposed-accepted) <= 0.0) return false;
    KinKal::TimeRange interval(std::min(accepted,proposed),std::max(accepted,proposed));
    if(!interval.restrict(ptraj.range()) || interval.null()) return false;
    size_t first = ptraj.nearestIndex(interval.begin()), last = ptraj.nearestIndex(interval.end());
    double offset = 2.0*intertol_ + std::fabs((plane.center()-prev.center()).Dot(prev.normal()));
    for(size_t ip = first; ip <= last; ++ip){
      auto const& piece = ptraj.piece(ip);
      auto clipped = piece.range();
      if(!clipped.restrict(interval) || clipped.null()) continue;
      auto witness = [&](double time){
        if(!std::isfinite(time) || time < clipped.begin() || time > clipped.end() ||
            sgn*(time-accepted) <= 0.0 || sgn*(proposed-time) <= 0.0) return false;
        auto pos = ptraj.position3(time);
        return std::fabs(plane.distance(pos)) > intertol_ && std::fabs(prev.distance(pos)) > intertol_;
      };
      if(witness(clipped.begin()) || witness(clipped.end())) return true;
      if(visitPlaneExtrema(piece,prev,clipped,witness)) return true;
      // Offset-plane roots also provide directly verified departure witnesses.
      for(double side : {-1.0,1.0}){
        KinKal::Plane outside(prev.normal(),prev.uDirection(),prev.center()+side*offset*prev.normal());
        auto root = KinKal::intersect(piece,outside,clipped,0.125*intertol_,tdir_);
        if(root.onsurface_ && witness(root.time_)) return true;
      }
    }
    return false;
  }

  template <class KTRAJ> bool ExtrapolateSweep::seam(KinKal::ParticleTrajectory<KTRAJ> const& ptraj,
      size_t itg, KinKal::Intersection const& inter) const {
    // Candidate searches are speculative; only accept() changes episodes, never a future departure witness.
    for(auto const& episode : episodes_[itg])
      if(!departed(ptraj,itg,episode.itarget_,episode.time_,inter.time_)) return true;
    return false;
  }

  template <class KTRAJ> void ExtrapolateSweep::ownBoundary(KTRAJ const& piece, SweepTarget const& tg, KinKal::Intersection& inter) const {
    double time = inter.time_;
    double begin = piece.range().begin(), end = piece.range().end();
    if(time != begin && time != end) return;
    // Encode exact piece-edge roots one ULP inside their source piece so lookup and material insertion agree.
    double inside = std::nextafter(time,time == begin ? end : begin);
    if(!(inside > begin && inside < end) || !std::isfinite(inside))
      throw std::runtime_error("KKExtrap: boundary root has no representable piece interior");
    double distance = std::fabs(inside-time)*piece.speed(time);
    if(!std::isfinite(distance) || distance > intertol_)
      throw std::runtime_error("KKExtrap: boundary-root displacement exceeds intersection tolerance");
    inter.time_ = inside;
    inter.pos_ = piece.position3(inside);
    inter.pdir_ = piece.direction(inside);
    inter.norm_ = tg.surf_->normal(inter.pos_);
    inter.onsurface_ = tg.surf_->onSurface(inter.pos_,intertol_);
    inter.inbounds_ = tg.surf_->inBounds(inter.pos_,intertol_);
    inter.inrange_ = piece.range().inRange(inside);
    if(!inter.onsurface_) throw std::runtime_error("KKExtrap: boundary-root displacement leaves surface tolerance");
  }

  template <class KTRAJ> bool ExtrapolateSweep::firstCrossing(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, KTRAJ const& piece, size_t itg, KinKal::TimeRange win, KinKal::Intersection& inter) const {
    auto const& tg = targets_[itg];
    double sgn = tdir_ == KinKal::TimeDir::forwards ? 1.0 : -1.0;
    double ts = tdir_ == KinKal::TimeDir::forwards ? win.begin() : win.end();
    double stop = tdir_ == KinKal::TimeDir::forwards ? win.end() : win.begin();
    double speed = std::max(piece.speed(ts),1.0);
    double step = 4.0*intertol_/speed;
    // a track at incidence |cos| stays within the tolerance for about tolerance/|cos| of path; skip past it before finding this target again
    double leave = tg.mincos_ > 0.0 ? 2.0*intertol_/(std::max(coslast_[itg],tg.mincos_)*speed) : 0.0;
    double guard = std::max(step,leave);
    if(count_[itg] > 0) ts = sgn > 0.0 ? std::max(ts,tlast_[itg]+guard) : std::min(ts,tlast_[itg]-guard);
    auto const* plane = dynamic_cast<KinKal::Plane const*>(tg.surf_);
    unsigned nroot(0);
    while(sgn*(stop - ts) > 0.0){
      if(++nroot > maxRootTrials_) throw std::runtime_error("KKExtrap: surface-root search exceeded retry limit");
      KinKal::TimeRange query = sgn > 0.0 ? KinKal::TimeRange(ts,stop) : KinKal::TimeRange(stop,ts);
      KinKal::VEC3 pos = piece.position3(ts);
      if((nroot == 1 || !plane) && tg.surf_->onSurface(pos,intertol_) &&
          std::fabs(tg.surf_->normal(pos).Dot(piece.direction(ts))) >= minv_){
        inter = KinKal::Intersection();
        inter.time_ = ts;
        inter.pos_ = pos;
        inter.pdir_ = piece.direction(ts);
        inter.norm_ = tg.surf_->normal(pos);
        inter.onsurface_ = true;
        inter.inbounds_ = tg.surf_->inBounds(pos,intertol_);
        inter.inrange_ = ts >= piece.range().begin() && ts <= piece.range().end();
      } else if(plane) {
        inter = planeCrossing(piece,*plane,query);
      } else {
        inter = KinKal::intersect(piece,*tg.surf_,query,intertol_,tdir_);
      }
      if(debug_ > 3) std::cout << "Sweep " << tg.sid_ << " " << inter << std::endl;
      if(!inter.onsurface_) return false;
      if(!std::isfinite(inter.time_)) throw std::runtime_error("KKExtrap: nonfinite surface-root time");
      if(!inter.inrange_) return false;
      if(sgn*(inter.time_ - ts) < 0.0) throw std::runtime_error("KKExtrap: surface-root search did not advance");
      if(sgn*(stop - inter.time_) < 0.0 || inter.time_ < piece.range().begin() || inter.time_ > piece.range().end())
        throw std::runtime_error("KKExtrap: surface root outside its piece/window");
      ownBoundary(piece,tg,inter);
      if(sgn*(inter.time_ - ts) < 0.0 || sgn*(stop - inter.time_) < 0.0) return false;
      if(!inter.inrange_) throw std::runtime_error("KKExtrap: surface root has no owning piece");
      bool grazing = std::fabs(inter.norm_.Dot(inter.pdir_)) < tg.mincos_;
      if((inter.inbounds_ || !tg.inbounds_) && !grazing && !seam(ptraj,itg,inter)) return true;
      // Roots within four intersection tolerances of path length are unresolved by this search; a grazing track first leaves the tolerance.
      double next = inter.time_ + sgn*(grazing ? std::max(step,2.0*intertol_/(tg.mincos_*speed)) : step);
      if(!std::isfinite(next) || sgn*(next - ts) <= 0.0)
        throw std::runtime_error("KKExtrap: surface-root rejection did not advance");
      ts = next;
    }
    return false;
  }

  template <class KTRAJ> bool ExtrapolateSweep::search(KinKal::ParticleTrajectory<KTRAJ> const& ptraj, double tend) const {
    double sgn = tdir_ == KinKal::TimeDir::forwards ? 1.0 : -1.0;
    double t0 = std::clamp(searched_,ptraj.range().begin(),ptraj.range().end());
    tend = std::clamp(tend,ptraj.range().begin(),ptraj.range().end());
    double dtwin = maxDtStep_ > 0.0 ? maxDtStep_ : std::fabs(tend - t0);
    while(sgn*(tend - t0) > 0.0){
      double t1 = t0 + sgn*std::min(dtwin,std::fabs(tend - t0));
      if(sgn*(t1 - t0) <= 0.0) throw std::runtime_error("KKExtrap: search window did not advance");
      KinKal::TimeRange window(std::min(t0,t1),std::max(t0,t1));
      bool found(false);
      Crossing best;
      size_t first = ptraj.nearestIndex(window.begin());
      size_t last = ptraj.nearestIndex(window.end());
      for(size_t ip = first; ip <= last; ++ip){
        auto const& piece = ptraj.piece(ip);
        auto win = piece.range();
        if(!win.restrict(window) || win.null()) continue;
        KinKal::VEC3 wpos = piece.position3(tdir_ == KinKal::TimeDir::forwards ? win.begin() : win.end());
        double reach = piece.speed(win.mid())*win.range()*(1.0+1.0e-9) + intertol_;
        for(size_t itg = 0; itg < targets_.size(); ++itg){
          if(!usable(itg) || !withinReach(itg,wpos,reach)) continue;
          KinKal::Intersection inter;
          if(!firstCrossing(ptraj,piece,itg,win,inter)) continue;
          bool earlier = !found || sgn*inter.time_ < sgn*best.inter_.time_;
          bool tie = found && inter.time_ == best.inter_.time_ && targets_[itg].sid_ < targets_[best.itarget_].sid_;
          if(earlier || tie){
            found = true;
            best.itarget_ = itg;
            best.inter_ = inter;
          }
        }
      }
      if(found){
        crossing_ = best;
        found_ = true;
        searched_ = t0;
        return true;
      }
      t0 = t1;
      searched_ = t1;
    }
    return false;
  }

  template <class KTRAJ> bool ExtrapolateSweep::reachable(KinKal::ParticleTrajectory<KTRAJ> const& ptraj) const {
    double sgn = tdir_ == KinKal::TimeDir::forwards ? 1.0 : -1.0;
    auto const& ktraj = tdir_ == KinKal::TimeDir::forwards ? ptraj.back() : ptraj.front();
    double time = tdir_ == KinKal::TimeDir::forwards ? ktraj.range().end() : ktraj.range().begin();
    KinKal::VEC3 pos = ktraj.position3(time);
    KinKal::VEC3 vel = ktraj.velocity(time)*sgn; // signed by the extrapolation direction
    double rho = std::sqrt(pos.X()*pos.X() + pos.Y()*pos.Y());
    double vrho = rho > 0.0 ? (pos.X()*vel.X() + pos.Y()*vel.Y())/rho : 0.0;
    for(size_t itg = 0; itg < targets_.size(); ++itg){
      auto const& tg = targets_[itg];
      if(!tg.drives_ || !usable(itg)) continue;
      if(rho > tg.rmax_ && vrho > 0.0) continue;
      bool reach(false);
      switch(tg.reach_) {
        case SweepTarget::radial:
          reach = true;
          break;
        case SweepTarget::plane: {
          double vn = vel.Dot(tg.pl_->normal());
          double pd = (tg.pl_->center()-pos).Dot(tg.pl_->normal());
          reach = std::fabs(vn) >= minv_ && pd*vn > 0.0;
          break;
        }
        case SweepTarget::zplane: {
          double dz = tg.zval_ - pos.Z();
          if(dz*vel.Z() >= 0.0) {
            reach = true;
          } else {
            // near a reflection vz is an unreliable direction estimate: use the linearized trajectory instead
            KinKal::Plane zplane(KinKal::VEC3(0.0,0.0,1.0),KinKal::VEC3(1.0,0.0,0.0),KinKal::VEC3(0.0,0.0,tg.zval_));
            double dist(0.0);
            zplane.intersect(ktraj.linearize(time),dist,tdir_ == KinKal::TimeDir::forwards,intertol_);
            reach = dist*sgn > 0.0;
          }
          break;
        }
        case SweepTarget::zrange: {
          double vb = ktraj.bnom().R() > 0.0 ? vel.Dot(ktraj.bnom().Unit()) : vel.Z();
          reach = (pos.Z() > tg.zreglo_ && pos.Z() < tg.zreghi_) || (pos.Z() >= tg.zreghi_ && vb < 0.0) || (pos.Z() <= tg.zreglo_ && vb > 0.0);
          break;
        }
        case SweepTarget::inside:
          reach = rho < tg.rreg_ && pos.Z() > tg.zreglo_ && pos.Z() < tg.zreghi_;
          break;
      }
      if(reach){
        if(debug_ > 2) std::cout << "Sweep: " << tg.sid_ << " reachable from time " << time << " z " << pos.Z() << " rho " << rho << std::endl;
        return true;
      }
    }
    if(debug_ > 1) std::cout << "Sweep: nothing reachable from time " << time << " z " << pos.Z() << " rho " << rho << std::endl;
    return false;
  }
}
#endif
