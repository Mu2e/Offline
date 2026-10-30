// Geometry and material-insertion regression checks for the extrapolation sweep.
#include "Offline/Mu2eKinKal/inc/ExtrapolateSweep.hh"
#include "Offline/Mu2eKinKal/inc/KKShellXing.hh"
#include "KinKal/Fit/Material.hh"
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/SimpleFileFinder.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
  using KinKal::TimeDir;
  using KinKal::TimeRange;
  using KinKal::VEC3;
  using mu2e::ExtrapolateSweep;
  using mu2e::SweepTarget;
  using mu2e::SweepTargetCollection;
  constexpr double tol = 1.0e-4;

  void require(bool ok, std::string const& message) {
    if(!ok) throw std::runtime_error(message);
  }

  template <class KTRAJ> KinKal::Intersection crossing(KinKal::ParticleTrajectory<KTRAJ> const& traj,
      SweepTargetCollection const& targets, TimeDir dir, double step=1.0e6) {
    ExtrapolateSweep sweep(targets,dir,1.0e6,step,1.0e-4,tol,1.0e-6,0);
    sweep.start(dir == TimeDir::forwards ? traj.range().begin() : traj.range().end());
    sweep.limit(dir == TimeDir::forwards ? traj.range().end() : traj.range().begin(),false);
    require(!sweep.needsExtrapolation(traj,dir),"limited search requests extension");
    require(sweep.hasCrossing(),"no crossing found");
    return sweep.crossing().inter_;
  }

  void lineRoots(TimeDir dir, bool split, bool requireBounds) {
    bool forward = dir == TimeDir::forwards;
    KinKal::ParticleState state(VEC3(-200.0,0.0,forward ? -150.0 : -50.0),VEC3(100.0,0.0,50.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(0.0,1.0));
    double vx = line.velocity(0.0).X();
    line.setRange(TimeRange(0.0,(split ? 400.0 : 500.0)/vx));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    if(split){
      auto tail = line;
      tail.setRange(TimeRange(400.0/vx,500.0/vx));
      traj.append(tail);
    }
    KinKal::Cylinder cyl(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,0.0),100.0,50.0);
    SweepTarget target;
    target.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::DS_CryoInner);
    target.surf_ = &cyl;
    target.inbounds_ = requireBounds;
    auto inter = crossing(traj,{target},dir);
    double expected = requireBounds ? (forward ? 100.0 : -100.0) : (forward ? -100.0 : 100.0);
    require(std::fabs(inter.pos_.X()-expected) < 10.0*tol,"incorrect cylinder root");
    require(inter.inbounds_ == requireBounds,"unexpected root bounds");
  }

  template <class HELIX> void tiltedZ() {
    KinKal::ParticleState state(VEC3(10.0,20.0,30.0),VEC3(15.0,80.0,40.0),0.0,105.658,1);
    HELIX helix(state,VEC3(1.0,0.0,0.0),TimeRange(0.0,1.0));
    double period = 2.0*std::acos(-1.0)*helix.energy(0.0)/(KinKal::BFieldMap::cbar()*CLHEP::c_light);
    double zmax = -1.0e100, tmax = 0.0;
    for(unsigned i = 0; i <= 512; ++i){
      double time = period*i/512.0;
      double z = helix.position3(time).Z();
      if(z > zmax){ zmax = z; tmax = time; }
    }
    double lo = tmax-0.2*period, hi = tmax+0.2*period;
    helix.setRange(TimeRange(lo,hi));
    double zend = std::max(helix.position3(lo).Z(),helix.position3(hi).Z());
    double planeZ = 0.5*(zend+zmax);
    require(planeZ-zend > 100.0*tol,"tilted fixture lacks an interior z extremum");
    KinKal::Plane plane(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,planeZ));
    SweepTarget target;
    target.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::TT_Front);
    target.surf_ = &plane;
    target.pl_ = &plane;
    target.reach_ = SweepTarget::plane;
    KinKal::ParticleTrajectory<HELIX> traj(helix);
    for(auto dir : {TimeDir::forwards,TimeDir::backwards}){
      double a = dir == TimeDir::forwards ? lo : tmax;
      double b = dir == TimeDir::forwards ? tmax : hi;
      double fa = helix.position3(a).Z()-planeZ;
      for(unsigned iter = 0; iter < 80; ++iter){
        double mid = 0.5*(a+b);
        if((helix.position3(mid).Z()-planeZ)*fa > 0.0) a = mid;
        else b = mid;
      }
      double expected = 0.5*(a+b);
      for(double step : {period,0.025*period}){
        auto inter = crossing(traj,{target},dir,step);
        require(std::fabs(inter.time_-expected)*helix.speed(expected) < 10.0*tol,"incorrect tilted-helix crossing");
      }
    }
  }

  void stableTieAndSeams() {
    KinKal::ParticleState state(VEC3(0.0,0.0,-1.0),VEC3(0.0,0.0,100.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(0.0,0.1));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    KinKal::Plane plane(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,0.0));
    KinKal::Plane nearby(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,0.01));
    SweepTarget a, b;
    a.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::DS_Front);
    b.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::DS_Back);
    a.effect_ = b.effect_ = SweepTarget::planeshell;
    a.surf_ = b.surf_ = &plane;
    a.pl_ = b.pl_ = &plane;
    for(bool reverse : {false,true}){
      SweepTargetCollection targets = reverse ? SweepTargetCollection{b,a} : SweepTargetCollection{a,b};
      ExtrapolateSweep sweep(targets,TimeDir::forwards,1.0,1.0,1.0e-4,tol,1.0e-6,0);
      sweep.start(0.0);
      sweep.limit(0.1,false);
      sweep.needsExtrapolation(traj,TimeDir::forwards);
      require(sweep.hasCrossing(),"tie has no crossing");
      require(targets[sweep.crossing().itarget_].sid_ == a.sid_,"tie depends on target order");
      sweep.accept();
      sweep.needsExtrapolation(traj,TimeDir::forwards);
      require(!sweep.hasCrossing(),"coplanar seam material counted twice");
    }
    b.surf_ = &nearby;
    b.pl_ = &nearby;
    SweepTargetCollection targets{a,b};
    ExtrapolateSweep sweep(targets,TimeDir::forwards,1.0,1.0,1.0e-4,tol,1.0e-6,0);
    sweep.start(0.0);
    sweep.limit(0.1,false);
    sweep.needsExtrapolation(traj,TimeDir::forwards);
    require(sweep.hasCrossing(),"first plane absent");
    sweep.accept();
    sweep.needsExtrapolation(traj,TimeDir::forwards);
    require(sweep.hasCrossing(),"distinct nearby plane incorrectly seam-suppressed");
    require(targets[sweep.crossing().itarget_].sid_ == b.sid_,"nearby plane identity differs");
  }

  void boundaries(TimeDir dir, bool split) {
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,0.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(0.0,1.0));
    double dt = 100.0/line.speed(0.0);
    line.setRange(TimeRange(0.0,split ? dt : 2.0*dt));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    if(split){
      auto tail = line;
      tail.setRange(TimeRange(dt,2.0*dt));
      traj.append(tail);
    }
    KinKal::Plane plane(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(100.0,0.0,0.0));
    SweepTarget target;
    target.surf_ = &plane;
    target.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::TT_Front);
    auto inter = crossing(traj,{target},dir,dt);
    require(std::fabs(inter.time_-dt)*line.speed(0.0) < 10.0*tol,"window/piece boundary root lost");
    require(inter.good() && traj.nearestTraj(inter.time_)->range().inRange(inter.time_),"geometry boundary root has no owning piece");
    require((inter.pos_-traj.nearestTraj(inter.time_)->position3(inter.time_)).R() < tol,"geometry boundary root position differs");
  }

  void planarLimits(TimeDir dir) {
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,0.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(-1.0,1.0));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    KinKal::Plane plane(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,0.0));
    SweepTarget target;
    target.surf_ = target.pl_ = &plane;
    SweepTargetCollection targets{target};
    ExtrapolateSweep sweep(targets,dir,1.0,1.0,1.0e-4,tol,1.0e-6,0);
    sweep.start(dir == TimeDir::forwards ? -1.0 : 1.0);
    sweep.limit(0.0,false);
    require(!sweep.needsExtrapolation(traj,dir),"terminal-limit search requests extension");
    require(sweep.hasCrossing() == (dir == TimeDir::backwards),"planar limit violates half-open query range");
    if(sweep.hasCrossing()) require(sweep.crossing().inter_.time_ == 0.0,"included lower-limit root moved");
  }

  void solverLowerBoundary() {
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,0.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(0.0,1.0));
    KinKal::Plane plane(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,0.0));
    auto raw = KinKal::intersect(line,plane,line.range(),tol,TimeDir::backwards);
    require(raw.good() && raw.time_ == line.range().begin(),"solver fixture does not return an exact lower boundary");
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    SweepTarget target;
    target.surf_ = &plane;
    auto inter = crossing(traj,{target},TimeDir::backwards);
    require(inter.good() && inter.time_ == std::nextafter(0.0,1.0),"solver lower boundary lacks interior ownership");
    require(traj.nearestTraj(inter.time_) == traj.frontPtr(),"solver lower boundary has wrong piece");
  }

  void closeSeam(TimeDir dir) {
    double sign = dir == TimeDir::forwards ? 1.0 : -1.0;
    KinKal::ParticleState state(VEC3(0.0,0.0,-1.0),VEC3(0.0,0.0,100.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(0.0,0.02));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    KinKal::Plane first(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,0.0));
    KinKal::Plane second(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,sign*0.5*tol));
    SweepTarget a, b;
    a.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::DS_Front);
    b.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::DS_Back);
    a.effect_ = b.effect_ = SweepTarget::planeshell;
    a.surf_ = a.pl_ = &first;
    b.surf_ = b.pl_ = &second;
    SweepTargetCollection targets{a,b};
    ExtrapolateSweep sweep(targets,dir,1.0,1.0,1.0e-4,tol,1.0e-6,0);
    sweep.start(dir == TimeDir::forwards ? 0.0 : 0.02);
    sweep.limit(dir == TimeDir::forwards ? 0.02 : 0.0,false);
    sweep.needsExtrapolation(traj,dir);
    require(sweep.hasCrossing(),"close seam lacks first plane");
    sweep.accept();
    sweep.needsExtrapolation(traj,dir);
    require(!sweep.hasCrossing(),"close coplanar seam not suppressed");
  }

  SweepTarget seamTarget(KinKal::Plane const& plane, mu2e::SurfaceId::enum_type id, SweepTarget::Effect pool,
      MatEnv::DetMaterial const& material) {
    SweepTarget target;
    target.sid_ = mu2e::SurfaceId(id);
    target.effect_ = pool;
    target.surf_ = target.pl_ = &plane;
    target.mat_ = &material;
    target.thick_ = 0.01;
    return target;
  }

  void tangentEpisode(TimeDir dir, bool split, bool windows, SweepTarget::Effect pool,
      MatEnv::DetMaterial const& material, MatEnv::DetMaterial const& differentMaterial, unsigned control) {
    constexpr double seamTol = 0.1;
    double sign = dir == TimeDir::forwards ? 1.0 : -1.0;
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,10.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(-1.0,1.0));
    double vx = line.velocity(0.0).X(), end = sign*4.0/vx, begin = -end;
    double lo = std::min(begin,end), hi = std::max(begin,end), join = sign*0.6/vx;
    line.setRange(TimeRange(lo,split ? join : hi));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    if(split){
      auto tail = line;
      tail.setRange(TimeRange(join,hi));
      traj.append(tail);
    }
    KinKal::Rectangle first(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,0.0),control == 5 ? 0.05 : 10.0,10.0);
    KinKal::Rectangle second(VEC3(0.0,0.0,-1.0),VEC3(1.0,0.0,0.0),VEC3(sign*0.5,0.0,sign*0.05),10.0,10.0);
    KinKal::Rectangle tilted(VEC3(0.02,0.0,1.0),VEC3(1.0,0.0,-0.02),VEC3(0.0,0.0,0.0),10.0,10.0);
    KinKal::Rectangle separate(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,sign*0.2),10.0,10.0);
    KinKal::Rectangle nonoverlap(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(sign*2.2,0.0,sign*0.05),1.8,10.0);
    KinKal::Plane intervening(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(sign*0.8,0.0,0.0));
    auto a = seamTarget(first,mu2e::SurfaceIdEnum::DS_Front,pool,material);
    auto b = seamTarget(second,mu2e::SurfaceIdEnum::DS_Back,pool,material);
    auto c = seamTarget(intervening,mu2e::SurfaceIdEnum::TT_Front,pool,material);
    if(control == 1) b.mat_ = &differentMaterial;
    if(control == 2) b.thick_ *= 2.0;
    if(control == 3) b.surf_ = b.pl_ = &tilted;
    if(control == 4) b.surf_ = b.pl_ = &separate;
    if(control == 5) b.surf_ = b.pl_ = &nonoverlap;
    if(control == 6) b.effect_ = pool == SweepTarget::planeshell ? SweepTarget::crvsector : SweepTarget::planeshell;
    if(control == 7) a.maxcross_ = 0;
    SweepTargetCollection targets{a,b,c};
    double step = windows ? 0.15/vx : hi-lo;
    ExtrapolateSweep sweep(targets,dir,1.0,step,1.0e-4,seamTol,1.0e-6,0);
    sweep.start(0.0);
    sweep.limit(end,false);
    require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing(),"tangent episode lacks first crossing");
    require(sweep.crossing().itarget_ == 0,"tangent first target differs");
    double accepted = sweep.crossing().inter_.time_;
    require(std::fabs(first.distance(sweep.crossing().inter_.pos_)) < seamTol,"tangent first point off plane");
    sweep.accept();
    sweep.start(accepted);
    if(control == 0 || control == 7){
      sweep.limit(sign*0.15/vx,false);
      require(!sweep.needsExtrapolation(traj,dir) && !sweep.hasCrossing(),"tangent seam escaped a short window");
      sweep.start(sign*0.15/vx);
      sweep.limit(end,false);
      require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing(),"intervening plane absent");
      require(sweep.crossing().itarget_ == 2,"tangent seam escaped before different-plane accept");
      accepted = sweep.crossing().inter_.time_;
      sweep.accept();
      sweep.start(accepted);
      require(!sweep.needsExtrapolation(traj,dir) && !sweep.hasCrossing(),"episode lost across start/piece/window/different-plane accept");
    } else {
      bool neighbor = false;
      for(unsigned trial = 0; trial < 3; ++trial){
        require(!sweep.needsExtrapolation(traj,dir),"control requests extrapolation");
        if(!sweep.hasCrossing()) break;
        if(sweep.crossing().itarget_ == 1) neighbor = true;
        sweep.accept();
      }
      require(neighbor,"ineligible seam neighbor incorrectly suppressed");
    }
  }

  void helixEpisodes(TimeDir dir, bool split, bool windows, bool self, MatEnv::DetMaterial const& material) {
    double sign = dir == TimeDir::forwards ? 1.0 : -1.0;
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,20.0),0.0,105.658,1);
    KinKal::LoopHelix helix(state,VEC3(0.0,0.0,1.0),TimeRange(-1.0,1.0));
    double period = 2.0*std::acos(-1.0)*helix.energy(0.0)/(KinKal::BFieldMap::cbar()*CLHEP::c_light);
    double begin = -sign*0.1*period, end = sign*1.2*period;
    double lo = std::min(begin,end), hi = std::max(begin,end), join = sign*0.25*period;
    helix.setRange(TimeRange(lo,split ? join : hi));
    KinKal::ParticleTrajectory<KinKal::LoopHelix> traj(helix);
    if(split){
      auto tail = helix;
      tail.setRange(TimeRange(join,hi));
      traj.append(tail);
    }
    KinKal::Rectangle first(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,0.0),2000.0,2000.0);
    KinKal::Rectangle second(VEC3(-1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(0.0,0.0,10.0),2000.0,2000.0);
    auto a = seamTarget(first,mu2e::SurfaceIdEnum::DS_Front,SweepTarget::planeshell,material);
    auto b = seamTarget(second,mu2e::SurfaceIdEnum::DS_Back,SweepTarget::planeshell,material);
    a.maxcross_ = self ? 0 : 1;
    SweepTargetCollection targets = self ? SweepTargetCollection{a} : SweepTargetCollection{a,b};
    ExtrapolateSweep sweep(targets,dir,100.0,windows ? 0.13*period : hi-lo,1.0e-4,0.1,1.0e-6,0);
    sweep.start(begin);
    sweep.limit(end,false);
    for(unsigned icross = 0; icross < (self ? 3u : 2u); ++icross){
      require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing(),"depart-return helix crossing absent");
      auto const inter = sweep.crossing().inter_;
      require(sweep.crossing().itarget_ == (self ? 0u : icross),"neighbor quota consumed by suppressed seam");
      require(std::fabs(inter.time_-sign*0.5*period*icross)*helix.speed(0.0) < 1.0,"helix recross time differs");
      if(icross > 0){
        auto pos = traj.position3(sign*(0.5*icross-0.25)*period);
        require(std::fabs(first.distance(pos)) > 0.1 && std::fabs(second.distance(pos)) > 0.1,"helix lacks out-of-both departure");
      }
      sweep.accept();
      sweep.start(inter.time_);
    }
    require(!sweep.needsExtrapolation(traj,dir) && !sweep.hasCrossing(),"helix episode/quota retained an extra crossing");
  }

  template <class HELIX> void shallowEpisode(TimeDir dir, bool split, bool windows, bool reverse, MatEnv::DetMaterial const& material) {
    double sign = dir == TimeDir::forwards ? 1.0 : -1.0;
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,20.0),0.0,105.658,1);
    HELIX helix(state,VEC3(0.0,0.0,1.0),TimeRange(-100.0,100.0));
    double omega = KinKal::BFieldMap::cbar()*CLHEP::c_light/helix.energy(), pi = std::acos(-1.0);
    double maximum = sign*0.5*pi/omega;
    double radius = std::fabs(helix.position3(maximum).X());
    double halfwidth = std::acos(1.0-0.15/radius)/omega;
    double begin = maximum-sign*halfwidth, end = maximum+sign*1.5*halfwidth;
    double lo = std::min(begin-sign*halfwidth,end), hi = std::max(begin-sign*halfwidth,end);
    helix.setRange(TimeRange(lo,split ? maximum : hi));
    KinKal::ParticleTrajectory<HELIX> traj(helix);
    if(split){
      auto tail = helix;
      tail.setRange(TimeRange(maximum,hi));
      traj.append(tail);
    }
    KinKal::Rectangle plane(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),helix.position3(begin),2000.0,2000.0);
    auto a = seamTarget(plane,mu2e::SurfaceIdEnum::DS_Front,SweepTarget::planeshell,material);
    auto b = seamTarget(plane,mu2e::SurfaceIdEnum::DS_Back,SweepTarget::planeshell,material);
    SweepTargetCollection targets = reverse ? SweepTargetCollection{b,a} : SweepTargetCollection{a,b};
    ExtrapolateSweep sweep(targets,dir,100.0,windows ? 0.3*halfwidth : 100.0,1.0e-4,0.1,1.0e-6,0);
    sweep.start(begin);
    sweep.limit(end,false);
    require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing() && targets[sweep.crossing().itarget_].sid_ == a.sid_,"shallow episode lacks first crossing");
    sweep.accept();
    require(std::fabs(plane.distance(helix.position3(maximum))) > 0.1 &&
        std::fabs(plane.distance(helix.position3(maximum))) < 0.2,"shallow fixture does not resolve only one tolerance");
    require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing() && targets[sweep.crossing().itarget_].sid_ == b.sid_,"shallow resolved depart-return was suppressed");
    require(sign*(sweep.crossing().inter_.time_-maximum) > 0.0,"shallow return accepted before resolved departure");
    require(std::fabs(plane.distance(sweep.crossing().inter_.pos_)) < 0.1,"shallow return outside plane tolerance");
    // An outgoing query must find the later root without requiring a preceding accepted seam.
    targets.front().maxcross_ = targets.back().maxcross_ = 1;
    ExtrapolateSweep outgoing(targets,dir,100.0,windows ? 0.3*halfwidth : 100.0,1.0e-4,0.1,1.0e-6,0);
    double queryStart = maximum-sign*0.5*halfwidth;
    require(std::fabs(plane.distance(helix.position3(queryStart))) > 0.1,"outgoing query starts inside tolerance slab");
    outgoing.start(queryStart);
    outgoing.limit(end,false);
    require(!outgoing.needsExtrapolation(traj,dir) && outgoing.hasCrossing(),"outgoing query lost bracketed later plane root");
    auto const& inter = outgoing.crossing().inter_;
    double expected = maximum+sign*halfwidth;
    require(sign*(inter.time_-queryStart) > 0.0 && sign*(inter.time_-maximum) > 0.0,"outgoing query returned earlier plane root");
    if(!windows) require(std::fabs(inter.time_-expected)*helix.speed(inter.time_) < 0.1,"bracketed shallow return time differs");
    require(plane.onSurface(inter.pos_,0.1),"outgoing return is outside plane tolerance");
    require((inter.pdir_-helix.direction(inter.time_)).R() < 1.0e-12,"backwards root reversed physical momentum direction");
  }

  template <class HELIX> void planeContacts(TimeDir dir) {
    double sign = dir == TimeDir::forwards ? 1.0 : -1.0, tolerance = 0.1;
    KinKal::ParticleState state(VEC3(10.0,20.0,30.0),VEC3(100.0,0.0,20.0),0.0,105.658,1);
    HELIX helix(state,VEC3(0.0,0.0,1.0),TimeRange(-100.0,100.0));
    double omega = KinKal::BFieldMap::cbar()*CLHEP::c_light/helix.energy();
    double peak = sign*0.5*std::acos(-1.0)/omega, width = 0.05/omega;
    helix.setRange(TimeRange(peak-width,peak+width));
    KinKal::ParticleTrajectory<HELIX> traj(helix);
    for(double offset : {0.0,0.5*tolerance,1.5*tolerance}){
      auto pos = helix.position3(peak)+VEC3(sign*offset,0.0,0.0);
      KinKal::Plane plane(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),pos);
      SweepTarget target;
      target.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::TT_Front);
      target.surf_ = target.pl_ = &plane;
      for(unsigned edge : {0u,1u}){
        for(double mincos : {0.0,0.01}){
          target.mincos_ = mincos;
          SweepTargetCollection targets{target};
          ExtrapolateSweep sweep(targets,dir,100.0,100.0,1.0e-4,tolerance,1.0e-6,0);
          double begin = edge == 1 ? peak : helix.range().begin();
          double end = helix.range().end();
          sweep.start(sign > 0.0 ? begin : end);
          sweep.limit(sign > 0.0 ? end : begin,false);
          require(!sweep.needsExtrapolation(traj,dir),"contact search requests extension");
          bool expected = offset < tolerance && mincos == 0.0;
          require(sweep.hasCrossing() == expected,"plane contact/incidence policy differs");
          if(expected){
            auto const& inter = sweep.crossing().inter_;
            require(plane.onSurface(inter.pos_,tolerance),"contact outside plane tolerance");
            require(std::fabs(inter.time_-peak)*helix.speed(peak) < tolerance,"contact time is not the plane-distance extremum");
          }
        }
      }
    }
  }

  template <class HELIX> void longPlaneQuery(TimeDir dir) {
    double sign = dir == TimeDir::forwards ? 1.0 : -1.0;
    KinKal::ParticleState state(VEC3(10.0,20.0,30.0),VEC3(100.0,0.0,20.0),0.0,105.658,1);
    HELIX helix(state,VEC3(0.0,0.0,1.0),TimeRange(-100.0,100.0));
    double period = 2.0*std::acos(-1.0)*helix.energy()/(KinKal::BFieldMap::cbar()*CLHEP::c_light);
    helix.setRange(TimeRange(sign > 0.0 ? 0.0 : -300.0*period,sign > 0.0 ? 300.0*period : 0.0));
    KinKal::ParticleTrajectory<HELIX> traj(helix);
    KinKal::Plane plane(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),helix.position3(sign*0.125*period));
    SweepTarget target;
    target.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::TT_Front);
    target.surf_ = target.pl_ = &plane;
    auto inter = crossing(traj,{target},dir);
    require(std::fabs(inter.time_-sign*0.125*period)*helix.speed(inter.time_) < tol,"long query lost early bracketed crossing");
  }

  void speculativeEpisode(TimeDir dir, MatEnv::DetMaterial const& material) {
    double sign = dir == TimeDir::forwards ? 1.0 : -1.0;
    auto makeLine = [](VEC3 const& pos, VEC3 const& momentum, double time){
      KinKal::ParticleState state(pos,momentum,time,105.658,1);
      return KinKal::KinematicLine(state,VEC3(0.0,0.0,0.0),TimeRange(time-1.0,time+1.0));
    };
    auto first = makeLine(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,2.0),0.0);
    double t1 = sign*2.0/first.velocity(0.0).X();
    auto middle = makeLine(first.position3(t1),VEC3(100.0,0.0,46.0),t1);
    double t2 = t1+sign/middle.velocity(t1).X();
    auto last = makeLine(middle.position3(t2),VEC3(100.0,0.0,-50.0),t2);
    double end = t2+sign*1.2/last.velocity(t2).X();
    first.setRange(TimeRange(std::min(0.0,t1),std::max(0.0,t1)));
    middle.setRange(TimeRange(std::min(t1,t2),std::max(t1,t2)));
    last.setRange(TimeRange(std::min(t2,end),std::max(t2,end)));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj;
    if(dir == TimeDir::forwards){ traj.append(first); traj.append(middle); traj.append(last); }
    else { traj.append(last); traj.append(middle); traj.append(first); }
    KinKal::Rectangle plane(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,0.0),10.0,10.0);
    KinKal::Plane other(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(sign*0.8,0.0,0.0));
    auto a = seamTarget(plane,mu2e::SurfaceIdEnum::DS_Front,SweepTarget::planeshell,material);
    auto b = seamTarget(plane,mu2e::SurfaceIdEnum::DS_Back,SweepTarget::planeshell,material);
    auto c = seamTarget(other,mu2e::SurfaceIdEnum::TT_Front,SweepTarget::planeshell,material);
    SweepTargetCollection targets{a,b,c};
    ExtrapolateSweep sweep(targets,dir,1.0,1.0,1.0e-4,0.1,1.0e-6,0);
    sweep.start(0.0);
    sweep.limit(end,false);
    require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing() && sweep.crossing().itarget_ == 0,
        "speculative fixture lacks initial seam crossing");
    double initial = sweep.crossing().inter_.time_;
    sweep.accept();
    sweep.setActive(2,false);
    require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing() && sweep.crossing().itarget_ == 1 &&
        sign*(sweep.crossing().inter_.time_-t2) > 0.0,"speculative fixture lacks future return candidate");
    sweep.setActive(2,true);
    sweep.start(initial);
    require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing() && sweep.crossing().itarget_ == 2,
        "earlier different plane did not win over future return");
    double accepted = sweep.crossing().inter_.time_;
    sweep.accept();
    auto changed = first;
    changed.setRange(TimeRange(std::min(0.0,end),std::max(0.0,end)));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> updated(changed);
    require(std::fabs(plane.distance(updated.position3(end))) < 0.1,"updated trajectory left the seam slab");
    sweep.start(accepted);
    require(!sweep.needsExtrapolation(updated,dir) && !sweep.hasCrossing(),
        "speculative future witness released an episode before trajectory mutation");
  }

  KinKal::LoopHelix materialPiece(double join, TimeRange const& range, VEC3 const& momentum) {
    KinKal::ParticleState state(VEC3(10.0,20.0,30.0),momentum,join,105.658,1);
    KinKal::DMAT covariance;
    for(unsigned row = 0; row < 6; ++row)
      for(unsigned col = 0; col <= row; ++col) covariance(row,col) = row == col ? 0.01 : 0.0;
    return KinKal::LoopHelix(KinKal::ParticleStateEstimate(state,covariance),VEC3(0.0,0.0,1.0),range);
  }

  void finiteState(KinKal::LoopHelix const& piece, double time) {
    auto state = piece.stateEstimate(time);
    for(unsigned row = 0; row < 6; ++row){
      require(std::isfinite(state.state()[row]),"nonfinite material state");
      require(state.stateCovariance()(row,row) > 0.0,"material state has nonpositive variance");
      for(unsigned col = 0; col < 6; ++col)
        require(std::isfinite(state.stateCovariance()(row,col)),"nonfinite material state covariance");
    }
  }

  void sameIntersection(KinKal::Intersection const& a, KinKal::Intersection const& b) {
    require(a.time_ == b.time_ && a == b && a.gap_ == b.gap_,"Material constructor changed crossing time/flags");
    require((a.pos_-b.pos_).R() == 0.0 && (a.pdir_-b.pdir_).R() == 0.0 && (a.norm_-b.norm_).R() == 0.0,
        "Material constructor changed crossing geometry");
  }

  void materialPieceEdge(TimeDir dir, unsigned npiece, double join, MatEnv::DetMaterial const& targetMaterial,
      bool boundary = true) {
    using Helix = KinKal::LoopHelix;
    using Shell = mu2e::KKShellXing<Helix,KinKal::Annulus>;
    bool forward = dir == TimeDir::forwards;
    unsigned joinIndex = forward ? npiece-1 : 1;
    KinKal::ParticleTrajectory<Helix> traj;
    for(unsigned ip = 0; ip < npiece; ++ip){
      double begin = join+(static_cast<int>(ip)-static_cast<int>(joinIndex))*0.1;
      double end = join+(static_cast<int>(ip)+1-static_cast<int>(joinIndex))*0.1;
      bool selected = ip == (forward ? npiece-1 : 0);
      auto piece = materialPiece(join,TimeRange(begin,end),
          selected ? VEC3(60.0,30.0,100.0) : VEC3(65.0+5.0*ip,35.0,110.0+5.0*ip));
      traj.append(piece);
    }
    // Shared nearestIndex assigns an exact join to the earlier piece.
    for(unsigned ip = 1; ip < npiece; ++ip){
      double edge = traj.piece(ip).range().begin();
      require(traj.nearestIndex(edge) == ip-1,"nearestIndex no longer chooses earlier join piece");
      require(traj.nearestTraj(edge) == traj.indexTraj(ip-1),"exact join pointer differs");
    }
    auto ref = forward ? traj.backPtr() : traj.frontPtr();
    auto neighbor = traj.indexTraj(forward ? npiece-2 : 1);
    require(ref != neighbor && (ref->momentum3(join)-neighbor->momentum3(join)).R() > 1.0,
        "material fixture does not distinguish neighboring pieces");
    double root = boundary ? join : join+(forward ? 0.05 : -0.05);
    double expected = boundary ? std::nextafter(join,forward ? ref->range().end() : ref->range().begin()) : root;
    auto annulus = std::make_shared<KinKal::Annulus>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),ref->position3(root),0.0,100.0);
    SweepTarget target;
    target.sid_ = mu2e::SurfaceId(mu2e::SurfaceIdEnum::TT_Front);
    target.effect_ = SweepTarget::stfoil;
    target.surf_ = target.pl_ = annulus.get();
    target.ann_ = annulus;
    target.mat_ = &targetMaterial;
    target.thick_ = 0.1056;
    target.maxcross_ = 1;
    SweepTargetCollection targets{target};
    ExtrapolateSweep sweep(targets,dir,1.0,1.0,1.0e-4,tol,1.0e-6,0);
    sweep.start(root);
    sweep.limit(forward ? ref->range().end() : std::nextafter(ref->range().begin(),ref->range().end()),false);
    require(!sweep.needsExtrapolation(traj,dir) && sweep.hasCrossing(),"piece-edge material crossing absent");
    auto inter = sweep.crossing().inter_;
    require(inter.time_ == expected,"piece-edge root is not exactly one ULP interior (or interior root moved)");
    require(inter.good() && !inter.gap_,"piece-edge crossing flags invalid");
    require(inter.time_ > ref->range().begin() && inter.time_ < ref->range().end(),"crossing not strictly inside end piece");
    require(traj.nearestTraj(inter.time_) == ref && ref->range().inRange(inter.time_),"crossing has wrong owning piece");
    require(std::fabs(inter.time_-root)*ref->speed(root) <= tol,"normalized time changed physical root");
    require((inter.pos_-ref->position3(inter.time_)).R() == 0.0 &&
        (inter.pdir_-ref->direction(inter.time_)).R() == 0.0 &&
        (inter.norm_-annulus->normal(inter.pos_)).R() == 0.0,"normalized crossing geometry is stale");
    require(annulus->onSurface(inter.pos_,tol) && annulus->inBounds(inter.pos_,tol),"normalized crossing left annulus");
    auto shell = std::make_shared<Shell>(annulus,target.sid_,targetMaterial,inter,ref,target.thick_,tol);
    require(&shell->referenceTrajectory() == ref.get() && shell->active(),"shell lost owning piece or active material");
    KinKal::Material<Helix> material(shell,traj);
    // Constructor reference lookup must retain the accepted root and source piece.
    require(material.elementXingPtr() == shell && &material.referenceTrajectory() == ref.get(),"Material constructor changed reference pointer");
    sameIntersection(shell->intersection(),inter);
    require(material.time() == expected && material.active() && shell->matXings().size() == 1,"Material constructor lost active crossing");
    require(std::isfinite(shell->matXings().front().plen_) && shell->matXings().front().plen_ > 0.0,"invalid shell material path");
    finiteState(*ref,inter.time_);
    KinKal::SVEC3 dmom;
    KinKal::SMAT dmomvar;
    shell->momentumChange(dmom,dmomvar);
    double dm, parallel, perpendicular;
    shell->materialEffects(dm,parallel,perpendicular);
    require(std::isfinite(dm) && dm < 0.0 && std::isfinite(parallel) && parallel >= 0.0 &&
        std::isfinite(perpendicular) && perpendicular > 0.0,"material effects inactive or nonfinite");
    for(unsigned row = 0; row < 3; ++row){
      require(std::isfinite(dmom[row]),"nonfinite material momentum change");
      for(unsigned col = 0; col < 3; ++col) require(std::isfinite(dmomvar(row,col)),"nonfinite material momentum variance");
    }
    auto originalRange = traj.range();
    auto originalPieces = traj.pieces();
    VEC3 expectedMomentum = ref->momentum3(inter.time_)+(forward ? 1.0 : -1.0)*VEC3(dmom[0],dmom[1],dmom[2]);
    material.extrapolate(traj,dir);
    require(traj.pieces().size() == npiece+1,"material extrapolation did not add one piece");
    require(traj.range().begin() == originalRange.begin() && traj.range().end() == originalRange.end(),"material extrapolation changed total range");
    for(unsigned ip = 0; ip < traj.pieces().size(); ++ip){
      require(traj.indexTraj(ip) && !traj.piece(ip).range().null(),"material insertion produced null piece/range");
      if(ip > 0) require(traj.piece(ip-1).range().end() == traj.piece(ip).range().begin(),"material insertion ranges not contiguous");
    }
    for(unsigned ip = 0; ip < npiece; ++ip)
      require(traj.indexTraj(ip+(forward ? 0 : 1)) == originalPieces[ip],"material insertion replaced existing piece pointer");
    auto added = forward ? traj.backPtr() : traj.frontPtr();
    require(added != ref && (forward ? added->range().begin() : added->range().end()) == expected,"new material piece boundary differs");
    require((added->position3(expected)-inter.pos_).R() < tol &&
        (added->momentum3(expected)-expectedMomentum).R() < tol,"material extrapolation state differs");
    finiteState(*added,expected);
    sameIntersection(shell->intersection(),inter);
    require(&shell->referenceTrajectory() == ref.get(),"extrapolation changed shell reference");
    sweep.accept();
    require(!sweep.hasCrossing(),"accept retained cached crossing");
    require(!sweep.needsExtrapolation(traj,dir) && !sweep.hasCrossing(),"maxcross=1 material re-found after accept");
  }

  void invalidBoundaryInterior(TimeDir dir, bool narrow) {
    bool forward = dir == TimeDir::forwards;
    double begin = narrow ? 90208.803182193456 : std::ldexp(1.0,40);
    double end = narrow ? std::nextafter(begin,std::numeric_limits<double>::infinity()) : begin+0.1;
    double root = forward ? begin : end;
    KinKal::ParticleState state(VEC3(10.0,20.0,30.0),VEC3(60.0,30.0,100.0),root,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(begin,end));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    KinKal::Plane plane(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),line.position3(root));
    SweepTarget target;
    target.surf_ = &plane;
    SweepTargetCollection targets{target};
    ExtrapolateSweep sweep(targets,dir,1.0,1.0,1.0e-4,tol,1.0e-6,0);
    sweep.start(root);
    sweep.limit(forward ? end : begin,false);
    bool failed = false;
    try {
      sweep.needsExtrapolation(traj,dir);
    } catch(std::runtime_error const& error) {
      auto expected = narrow ? "no representable piece interior" : "displacement exceeds intersection tolerance";
      failed = std::string(error.what()).find(expected) != std::string::npos;
    }
    require(failed && !sweep.hasCrossing(),"invalid boundary interior must fail explicitly without accepting geometry");
  }

  // The per-window reach bound must keep crossings just before a window end and in later windows.
  void reachBound(TimeDir dir) {
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(0.0,0.0,100.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(0.0,10.0));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    bool forward = dir == TimeDir::forwards;
    for(double dt : {0.999999,1.5,3.25}){
      double tcross = forward ? dt : 10.0-dt;
      KinKal::Plane plane(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),line.position3(tcross));
      SweepTarget target;
      target.surf_ = target.pl_ = &plane;
      auto inter = crossing(traj,{target},dir,1.0);
      require(std::fabs(inter.time_-tcross)*line.speed(tcross) <= tol,"reach bound dropped a crossing near a window end");
    }
  }

  // A crossing ranges out forwards in time exactly when its mean energy loss reaches the kinetic energy, never backwards.
  void rangeOut(MatEnv::DetMaterial const& material) {
    using Helix = KinKal::LoopHelix;
    using Shell = mu2e::KKShellXing<Helix,KinKal::Plane>;
    double const mass = 105.658, thick = 10.0;
    for(double pmag : {30.0,1000.0}){
      KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(0.0,0.2*pmag,pmag),0.0,mass,1);
      auto helix = std::make_shared<Helix>(state,VEC3(0.0,0.0,1.0),TimeRange(0.0,1.0));
      auto plane = std::make_shared<KinKal::Plane>(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),helix->position3(0.5));
      auto inter = KinKal::intersect(*helix,*plane,helix->range(),tol);
      require(inter.good(),"range-out fixture has no crossing");
      Shell shell(plane,mu2e::SurfaceId(mu2e::SurfaceIdEnum::TT_Front),material,inter,helix,thick,tol);
      double p = helix->momentum(), ke = helix->energy()-mass;
      double eloss = std::fabs(material.energyLoss(p,shell.matXings().front().plen_,mass));
      bool stops = eloss >= ke;
      require(stops == (pmag < 100.0),"range-out fixture does not bracket the kinetic energy");
      require(mu2e::rangesOut(shell,TimeDir::forwards) == stops,"forwards range-out disagrees with energy loss vs kinetic energy");
      require(!mu2e::rangesOut(shell,TimeDir::backwards),"backwards extrapolation must never range out");
    }
  }

  void retryBudget() {
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,20.0),0.0,105.658,1);
    KinKal::LoopHelix helix(state,VEC3(0.0,0.0,1.0),TimeRange(0.0,1.0));
    double period = 2.0*std::acos(-1.0)/std::fabs(helix.omega());
    helix.setRange(TimeRange(0.0,300.0*period));
    KinKal::ParticleTrajectory<KinKal::LoopHelix> traj(helix);
    KinKal::Plane plane(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),VEC3(2.0*helix.bendRadius(),0.0,0.0));
    SweepTarget target;
    target.surf_ = target.pl_ = &plane;
    bool failed = false;
    try {
      crossing(traj,{target},TimeDir::forwards);
    } catch(std::runtime_error const& error) {
      failed = std::string(error.what()).find("extrema budget") != std::string::npos;
    }
    require(failed,"unresolved long plane search must fail explicitly at the extrema budget");
  }

  void emptySearch() {
    KinKal::ParticleState state(VEC3(0.0,0.0,0.0),VEC3(0.0,0.0,100.0),0.0,105.658,1);
    KinKal::KinematicLine line(state,VEC3(0.0,0.0,0.0),TimeRange(0.0,0.1));
    KinKal::ParticleTrajectory<KinKal::KinematicLine> traj(line);
    SweepTargetCollection targets;
    ExtrapolateSweep sweep(targets,TimeDir::forwards,1.0,1.0,1.0e-4,tol,1.0e-6,0);
    sweep.start(0.0);
    sweep.limit(0.1,false);
    require(!sweep.needsExtrapolation(traj,TimeDir::forwards) && !sweep.hasCrossing(),"empty target search");
  }
}

int main() {
  try {
    auto workdir = std::getenv("MUSE_WORK_DIR");
    require(workdir && *workdir,"MUSE_WORK_DIR is required for TrackerConditions material data");
    MatEnv::SimpleFileFinder finder("MUSE_WORK_DIR","/Offline/TrackerConditions/data/");
    MatEnv::MatDBInfo materials(finder);
    auto targetMaterial = materials.findDetMaterial("Target");
    require(static_cast<bool>(targetMaterial),"Target material unavailable");
    auto differentMaterial = materials.findDetMaterial("Mylar");
    require(static_cast<bool>(differentMaterial),"Mylar material unavailable");
    for(auto dir : {TimeDir::forwards,TimeDir::backwards}){
      for(unsigned npiece : {2u,3u})
        for(double join : {90208.803182193456,-90208.803182193456,0.0}){
          materialPieceEdge(dir,npiece,join,*targetMaterial);
          materialPieceEdge(dir,npiece,join,*targetMaterial,false);
        }
      reachBound(dir);
      invalidBoundaryInterior(dir,true);
      invalidBoundaryInterior(dir,false);
      closeSeam(dir);
      planarLimits(dir);
      planeContacts<KinKal::LoopHelix>(dir);
      planeContacts<KinKal::CentralHelix>(dir);
      longPlaneQuery<KinKal::LoopHelix>(dir);
      longPlaneQuery<KinKal::CentralHelix>(dir);
      speculativeEpisode(dir,*targetMaterial);
      for(bool split : {false,true}){
        for(bool windows : {false,true}){
          for(bool reverse : {false,true}){
            shallowEpisode<KinKal::LoopHelix>(dir,split,windows,reverse,*targetMaterial);
            shallowEpisode<KinKal::CentralHelix>(dir,split,windows,reverse,*targetMaterial);
          }
          for(auto pool : {SweepTarget::planeshell,SweepTarget::crvsector,SweepTarget::stfoil}){
            tangentEpisode(dir,split,windows,pool,*targetMaterial,*differentMaterial,0);
            tangentEpisode(dir,split,windows,pool,*targetMaterial,*differentMaterial,7);
          }
          helixEpisodes(dir,split,windows,false,*targetMaterial);
          helixEpisodes(dir,split,windows,true,*targetMaterial);
        }
      }
      for(unsigned control = 1; control <= 6; ++control)
        for(bool windows : {false,true}) tangentEpisode(dir,true,windows,SweepTarget::planeshell,*targetMaterial,*differentMaterial,control);
      for(bool split : {false,true}){
        boundaries(dir,split);
        for(bool bounds : {false,true}) lineRoots(dir,split,bounds);
      }
    }
    rangeOut(*targetMaterial);
    solverLowerBoundary();
    tiltedZ<KinKal::CentralHelix>();
    tiltedZ<KinKal::LoopHelix>();
    stableTieAndSeams();
    retryBudget();
    emptySearch();
    std::cout << "ExtrapolateSweepTest: all checks passed" << std::endl;
    return 0;
  } catch(std::exception const& error) {
    std::cerr << "ExtrapolateSweepTest: " << error.what() << std::endl;
    return 1;
  }
}
