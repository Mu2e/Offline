//
// Test intersections with KinKal objects
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/TimeRange.hh"
#include "Offline/RecoGeom/inc/Intersection.hh"
#include "Offline/RecoGeom/inc/Cylinder.hh"
#include "Offline/RecoGeom/inc/Disk.hh"
#include "Offline/RecoGeom/inc/Annulus.hh"
#include "Offline/RecoGeom/inc/Rectangle.hh"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

using VEC3 = ROOT::Math::XYZVectorD;
using mu2e::RecoGeom::Ray;
using mu2e::RecoGeom::Cylinder;
using mu2e::RecoGeom::Annulus;
using mu2e::RecoGeom::Rectangle;
using mu2e::RecoGeom::Disk;
using mu2e::RecoGeom::IntersectFlag;
using KinKal::ParticleState;
using KinKal::TimeRange;
static struct option long_options[] = {
  {"scost",     required_argument, 0, 'c'  },
  {"sphi",        required_argument, 0, 'p'  },
  {"slen1",        required_argument, 0, 'r'  },
  {"slen2",        required_argument, 0, 'l'  },
  {"pcost",     required_argument, 0, 'C'  },
  {"pphi",        required_argument, 0, 'P'  },
  {"pmom",        required_argument, 0, 'M'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: IntersectionTest  --slen1 --slen2 -(cylinder radius, half length, disk u, v 1/2 size, or Annulus inner/outer radius), -scost --sphi (surface axis direction)  --pmom, --pcost, --pphi  (particle momentum in MeV/c, direction angles) ");
}

int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  VEC3 point(0.0,0.0,0.0);
  double scost(1.0), sphi(0.0), slen1(400), slen2(1000);
  double pcost(0.5), pphi(0.0), pmom(400);
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'c' : scost = atof(optarg);
                 break;
      case 'p' : sphi = atof(optarg);
                 break;
      case 'r' : slen1 = atof(optarg);
                 break;
      case 'l' : slen2 = atof(optarg);
                 break;
      case 'C' : pcost = atof(optarg);
                 break;
      case 'P' : pphi = atof(optarg);
                 break;
      case 'M' : pmom = atof(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  VEC3 origin(0.0,0.0,0.0);

  double ssint = sqrt(1.0-scost*scost);
  VEC3 axis(ssint*cos(sphi), ssint*sin(sphi), scost);

  double psint = sqrt(1.0-pcost*pcost);
  VEC3 momvec(psint*cos(pphi), psint*sin(pphi), pcost);
  momvec *= pmom;
  ParticleState pstate(origin,momvec,0.0,0.5,-1);
//  std::cout << "Test " << pstate << std::endl;
  std::cout << "Test particle momentum " << pstate.momentum3() << std::endl;
  double speed = pstate.speed();
  double tmax = 2*sqrt(slen1*slen1 + slen2*slen2)/speed;

  VEC3 bnom(0.0,0.0,1.0);
  TimeRange trange(0.0,tmax);
  KinKal::KinematicLine ktraj(pstate,bnom,trange);
  // intersect with various surfaces
  Cylinder cyl(axis,origin,slen1,slen2);
  std::cout << "Test " << cyl << std::endl;

  mu2e::RecoGeom::Intersect<KinKal::KinematicLine, mu2e::RecoGeom::Cylinder> lc_isect;
  auto lc_inter = lc_isect.intersect(ktraj,cyl, 0.0, 1.0e-8);
  std::cout << "KinematicLine Cylinder Intersection status " << lc_inter.flag_ << " position " << lc_inter.pos_ << " time " << lc_inter.time_ << std::endl;

  Disk disk(axis,origin,slen1);
  std::cout << "Test " << disk << std::endl;

  mu2e::RecoGeom::Intersect<KinKal::KinematicLine, mu2e::RecoGeom::Disk> ld_isect;
  auto ld_inter = ld_isect.intersect(ktraj,disk, 0.0, 1.0e-8);
  std::cout << "KinematicLine Disk Intersection status " << ld_inter.flag_ << " position " << ld_inter.pos_ << " time " << ld_inter.time_ << std::endl;

  Annulus ann(axis,origin,slen1, slen2);
  std::cout << "Test " << ann << std::endl;

  mu2e::RecoGeom::Intersect<KinKal::KinematicLine, mu2e::RecoGeom::Annulus> la_isect;
  auto la_inter = la_isect.intersect(ktraj,ann, 0.0, 1.0e-8);
  std::cout << "KinematicLine Annulus Intersection status " << la_inter.flag_ << " position " << la_inter.pos_ << " time " << la_inter.time_ << std::endl;

  VEC3 udir(scost*cos(sphi), scost*sin(sphi), -ssint);
  Rectangle rect(axis,origin,udir,slen1, slen2);
  std::cout << "Test " << rect << std::endl;

  mu2e::RecoGeom::Intersect<KinKal::KinematicLine, mu2e::RecoGeom::Rectangle> lr_isect;
  auto lr_inter = lr_isect.intersect(ktraj,rect, 0.0, 1.0e-8);
  std::cout << "KinematicLine Rectangle Intersection status " << lr_inter.flag_ << " position " << lr_inter.pos_ << " time " << lr_inter.time_ << std::endl;
  return 0;
}
