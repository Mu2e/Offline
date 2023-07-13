//
// Test intersections with KinKal objects
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/TimeRange.hh"
#include "Offline/RecoGeom/inc/Intersection.hh"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

using VEC3 = ROOT::Math::XYZVectorD;
using mu2e::RecoGeom::Ray;
using mu2e::RecoGeom::Cylinder;
using mu2e::RecoGeom::Annulus;
using mu2e::RecoGeom::Rectangle;
using mu2e::RecoGeom::IntersectFlag;
using KinKal::ParticleState;
using KinKal::TimeRange;
static struct option long_options[] = {
  {"ccost",     required_argument, 0, 'c'  },
  {"cphi",        required_argument, 0, 'p'  },
  {"crad",        required_argument, 0, 'r'  },
  {"clen",        required_argument, 0, 'l'  },
  {"pcost",     required_argument, 0, 'C'  },
  {"pphi",        required_argument, 0, 'P'  },
  {"pmom",        required_argument, 0, 'M'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: IntersectionTest  --crad --clen --ccost --cphi (cylinder radius, half length, axis direction) --pmom, --pcost, --pphi  (particle momentum in MeV/c, direction) ");
}

int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  VEC3 point(0.0,0.0,0.0);
  double ccost(1.0), cphi(0.0), crad(400), clen(1000);
  double pcost(0.5), pphi(0.0), pmom(400);
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'c' : ccost = atof(optarg);
                 break;
      case 'p' : cphi = atof(optarg);
                 break;
      case 'r' : crad = atof(optarg);
                 break;
      case 'l' : clen = atof(optarg);
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

  double csint = sqrt(1.0-ccost*ccost);
  VEC3 axis(csint*cos(cphi), csint*sin(cphi), ccost);
  Cylinder cyl(axis,origin,crad,clen);
  std::cout << "Test " << cyl << std::endl;

  double psint = sqrt(1.0-pcost*pcost);
  VEC3 momvec(psint*cos(pphi), psint*sin(pphi), pcost);
  momvec *= pmom;
  ParticleState pstate(origin,momvec,0.0,0.5,-1);
//  std::cout << "Test " << pstate << std::endl;
  std::cout << "Test particle momentum " << pstate.momentum3() << std::endl;
  double speed = pstate.speed();
  double tmax = 2*sqrt(crad*crad + clen*clen)/speed;

  VEC3 bnom(0.0,0.0,1.0);
  TimeRange trange(0.0,tmax);
  KinKal::KinematicLine ktraj(pstate,bnom,trange);
  mu2e::RecoGeom::Intersect<KinKal::KinematicLine, mu2e::RecoGeom::Cylinder> isect;
  auto inter = isect.intersect(ktraj,cyl, 0.0, 1.0e-8);
  std::cout << "Intersection status " << inter.flag_ << " position " << inter.pos_ << " time " << inter.time_ << std::endl;

  return 0;
}
