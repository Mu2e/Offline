//#include "Offline/TrkHitReco/inc/StereoLine.hh"
#include "Offline/TrkHitReco/inc/StereoPoint.hh"
#include "Offline/TrkHitReco/inc/CombineStereoPoints.hh"
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <getopt.h>

using mu2e::StereoPoint;
//using mu2e::StereoLine;
using mu2e::CombineStereoPoints;
using VEC3 = ROOT::Math::XYZVectorF;

static struct option long_options[] = {
  {"npts",     required_argument, 0, 'n' },
  {"dz",     required_argument, 0, 'z' },
  {"px",     required_argument, 0, 'X'  },
  {"py",     required_argument, 0, 'Y'  },
  {"rx",     required_argument, 0, 'x'  },
  {"ry",     required_argument, 0, 'y'  },
  {"ures",     required_argument, 0, 'u'  },
  {"vres",     required_argument, 0, 'v'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: StereoLineTest  --npts --dz (z range) --px --py (point) --rx --ry (slope)--ures --vres (transverse resolutions)");
}

int main(int argc, char** argv) {

  int opt;
  unsigned npts(5);
  VEC3 pos, dir;
  int long_index =0;
  double ures(1.0), vres(1.0), deltaz(10.0);
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'x' : dir.SetX(atof(optarg));
                 break;
      case 'X' : pos.SetX(atof(optarg));
                 break;
      case 'y' : dir.SetY(atof(optarg));
                 break;
      case 'Y' : pos.SetY(atof(optarg));
                 break;
      case 'z' : deltaz = atof(optarg);
                 break;
      case 'u' : ures = atof(optarg);
                 break;
      case 'v' : vres = atof(optarg);
                 break;
      case 'n' : npts = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }

  dir.SetZ(1.0);
  dir = dir.Unit();
  double dz = deltaz/(npts-1);
  double z0 = -npts/2*dz;
  std::random_device r;
  std::default_random_engine e1(r());
  std::uniform_int_distribution<int> phirange(-M_PI,M_PI);
  std::normal_distribution urand{0.0, ures};
  std::normal_distribution vrand{0.0, vres};

  CombineStereoPoints cp;
  for(unsigned ipt = 0; ipt < npts; ++ipt) {
    // set U direction
    double du = urand(e1);
    double dv = urand(e1);
    double eta = phirange(e1);
    VEC3 udir(cos(eta),sin(eta),0.0);
    VEC3 vdir(-sin(eta),cos(eta),0.0);
    double z = z0 + dz*ipt;
    VEC3 spos = pos + dir*z + du*udir + dv*vdir;
    StereoPoint spoint(spos,udir,ures*ures,vres*vres);
    cp.addPoint(spoint,ipt);
  }


  return 0;
}
