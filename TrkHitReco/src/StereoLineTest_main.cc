#include "Offline/TrkHitReco/inc/StereoLine.hh"
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
  {"pz",     required_argument, 0, 'Z' },
  {"px",     required_argument, 0, 'X'  },
  {"py",     required_argument, 0, 'Y'  },
  {"rx",     required_argument, 0, 'x'  },
  {"ry",     required_argument, 0, 'y'  },
  {"ures",     required_argument, 0, 'u'  },
  {"vres",     required_argument, 0, 'v'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: StereoLineTest  --npts --dz (z range) --px --py --pz (point) --rx --ry (slope) --ures --vres (transverse resolutions) \n");
}

int main(int argc, char** argv) {

  int opt;
  unsigned npts(5);
  VEC3 pos, dir(0.2,0.2,1.0);
  int long_index =0;
  double ures(30.0), vres(2.0), deltaz(150.0);
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
      case 'Z' : pos.SetZ(atof(optarg));
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
  std::cout << "Simulating " << npts << " points " << std::endl;
  std::cout << "Position " << pos << std::endl;
  std::cout << "Direction " << dir.Unit() << std::endl;
  std::cout << "U, V resolution " << ures  << " , " << vres << std::endl;
  std::cout << "Z range " << deltaz << std::endl;

  double dz = deltaz/(npts-1);
  std::random_device r;
  std::default_random_engine eng(r());
  std::uniform_int_distribution<int> phirange(-M_PI,M_PI);
  std::normal_distribution urand{0.0, ures};
  std::normal_distribution vrand{0.0, vres};

  CombineStereoPoints cp;
  for(unsigned ipt = 0; ipt < npts; ++ipt) {
    // set U direction
    double du = urand(eng);
    double dv = vrand(eng);
    double eta = phirange(eng);
    VEC3 udir(cos(eta),sin(eta),0.0);
    VEC3 vdir(-sin(eta),cos(eta),0.0);
    double z = dz*ipt - 0.5*deltaz;
    VEC3 spos = pos + dir*z + du*udir + dv*vdir;
    StereoPoint spoint(spos,udir,ures*ures,vres*vres);
    cp.addPoint(spoint,ipt);
    std::cout << "Creating point " << spos << std::endl;
  }
  std::cout << "CombineStereoPoints with " << cp.nPoints() << " points " << cp.point() << " chisq " << cp.chisquared() << " consistency " << cp.consistency() << std::endl;

  mu2e::StereoLine sline;
  bool doline = cp.stereoLine(sline);

  std::cout << "Stereo Line status " << doline << std::endl;
  if(doline) std::cout << sline;
  std::cout << std::endl;




  return 0;
}
