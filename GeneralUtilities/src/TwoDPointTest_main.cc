#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

using mu2e::TwoDPoint;
using mu2e::CombinedTwoDPoints;
using VEC2 = ROOT::Math::XYVectorF;
static struct option long_options[] = {
  {"p1x",     required_argument, 0, 'x' },
  {"p1y",     required_argument, 0, 'y'  },
  {"p2x",     required_argument, 0, 'X'  },
  {"p2y",     required_argument, 0, 'Y'  },
  {"p1ures",     required_argument, 0, 'u'  },
  {"p1vres",     required_argument, 0, 'v'  },
  {"p2ures",     required_argument, 0, 'U'  },
  {"p2vres",     required_argument, 0, 'V'  },
  {"p1cos",     required_argument, 0, 'c'  },
  {"p2cos",     required_argument, 0, 'C'  },
  {NULL, 0,0,0}
};
int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  float p1x(2.0), p1y(0.0), p2x(1.0), p2y(0.0), p1ures(1.0), p2ures(1.0);
  float p1vres(1.0), p2vres(1.0), p1cos(1.0), p2cos(1.0);
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'x' : p1x = atof(optarg);
                 break;
      case 'X' : p2x = atof(optarg);
                 break;
      case 'y' : p1y = atof(optarg);
                 break;
      case 'Y' : p2y = atof(optarg);
                 break;
      case 'u' : p1ures = atof(optarg);
                 break;
      case 'U' : p2ures = atof(optarg);
                 break;
      case 'v' : p1vres = atof(optarg);
                 break;
      case 'V' : p2vres = atof(optarg);
                 break;
      case 'c' : p1cos = atof(optarg);
                 break;
      case 'C' : p2cos = atof(optarg);
                 break;
      default: exit(EXIT_FAILURE);
    }
  }

  VEC2 p1(p1x,p1y);
  VEC2 p2(p2x,p2y);
  std::vector<TwoDPoint> points;
  points.push_back(TwoDPoint(p1,p1cos, p1ures*p1ures,p1vres*p1vres));
  points.push_back(TwoDPoint(p2,p2cos, p2ures*p2ures,p2vres*p2vres));
  CombinedTwoDPoints cp(points);
  cp.print(std::cout);
  CombinedTwoDPoints cp2(points.front());
  cp2.addPoint(points.back());
  cp2.print(std::cout);

  CombinedTwoDPoints cp3(points.front());
  double dchi2 = cp3.dChi2(points.back());
  std::cout << "DChisquared = " << dchi2 << std::endl;

  return 0;
}
