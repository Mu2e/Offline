#include "Offline/GeneralUtilities/inc/TwoDPoint.hh"
#include "Offline/GeneralUtilities/inc/CombineTwoDPoints.hh"
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

using mu2e::TwoDPoint;
using mu2e::CombineTwoDPoints;
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
  {"inres",     required_argument, 0, 'i'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: TwoDPointTest  --p1x --p1y (point1 (x,y) --p1ures --p1vres (point1 u,v resolution) --p1cos (point1 u_x) (same for point2) --inres (circular intrinsic resolution) ");
}

int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  float p1x(2.0), p1y(0.0), p2x(1.0), p2y(0.0), p1ures(1.0), p2ures(1.0), inres(0.0);
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
      case 'i' : inres = atof(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }

  VEC2 p1(p1x,p1y);
  VEC2 p2(p2x,p2y);
  float invar = inres*inres;
  float p1uvar = p1ures*p1ures;
  float p2uvar = p2ures*p2ures;
  float p1vvar = p1vres*p1vres;
  float p2vvar = p2vres*p2vres;
  std::vector<TwoDPoint> points;
  points.push_back(TwoDPoint(p1,p1cos, p1uvar, p1vvar));
  points.push_back(TwoDPoint(p2,p2cos, p2uvar, p2vvar));
  CombineTwoDPoints cp(points,invar);
  cp.print(std::cout);
  std::cout << "DChi2 point1 " << cp.dChi2(0) << " DChi2 point2 " << cp.dChi2(1) << std::endl;

  CombineTwoDPoints cp2(points.front(),invar);
  cp2.addPoint(points.back());
  cp2.print(std::cout);

  CombineTwoDPoints cp3(points.front(),invar);
  double dchi2 = cp3.dChi2(points.back());
  std::cout << "DChisquared = " << dchi2 << std::endl;
  // outlier search:e
  //
  VEC2 p3(2*p2x,2*p2y);
  cp.addPoint(TwoDPoint(p3,p2cos,p2uvar,p2vvar));
  std::cout << "After Adding point" << std::endl;
  cp.print(std::cout);
  auto const& wts = cp.weights();
  double maxdchi2(-1.0);
  size_t maxkey(0);
  for(auto iwt = wts.begin(); iwt != wts.end(); ++iwt){
    double dchi2 = cp.dChi2(iwt->first);
    std::cout << "Key " << iwt->first << " Weight " << iwt->second.wt_ << " dchisq " << dchi2 << std::endl;
    if(dchi2 > maxdchi2){
      maxkey = iwt->first;
      maxdchi2 = dchi2;
    }
  }
  cp.removePoint(maxkey);
  std::cout << "After removing point " << maxkey << std::endl;
  cp.print(std::cout);


  return 0;
}
