#include "Offline/RecoGeom/inc/Cylinder.hh"
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
using mu2e::RecoGeom::IntersectFlag;
static struct option long_options[] = {
  {"x",     required_argument, 0, 'x' },
  {"y",     required_argument, 0, 'y'  },
  {"z",     required_argument, 0, 'z'  },
  {"costheta",     required_argument, 0, 't'  },
  {"phi",        required_argument, 0, 'p'  },
  {"tol",     required_argument, 0, 'T'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: RecoGeomTest  --x --y --z (point) --costheta, --phi (direction angles) --tol (tolerance)");
}

int main(int argc, char** argv) {

  int opt;
  int long_index =0;
  VEC3 point(0.0,0.0,0.0);
  double cost(0.5), phi(0.0), tol(1e-8);
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'x' : point.SetX(atof(optarg));
                 break;
      case 'y' : point.SetY(atof(optarg));
                 break;
      case 'z' : point.SetZ(atof(optarg));
                 break;
      case 't' : cost = atof(optarg);
                 break;
      case 'p' : phi = atof(optarg);
                 break;
      case 'T' : tol = atof(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }

  double sint = sqrt(1.0-cost*cost);
  VEC3 rdir(sint*cos(phi), sint*sin(phi), cost);
  Ray ray(rdir,point);

  VEC3 zdir(0.0,0.0,1.0);
  VEC3 origin(0.0,0.0,0.0);
  Annulus ann(zdir,origin,1.0,2.0);
  VEC3 udir(1.0,0.0,0.0);
  Rectangle rect(zdir,origin,udir,1.0,2.0);
  Cylinder cyl(zdir,origin,2.0,10.0);

  std::cout << "Test " << ray << std::endl;
  std::cout << "Test " << ann << std::endl;
  std::cout << "Test " << rect << std::endl;
  std::cout << "Test " << cyl << std::endl;

  if(ann.onSurface(point))
    std::cout << "On Annulus "<< std::endl;
  else
    std::cout <<  "Not On Annulus "<< std::endl;
  if(rect.onSurface(point))
    std::cout << "On Rectangle "<< std::endl;
  else
    std::cout << "Not On Rectangle "<< std::endl;
  if(cyl.onSurface(point))
    std::cout << "On Cylinder "<< std::endl;
  else
    std::cout << "Not On Cylinder "<< std::endl;

  double dist;
  auto iflag = ann.intersect(ray,dist,tol);
  if(iflag.onsurface_)
    std::cout << "Annulus intersect " << iflag << " at distance " << dist << " point " << ray.position(dist) << std::endl;
  else
    std::cout << "No Annulus intersection" << std::endl;

  iflag = rect.intersect(ray,dist,tol);
  if(iflag.onsurface_)
    std::cout << "Rectangle intersect " << iflag << " at distance " << dist << " point " << ray.position(dist) << std::endl;
  else
    std::cout << "No Rectangle intersection" << std::endl;

  iflag = cyl.intersect(ray,dist,tol);
  if(iflag.onsurface_)
    std::cout << "Cylinder intersect " << iflag << " at distance " << dist << " point " << ray.position(dist) << std::endl;
  else
    std::cout << "No Cylinder intersection" << std::endl;

  return 0;
}
