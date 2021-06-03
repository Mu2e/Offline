//
// Create a two dimensional histogram of the magnetic field
//
// Original author Michael MacKenzie (March 2020)
//

// Mu2e includes
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/GeomHandle.hh"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// ROOT includes
#include "TFile.h"
#include "TH1F.h"
#include "TH2.h"

// cetlib includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

namespace mu2e {

  class BFieldPlotter : public art::EDAnalyzer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> plane     {Name("plane"     ), Comment("Axis plane is defined by (x, y, or z)")};
      fhicl::Atom<double>      planeValue{Name("planeValue"), Comment("Value on the axis the plane intersects (mm)")};
      fhicl::Atom<double>      axisOneMin{Name("axisOneMin"), Comment("Axis one lower edge (bin centered) for plotting (mm) (plane = x/y/z --> axis one = y/x/x)")};
      fhicl::Atom<double>      axisOneMax{Name("axisOneMax"), Comment("Axis one upper edge (bin centered) for plotting (mm) (plane = x/y/z --> axis one = y/x/x)")};
      fhicl::Atom<double>      axisTwoMin{Name("axisTwoMin"), Comment("Axis two lower edge (bin centered) for plotting (mm) (plane = x/y/z --> axis two = z/z/y)")};
      fhicl::Atom<double>      axisTwoMax{Name("axisTwoMax"), Comment("Axis two upper edge (bin centered) for plotting (mm) (plane = x/y/z --> axis two = z/z/y)")};
      fhicl::Atom<double>      mapBinSize{Name("mapBinSize"), Comment("Map bin size (mm) (must be a divisor of both axis lengths)"), 10.};
    };
    typedef art::EDAnalyzer::Table<Config> Parameters;

    explicit BFieldPlotter(const Parameters& pset);
    virtual ~BFieldPlotter() { }

    void analyze(const art::Event& e);

    void fillHistogram(BFMap const *map, art::ServiceHandle<art::TFileService> tfs);

  private:

    std::string _plane     ; // x, y, or z map plane
    double      _planeValue; // value to set the plane at (mm)
    double      _axisOneMin; // sampling axis values (mm)
    double      _axisOneMax; // sampling axis values (mm)
    double      _axisTwoMin; // sampling axis values (mm)
    double      _axisTwoMax; // sampling axis values (mm)
    double      _mapBinSize; // histogram bin size (mm)

    std::map<std::string,TH2F*> _hMap; //histogram of the map
    
  };

  BFieldPlotter::BFieldPlotter(const Parameters& pset) :
    art::EDAnalyzer(pset)
    , _plane     (pset().plane     ())
    , _planeValue(pset().planeValue())
    , _axisOneMin(pset().axisOneMin())
    , _axisOneMax(pset().axisOneMax())
    , _axisTwoMin(pset().axisTwoMin())
    , _axisTwoMax(pset().axisTwoMax())
    , _mapBinSize(pset().mapBinSize())
  {
    if(_axisOneMin >= _axisOneMax || _axisTwoMin >= _axisTwoMax) {
      throw cet::exception("BADCONFIG") << "BField mapping plane ill defined: "
					<< _axisOneMin << " < axisOneValues < " << _axisOneMax << ", "
					<< _axisTwoMin << " < axisTwoValues < " << _axisTwoMax;
    }
    if(_mapBinSize <= 0.) {
      throw cet::exception("BADCONFIG") << "BField map binning should be >= 0 but given as " << _mapBinSize;
    }

    if(!(_plane == "x" || _plane == "y" || _plane == "z")) {
      throw cet::exception("BADCONFIG") << "BField map plane not recognized! Options are x, y, or z but given " 
					<< _plane.c_str();
    }

  } //end constructor

  void BFieldPlotter::analyze(const art::Event& event) {
    art::ServiceHandle<art::TFileService> tfs;
    // Get the magnetic field manager and the map lists
    GeomHandle<BFieldManager> bf;
    const BFieldManager::MapContainerType& innerMaps = bf->getInnerMaps();
    const BFieldManager::MapContainerType& outerMaps = bf->getOuterMaps();
 
    //plot inner maps
    for(std::shared_ptr<BFMap> const & map : innerMaps) { 
      fillHistogram(map.get(), tfs);
    } 
    
    //plot outer maps
    for(std::shared_ptr<BFMap> const & map : outerMaps) { 
      fillHistogram(map.get(), tfs);
    } 

    //plot the default field returned from the manager
    fillHistogram(NULL, tfs);
  } // end analyze


  //fill a histogram with the magnetic field
  void BFieldPlotter::fillHistogram(BFMap const *map, art::ServiceHandle<art::TFileService> tfs) {
    //define histogram binning such that map edges are bin centers and inside histogram
    long nbinsOne = std::lround((_axisOneMax - _axisOneMin)/_mapBinSize) + 1;
    long nbinsTwo = std::lround((_axisTwoMax - _axisTwoMin)/_mapBinSize) + 1;
    //check that the map bin step size works with the edges given
    if(std::abs(_axisOneMin + (nbinsOne-1)*_mapBinSize - _axisOneMax > (_axisOneMax-_axisOneMin)/(100.*(nbinsOne-1))))
      throw cet::exception("BADCONFIG") << "BField mapping axis values not steppable with step size given: "
					<< _axisOneMin << " < axisOneValues < " << _axisOneMax << ", "
					<< "step size = " << _mapBinSize;
	
    if(std::abs(_axisTwoMin + (nbinsTwo-1)*_mapBinSize - _axisTwoMax > (_axisTwoMax-_axisTwoMin)/(100.*(nbinsTwo-1))))
      throw cet::exception("BADCONFIG") << "BField mapping axis values not steppable with step size given: "
					<< _axisTwoMin << " < axisTwoValues < " << _axisTwoMax << ", "
					<< "step size = " << _mapBinSize;

    //if a null map, plot the default field from the manager
    const std::string name = (map) ? map->getKey() : "default";

    if(_hMap[name]) return; //already filled the maps in a previous event

    //make a directory for the map
    art::TFileDirectory tfdir = tfs->mkdir( ("BFieldMapper_"+name).c_str() );

    //define a histogram
    _hMap[name] = tfdir.make<TH2F>(("hMap"+name).c_str()  , (name + " Magnetic field map").c_str(), 
				   // add offsets so all values are bin centers and fit edge into the map
				   nbinsOne, _axisOneMin - _mapBinSize/2.,_axisOneMax + _mapBinSize/2.,
				   nbinsTwo, _axisTwoMin - _mapBinSize/2.,_axisTwoMax + _mapBinSize/2.);
    
    GeomHandle<BFieldManager> bf; //only needed in null map case
    //Loop through the points in the magnetic field
    double axisOne, axisTwo;
    for(int binOne = 0; binOne < nbinsOne; ++binOne) {
      axisOne = _axisOneMin + binOne*_mapBinSize;
      if(binOne == 0 || binOne == nbinsOne-1) //if first or last bin, at slight step into map to avoid edge issues
	axisOne += (binOne) ? -_mapBinSize/100. : _mapBinSize/100.;
      for(int binTwo = 0; binTwo < nbinsTwo; ++binTwo) {
	axisTwo = _axisTwoMin + binTwo*_mapBinSize;
	if(binTwo == 0 || binTwo == nbinsTwo-1) //if first or last bin, at slight step into map to avoid edge issues
	  axisTwo += (binTwo) ? -_mapBinSize/100. : _mapBinSize/100.;
	CLHEP::Hep3Vector point;
	CLHEP::Hep3Vector field;
	if(_plane == "x") 
	  point = CLHEP::Hep3Vector(_planeValue, axisOne, axisTwo);
	else if(_plane == "y") 
	  point = CLHEP::Hep3Vector(axisOne, _planeValue, axisTwo);
	else if(_plane == "z") 
	  point = CLHEP::Hep3Vector(axisOne, axisTwo, _planeValue);
	if(map)
	  map->getBFieldWithStatus(point,field);
	else //if null, get the manager and plot the field returned by it
	  bf->getBFieldWithStatus(point,field);

	_hMap[name]->Fill(axisOne, axisTwo, field.mag()); //fill with weight of the field magnitude
      } //end axis two loop
    } //end axis one loop
  } //end fillHistorgram


}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::BFieldPlotter);
