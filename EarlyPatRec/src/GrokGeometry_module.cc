//=============================================================================
//
// Plugin to look at the straws in the TTracker.
// This is a geometry understander so that pattern recognition can work
// with a convenient set of concepts, yet those will have the proper data.
//
// $Id: GrokGeometry_module.cc,v 1.14 2013/10/21 20:34:14 gandr Exp $
// $Author: gandr $
// $Date: 2013/10/21 20:34:14 $
//
// Original author: Mark Fischler
//
//=============================================================================
// C++ includes
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <cassert>

#include "CLHEP/Vector/TwoVector.h"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// Mu2e includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "DataProducts/inc/StrawId.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "TTrackerGeom/inc/Station.hh"

using namespace std;

namespace mu2e {;
  enum PrintLevel { quiet  =-1,
                    normal = 0,
                    verbose= 1};
 //
 // Abstracted essence of geometric information about one straw,
 // for purposes of the PRF pattern regognition algorithm F
 struct PRF_Straw{
 public:
   PRF_Straw() : index(-1) {}
   PRF_Straw(Straw const& str)
   : sid        (str.id())
   , deviceId   (sid.getDeviceId())
   , layer      (sid.getLayerId().getLayer())
   , sector     (sid.getLayerId().getSector())
   , halfLength (str.getHalfLength())
   , midpoint   (str.getMidPoint())
   , direction  (str.getDirection())
   , index      (str.index())
   , nextOuterL (str.nextOuterSameLayer().asInt())
   , nextInnerL (str.nextInnerSameLayer().asInt())
   , nextOuterP (str.nextOuterInPanel().asInt())
   , nextInnerP (str.nextInnerInPanel().asInt())
   { assert(std::abs(direction.z()) < 0.0001); }

   bool operator<(const PRF_Straw other) const {
   // Here so that PRF_Straw can be put into a standard container
      return (index.asInt() < other.index.asInt());
   }

  public:
   StrawId sid;
   int    deviceId;
   int    layer;
   int    sector;
   double halfLength;
   CLHEP::Hep3Vector midpoint;
   CLHEP::Hep3Vector direction;
   StrawIndex index;
   int nextOuterL;
   int nextInnerL;
   int nextOuterP;
   int nextInnerP;
   
   std::ostream & Print(std::ostream & os) const
    {
      os  << "Straw " << index << ": "
          << " Device Id " << deviceId << " Layer " << layer
          << " Sector " << sector  << "\n"
          << "half-length " << halfLength
          << " mid-point " << midpoint << "\n"
          << "direction  " << direction 
          << " radius " 
          << std::sqrt(midpoint.x()*midpoint.x()+midpoint.y()*midpoint.y()) 
          << "\n"
          << "nextOuterSameLayer: "   << nextOuterL 
          << "  nextInnerSameLayer: " << nextInnerL 
          << "\nnextOuterInPanel: " << nextOuterP 
          << "  nextInnerInPanel: " << nextInnerP << "\n";          
      return os;
   }
 }; // end of definition of struct PRF_Straw
  std::ostream & operator<<(std::ostream & os, PRF_Straw const & straw)
  { return straw.Print(os); }

  //--------------------------------------------------------------------
  //
  //
  class GrokGeometry : public art::EDAnalyzer {
  public:
    explicit GrokGeometry(fhicl::ParameterSet const& pset):
      art::EDAnalyzer(pset),
      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel"))
    {
    }
    virtual ~GrokGeometry() { }

    virtual void beginJob();
    virtual void beginRun(art::Run const &);

    void analyze( art::Event const& e);
  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;
    // Module label of the geerator module.
    std::string _generatorModuleLabel;
    // Module label of the g4 module that made the hits.
    std::string _g4ModuleLabel;

  }; // end of GrokGeometry class definition


  void GrokGeometry::beginJob(){
    cout << "Diaglevel: "
         << _diagLevel << " "
         << _maxFullPrint<<endl;

    art::ServiceHandle<art::TFileService> tfs;

  } // end of GrokGeometry beginJob

  void GrokGeometry::beginRun(art::Run const &){

    cout << "****************** beginRun ****************\n\n";
    const Tracker& tracker = getTrackerOrThrow();
    const std::deque<Straw>&  trackerStraws = tracker.getAllStraws();
    typedef std::deque<Straw>::const_iterator DecIt;
    for (DecIt s = (DecIt)trackerStraws.begin(); s != trackerStraws.end(); ++s)
    {
      Straw const& str           = *s;
      PRF_Straw straw ( str );
      cout << straw;
    }
    
    cout << "\n\n****************** Stations viewpoint ****************\n\n";

    GeomHandle<TTracker> ttracker_handle;
    const TTracker& ttracker = *ttracker_handle;
    std::vector<Station> stations = ttracker.getStations();
    typedef std::vector<Station>::const_iterator StationsIt;
    typedef std::vector<Face const *>::const_iterator FacePtrIt;
    StationsIt endStations = stations.end();
    for (StationsIt s = stations.begin(); s != endStations; ++s) {
      cout << "**** Station " << s->id() << "  z = " << s->midZ() << "\n";
      cout << "  Faces: \n";
      FacePtrIt endFaces = s->getFaces().end();
      for (FacePtrIt f = s->getFaces().begin(); f != endFaces; ++f) {
        cout << "    " << (*f)->id() << "\n";
      } 
      typedef std::vector<Plane>::const_iterator PlanesIt;
      PlanesIt endPlanes = s->getPlanes().end();
      for (PlanesIt p = s->getPlanes().begin(); p != endPlanes; ++p) {
        cout << "    Plane " << p->id() << " z = " << p->midZ() << "\n";
        typedef std::vector<Face>::const_iterator FacesIt;
        FacesIt endFaces = p->getFaces().end();
        for (FacesIt f = p->getFaces().begin(); f != endFaces; ++f) {
          cout << "      Face " << f->id() << " z = " << f->midZ() << "\n";
          typedef std::vector<Panel>::const_iterator PanelsIt;
          PanelsIt endPanels = f->getPanels().end();
          for (PanelsIt pn = f->getPanels().begin(); pn != endPanels; ++pn) {
            cout << "        Panel " << pn->id() 
                 << " z = " << pn->midZ() 
                 << " view = " <<  pn->view() 
                 << " [" << pn->phi()*180.0/M_PI << "]\n";
            typedef std::vector<ZLayer>::const_iterator ZLayersIt;
            ZLayersIt endZLayers = pn->getZLayers().end();
            for (ZLayersIt z = pn->getZLayers().begin(); z != endZLayers; ++z) {
              cout << "          ZLayer " << z->id() 
                   << " z = " << z->midZ() << "\n";
              PRF_Straw straw0 ( *(z->getStraws().front()) );
              PRF_Straw strawN ( *(z->getStraws().back())  );
              cout << straw0 << "\n";
              cout << strawN << "\n";
            } // end of ZLayers output loop
          } // end of panels output loop
        } // end of faces output loop
      } // end of planes output loop
    } // end of stations output loop     
  } // end of GrokGeometry beginRun

  void GrokGeometry::analyze(art::Event const& evt)
  {
    if ( _diagLevel > 2 ) cout << "GrokGeometry: analyze() begin"<<endl;

  } // end of ::analyze.

}

using mu2e::GrokGeometry;
DEFINE_ART_MODULE(GrokGeometry);
