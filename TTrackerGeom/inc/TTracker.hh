#ifndef TTrackerGeom_TTracker_hh
#define TTrackerGeom_TTracker_hh
//
// Hold all geometry and identifier information about
// a TTracker.  This is intended as a "data only"
// class.
//
// $Id: TTracker.hh,v 1.15 2014/04/11 04:39:44 genser Exp $
// $Author: genser $
// $Date: 2014/04/11 04:39:44 $
//
// Original author Rob Kutschke
//

//#include <deque>
#include <vector>
#include <array>
#include <limits>
//#include <iomanip>
//#include <iostream>

#include "cetlib_except/exception.h"

//#include "DataProducts/inc/StrawId2.hh"  included via Straw via Tracker (-> Straw)

#include "TTrackerGeom/inc/Manifold.hh"
#include "TTrackerGeom/inc/Support.hh"
#include "TTrackerGeom/inc/SupportModel.hh"
#include "TTrackerGeom/inc/SupportStructure.hh"

#include "TrackerGeom/inc/Plane.hh"
#include "DataProducts/inc/PanelId.hh"
#include "TrackerGeom/inc/StrawDetail.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "GeomPrimitives/inc/TubsParams.hh"
#include "GeomPrimitives/inc/PlacedTubs.hh"

#include "TTrackerGeom/inc/Station.hh"

namespace mu2e {

  class TTracker: public Tracker{

    friend class TTrackerMaker;

  public:

    // =============== NewTracker Public Objects Start ==============

    // constexpr static int _nttstraws = StrawId2::_nplanes *
    //                                   StrawId2::_npanels *
    //                                   StrawId2::_nstraws; // uncomment after eliminating Tracker

    // constexpr static int _maxRedirect = std::numeric_limits<uint16_t>::max();
    constexpr static uint16_t _maxRedirect =
      ((StrawId2::_nplanes -1) << StrawId2::_planesft) +
      ((StrawId2::_npanels -1) << StrawId2::_panelsft) +
      StrawId2::_nstraws;

    // =============== NewTracker Public Objects End   ==============

    TTracker(){
      if (StrawId2::_nlayers != 2)
        throw cet::exception("GEOM")
          << "Expect configuration with 2 layers per panel\n";
    }  // TODO: insert proper initializer list, starting w/ base class

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    void fillPointers () const;
    void fillPointers2 () const;

    double rOut() const { return _rOut;}
    double z0()   const { return _z0;}
    double zHalfLength() const;

    double strawRadius() const{
      return getStraw(StrawId(0,0,0,0)).getDetail().outerRadius();
    }

    std::string const& envelopeMaterial() const { return _envelopeMaterial; }

    // Check for legal identifiers.
    bool isLegal(PlaneId d) const{
      return ( d>-1 &&
               std::vector<Plane>::size_type(d) <_planes.size()
               );
    };

    bool isLegal(const PanelId& pnlid) const{
      return (isLegal(pnlid.getPlaneId()) &&
              pnlid.getPanel() >-1   &&
              std::vector<Panel>::size_type(pnlid.getPanel()) < getPlane(pnlid.getPlaneId()).getPanels().size()
              );
    }

    // typedef std::vector<Panel>::size_type stypeLayer;
    // bool isLegal(const LayerId& lid ) const{
    //   return ( isLegal(lid.getPanelId()) &&
    //            lid.getLayer() > -1   &&
    //            std::vector<Layer>::size_type(lid.getLayer()) < getPanel(lid.getPanelId()).getLayers().size()
    //            );
    // }

    // bool isLegal(const StrawId& strid) const{
    //   return ( isLegal(strid.getLayerId()) &&
    //            strid.getStraw() > -1       &&
    //            strid.getStraw() < getLayer(strid.getLayerId()).nStraws()
    //            );
    // }

    bool isLegal(const StrawId& strid) const{
      return ( isLegal(strid.getPanelId()) &&
               strid.getStraw() > -1       &&
               strid.getStraw() < getPanel(strid.getPanelId()).nStraws()
               );
    }

    // Accessors
    int nPlanes() const{
      return _planes.size();
    }

    const std::vector<Plane>& getPlanes() const{
      return _planes;
    }

    const Plane& getPlane ( PlaneId id) const{
      return _planes.at(id);
    }

    const Panel& getPanel ( const PanelId& pnlid ) const{
      return _planes.at(pnlid.getPlane()).getPanel(pnlid);
    }

    const std::vector<double>& getManifoldHalfLengths () const{
      return _manifoldHalfLengths;
    }

    double panelOffset() const { return _panelZOffset; }

    // const Layer& getLayer ( const LayerId& lid ) const{
    //   return _planes.at(lid.getPlane()).getLayer(lid);
    // }

    const Straw& getStraw ( const StrawId& strid ) const{
      return _planes.at(strid.getPlane()).getStraw(strid);
    }

    // const Straw& getStraw ( StrawIndex i ) const{
    //   return _allStraws.at(i.asInt());
    // }

    int nStations() const{
      return _stations.size();
    }

    const std::vector<Station>& getStations() const{
      return _stations;
    }

    const Station& getStation ( StationId id) const{
      return _stations.at(id);
    }

    //    const std::deque<Straw>& getAllStraws() const {return _allStraws;}
    const std::array<Straw,TTracker::_nttstraws>& getAllStraws() const 
    {return _allStraws2;}

    const std::vector<StrawDetail>& getStrawDetails() const{
      return _strawDetails;
    }

    SupportModel getSupportModel() const{
      return _supportModel;
    }

    const Support& getSupportParams () const{
      return _supportParams;
    }

    const SupportStructure& getSupportStructure() const{
      return _supportStructure;
    }


    TubsParams getPlaneEnvelopeParams() const{
      return _planeEnvelopeParams;
    }

    TubsParams getPanelEnvelopeParams() const{
      return _panelEnvelopeParams;
    }

    const TubsParams& getInnerTrackerEnvelopeParams() const{
      return _innerTrackerEnvelopeParams;
    }

    PlacedTubs mother() const{
      return _mother;
    }

    // presence info for each straw.
    bool strawExists(StrawIndex const index) const {
      //      return _strawExists[index.asInt()];
      return _allStraws2_p.at(((_allStraws2.at(index.asInt())).id2()).asUint16()) != nullptr;
    }

    // =============== NewTracker Accessors Start ==============

    Straw const& getStraw( StrawId2 const id) const{
      return *(_allStraws2_p.at(id.asUint16()));
    }

    Plane const& getPlane( StrawId2 id ) const{
      return _allPlanes2.at(id.getPlane());
    }

    bool strawExists( StrawId2 const id) const{
      return _allStraws2_p.at(id.asUint16()) != nullptr;
    }

    const Straw& getStraw ( StrawIndex i ) const{
      // shold be correct by construction
      return _allStraws2.at(i.asInt());
    }

    const Straw& getStraw3 ( StrawIndex i ) const{
      // recalculation from the old to new straw layout
      uint16_t seqPanelNumber = i.asInt()/StrawId2::_nstraws;
      uint16_t panelNumber = seqPanelNumber%StrawId2::_npanels;
      uint16_t planeNumber = seqPanelNumber/StrawId2::_npanels;
      uint16_t panelNumberShifted = panelNumber << StrawId2::_panelsft;
      uint16_t planeNumberShifted = planeNumber << StrawId2::_planesft;
      uint16_t strawNumberInPanel = i.asInt()%StrawId2::_nstraws;
      constexpr static uint16_t strawsPerLayer =
        StrawId2::_nstraws/StrawId2::_nlayers;
      uint16_t sn = (strawNumberInPanel<strawsPerLayer) ?
        (strawNumberInPanel << 1 ) :
        (( strawNumberInPanel - strawsPerLayer) << 1 ) + 1;
      uint16_t i2 = planeNumberShifted + panelNumberShifted + sn;
      // std::cout << __func__ << " i, sn, i2 "
      //           << std::setw(6) << i
      //           << std::setw(6) << sn
      //           << std::setw(6) << i2
      //           << std::endl;
      return *(_allStraws2_p.at(i2));
    }

    const StrawId2 getStrawId2 ( StrawIndex i ) const{
      return (_allStraws2.at(i.asInt())).id2();
    }

    const StrawIndex getStrawIndex (  const StrawId2& id ) const{
      return (_allStraws2_p.at(id.asUint16()))-> index();
    }

    // =============== NewTracker Accessors End   ==============

#ifndef __CINT__

    // Loop over all straws and call F.
    // F can be a class with an operator() or a free function.
    template <class F>
    inline void forAllStraws ( F& f) const{
      for ( std::vector<Plane>::const_iterator i=_planes.begin(), e=_planes.end();
            i !=e; ++i){
        i->forAllStraws(f);
      }
    }

    template <class F>
    inline void forAllLayers ( F& f) const{
      for ( std::vector<Plane>::const_iterator i=_planes.begin(), e=_planes.end();
            i !=e; ++i){
        i->forAllLayers(f);
      }
    }

    template <class F>
    inline void forAllPanels ( F& f) const{
      for ( std::vector<Plane>::const_iterator i=_planes.begin(), e=_planes.end();
            i !=e; ++i){
        i->forAllPanels(f);
      }
    }

    template <class F>
    inline void forAllPlanes ( F& f) const{
      for ( std::vector<Plane>::const_iterator i=_planes.begin(), e=_planes.end();
            i !=e; ++i){
        f(*i);
      }
    }

#endif


  protected:

    // Position of the center of the tracker, in the Mu2e coordinate system.
    double _z0;

    // Outer radius of a logical volume that will just contain the entire tracker.
    double _rOut;

    // All envelope volumes are made of this.
    std::string _envelopeMaterial;

    // Detailed info about each type of straw.
    std::vector<StrawDetail> _strawDetails;

    // An TTracker is made of two planes, sides and vanes.
    std::vector<Plane> _planes;

    // An alternative viewpoint:
    // A TTracker is made of a collection of Stations.
    std::vector<Station> _stations;

    // There will be pointers to the objects in this container.
    //    std::deque<Straw>  _allStraws;

    // Deprecated: part of the ancient MECO TTracker design.
    // A few vestiges not yet removed.
    std::vector<Manifold> _allManifolds;

    // Outer envelope that holds the new style support structure.
    PlacedTubs _mother;

    // The envelope that holds all of the planes in the tracker,
    // including the plane supports.
    TubsParams _innerTrackerEnvelopeParams;

    // The envelope that holds all of the pieces in one plane, including supports.
    TubsParams _planeEnvelopeParams;

    // Ditto for Panel
    TubsParams _panelEnvelopeParams;

    // Which level of detail is present in the model of the support structure?
    SupportModel _supportModel;

    // All supports are the same shape; only relevant for _supportModel=="simple"
    Support _supportParams;
    double _panelZOffset; // introduced for version 5

    // A more detailed model of the supports; again each plane has identical supports.
    // only relevant for _supportModel == "detailedv0".
    SupportStructure _supportStructure;

    // All manifolds are the same shape.
    // Deprecated: these will go away soon.
    std::vector<double> _manifoldHalfLengths;

    // Inner radius of inside edge of innermost straw.
    double _envelopeInnerRadius;

    // presence info for each straw.
    //    std::vector<bool> _strawExists;

    // =============== NewTracker Private Objects Start ==============

    // Dense array.
    std::array<Plane,StrawId2::_nplanes> _allPlanes2;

    // Dense array.
    std::array<Straw,TTracker::_nttstraws> _allStraws2;

    // Sparse array: designed for indexing by StrawId2.
    // For all legal entries in StrawId2, this points to a straw in _straws2;
    // All other entries are null.
    std::array<Straw const*,TTracker::_maxRedirect> _allStraws2_p;

    // =============== NewTracker Private Objects End ==============

  };

} //namespace mu2e

#endif /* TTrackerGeom_TTracker_hh */
