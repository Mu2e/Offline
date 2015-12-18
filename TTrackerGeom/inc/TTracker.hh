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

#include <deque>
#include <vector>

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
    TTracker(){}  // TODO: insert proper initializer list, starting w/ base class

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    void fillPointers () const;

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

    typedef std::vector<Panel>::size_type stypeLayer;
    bool isLegal(const LayerId& lid ) const{
      return ( isLegal(lid.getPanelId()) &&
               lid.getLayer() > -1   &&
               std::vector<Layer>::size_type(lid.getLayer()) < getPanel(lid.getPanelId()).getLayers().size()
               );
    }

    bool isLegal(const StrawId& strid) const{
      return ( isLegal(strid.getLayerId()) &&
               strid.getStraw() > -1       &&
               std::vector<Straw>::size_type(strid.getStraw()) < getLayer(strid.getLayerId()).getStraws().size()
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

    const Layer& getLayer ( const LayerId& lid ) const{
      return _planes.at(lid.getPlane()).getLayer(lid);
    }

    const Straw& getStraw ( const StrawId& strid ) const{
      return _planes.at(strid.getPlane()).getStraw(strid);
    }

    const Straw& getStraw ( StrawIndex i ) const{
      return _allStraws.at(i.asInt());
    }

    int nStations() const{
      return _stations.size();
    }

    const std::vector<Station>& getStations() const{
      return _stations;
    }

    const Station& getStation ( StationId id) const{
      return _stations.at(id);
    }

    const std::deque<Straw>& getAllStraws() const {return _allStraws;}

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

    const std::vector<double>& getManifoldHalfLengths () const{
      return _manifoldHalfLengths;
    }

    TubsParams getPlaneEnvelopeParams() const{
      return _planeEnvelopeParams;
    }

    const TubsParams& getInnerTrackerEnvelopeParams() const{
      return _innerTrackerEnvelopeParams;
    }

    PlacedTubs mother() const{
      return _mother;
    }

    // presence info for each straw.
    bool strawExists(StrawIndex const index) const {
      return _strawExists[index.asInt()];
    }

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
    std::deque<Straw>  _allStraws;

    // Deprecated: part of the ancient MECO TTracker design.  A few vestiges not yet removed.
    std::vector<Manifold> _allManifolds;

    // Outer envelope that holds the new style support structure.
    PlacedTubs _mother;

    // The envelope that holds all of the planes in the tracker, including the plane supports.
    TubsParams _innerTrackerEnvelopeParams;

    // The envelope that holds all of the pieces in one plane, including supports.
    TubsParams _planeEnvelopeParams;

    // Which level of detail is present in the model of the support structure?
    SupportModel _supportModel;

    // All supports are the same shape; only relevant for _supportModel=="simple"
    Support _supportParams;

    // A more detailed model of the supports; again each plane has identical supports.
    // only relevant for _supportModel == "detailedv0".
    SupportStructure _supportStructure;

    // All manifolds are the same shape.
    // Deprecated: these will go away soon.
    std::vector<double> _manifoldHalfLengths;

    // Inner radius of inside edge of innermost straw.
    double _envelopeInnerRadius;

    // presence info for each straw.
    std::vector<bool> _strawExists;

  };

} //namespace mu2e

#endif /* TTrackerGeom_TTracker_hh */
