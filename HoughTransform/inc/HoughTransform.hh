#ifndef HOUGHTRANSFORM_HH
#define HOUGHTRANSFORM_HH
//
// $Id: HoughTransform.hh,v 1.6 2010/05/18 20:28:11 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 20:28:11 $
//
// helps perform Hough Transform looking for circles in the L-Tracker,
// closely tied to HitCluster algorithms.  
//original author R. Bernstein. Rework by P. Shanahan
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

// CLHEP inlcudes
#include "CLHEP/Vector/TwoVector.h"
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Services/interface/TFileService.h"
#include "FWCore/Framework/interface/TFileDirectory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "LTrackerGeom/inc/LTracker.hh"
#include "ToyDP/inc/StepPointMCCollection.hh"
#include "Mu2eUtilities/inc/sqrtOrThrow.hh"
#include "GeneralUtilities/inc/pow.hh"
#include "HitCluster/inc/HitCluster.hh"





namespace mu2e{

  namespace houghtransform{
    //---------------------------------------------------------
    //
    //
    class HoughTransform{

    public:
      typedef std::vector<mu2e::hitcluster::HitCluster> ClusterList;

      HoughTransform(ClusterList& clusters) :
       _hitClusters(clusters), goodHoughTracks(false), numberOfHoughTracks(0){
           FindCenters(); };

      ~HoughTransform(){};

      virtual std::string name() const {return "HoughTransform";}

      static void MakeClusters(StepPointMCCollection const* hits,
                                  ClusterList& clusters);


      struct houghCircleStruct
      {
        double radius;
        double x0;
        double y0;
        double dca;
        int numberOfStraws;
        houghCircleStruct(double radius_,double x0_, double y0_ ,double dca_, int numberOfStraws_):
          radius(radius_),x0(x0_),y0(y0_),dca(dca_),numberOfStraws(numberOfStraws_){}
      };

      typedef std::vector< houghCircleStruct > houghCandidates;

      //data
      bool goodHoughTracks;
      int numberOfHoughTracks;
      int minHitsforHoughFit;


      //methods
      //      void foundHoughTracks(GeomHandle<LTracker>& ltracker,edm::Handle<StepPointMCCollection>& hits,
      void foundHoughTracks(double radius, GeomHandle<LTracker>& ltracker,
                            houghCandidates&);
      void foundHoughTracks(GeomHandle<LTracker>& ltracker,
                            houghCandidates&);


      //      int countHitNeighbours( Straw const& straw, edm::Handle<StepPointMCCollection>& hits );
      int countHitNeighbours( Straw const& straw, StepPointMCCollection const* hits );
      //and a numbering scheme for the returned vector of radius, center x, center y no one needs to know about

      CLHEP::Hep3Vector computeClusterXYZ(std::vector<mu2e::hitcluster::Candidate>& candClust);

     private:

      static const bool _use3P=false;// use geometric cicle finder, not algebraic

      ClusterList& _hitClusters; // list of clusters input to Hough finding

      // regardless of where clusters come from, it's a good idea to 
      // make a list of centers so we don't keep recalculating them in
      // nested loops. 
      // we only use x and y, so it could probably be Hep2Vector, but
      // why tempt fate?
      typedef std::vector<CLHEP::Hep3Vector> clusterCenterVector;
      clusterCenterVector _clusterCenters;

      std::vector<int> _clusterSizes;
      houghCandidates _houghCircles;

      //these are index numbers of hit straws associated with the found 
      //Hough tracks, and then someone else
      //knows which straw goes with which index
      // PSH - not used? std::vector<int> strawHitIndices; 

      // PSH - not used? double _x0;
      // PSH - not used? double _y0;
      // PSH - not used? double _radius;

      bool solveForCircle(const CLHEP::Hep3Vector& v1,
                          const CLHEP::Hep3Vector& v2,
                          const CLHEP::Hep3Vector& v3,
                          double& radius, double& x0,double& y0, double& dca);

      bool solveForCircle3P(const CLHEP::Hep3Vector& v1, 
                            const CLHEP::Hep3Vector& v2, 
                            const CLHEP::Hep3Vector& v3, 
                            double& radius, double& x0,double& y0, double& dca);

      bool solveForCircle2P(const CLHEP::Hep3Vector& v1,
                            const CLHEP::Hep3Vector& v2, double radius,
                            double& x0m, double& y0m, double& dcam,
                            double& x0p, double& y0p, double& dcap);

      void FindCenters(); // find centers of all the clusters and store them


    };

  } //namespace HoughTransform
}   //namespace mu2e

#endif
