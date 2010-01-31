#ifndef HOUGHTRANSFORM_HH
#define HOUGHTRANSFORM_HH
//
// $Id: HoughTransform.hh,v 1.1 2009/12/09 17:36:36 rhbob Exp $
// $Author: rhbob $ 
// $Date: 2009/12/09 17:36:36 $
//
// performs Hough Transform looking for circles in the L-Tracker,
// closely tied to HitCluster algorithms
//original author R. Bernstein
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <vector>


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
      HoughTransform():
      	goodHoughTracks(false), numberOfHoughTracks(0){};

      ~HoughTransform(){};

      virtual std::string name() const {return "HoughTransform";}

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
      typedef std::vector<Hep3Vector> clusterCenterVector;
      std::vector< Hep3Vector > clusterCenters;
      std::vector<int> clusterSize;
      houghCandidates houghCircles;

      //methods
      //      void foundHoughTracks(GeomHandle<LTracker>& ltracker,edm::Handle<StepPointMCCollection>& hits,
      void foundHoughTracks(GeomHandle<LTracker>& ltracker,StepPointMCCollection const* hits,
			    houghCandidates&);
      void solveForCircle(double& x1,double& y1, double& x2, double& y2, double& x3, double& y3,double& radius,
			  double& x0,double& y0, double& dca);

      //      int countHitNeighbours( Straw const& straw, edm::Handle<StepPointMCCollection>& hits );
      int countHitNeighbours( Straw const& straw, StepPointMCCollection const* hits );
      //and a numbering scheme for the returned vector of radius, center x, center y no one needs to know about

      Hep3Vector computeClusterXYZ(std::vector<mu2e::hitcluster::Candidate>& candClust);

     private:
      //these are index numbers of hit straws associated with the found 
      //Hough tracks, and then someone else
      //knows which straw goes with which index
      std::vector<int> strawHitIndices; 

      double _x0;
      double _y0;
      double _radius;


    };

  } //namespace HoughTransform
}   //namespace mu2e

#endif
