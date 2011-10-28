#ifndef HitCluster_HitCluster_hh
#define HitCluster_HitCluster_hh
//forms clusters of adjacent straws in the L-Tracker for pattern recognition
//
// $Id: HitCluster.hh,v 1.11 2011/10/28 18:47:06 greenc Exp $
// $Author: greenc $
// $Date: 2011/10/28 18:47:06 $
//
//original author R. Bernstein
//

// C++ includes.
#include <cmath>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

// Framework includes.
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TNtuple.h"

// Mu2e includes.
#include "LTrackerGeom/inc/LTracker.hh"
#include "Mu2eUtilities/inc/sqrtOrThrow.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"


namespace mu2e{

  namespace hitcluster{
    //---------------------------------------------------------
    //
    //


    struct Candidate{

      //Identifier of the hit straw ( StrawIndex or VolumeId ).
      int id;

      //Index in the original hit container.
      StepPointMC const* hitPointer;

      explicit Candidate( int id_=-1, StepPointMC const* hitPointer_= 0):
        id(id_),
        hitPointer(hitPointer_){
      }

      //Compiler generated versions are OK for:
      //copy c'tor, destructor, operator=

      //The required behaviour of the comparison operators is that they
      //compare only the values of the id belonging to the seed hit of this candidate.
      bool operator==( Candidate const& rhs) const{
        return ( hitPointer->volumeId() == rhs.hitPointer->volumeId()  );
      }

      bool operator!=( Candidate const& rhs) const{
        return !( *this == rhs );
      }

      bool operator<( Candidate const& rhs) const{
        return ( hitPointer->volumeId() < rhs.hitPointer->volumeId() );
      }
      bool operator>( Candidate const& rhs) const{
        return ( hitPointer->volumeId() > rhs.hitPointer->volumeId());
      }

      bool operator>=( Candidate const& rhs) const{
        return ( hitPointer->volumeId() >= rhs.hitPointer->volumeId());
      }

      bool operator<=( Candidate const& rhs) const{
        return ( hitPointer->volumeId() <= rhs.hitPointer->volumeId());
      }

    };

    inline std::ostream& operator<<( std::ostream& ost,
                                     Candidate const& c){
      ost << "( "
          << c.id << ","
          << c.hitPointer->position() <<")";
      return ost;
    }

    class HitCluster{

    public:

      HitCluster(const StepPointMC& hit,Straw const& straw,StepPointMCCollection const* hits)
        :_hit(&hit),
         _straw(&straw),
         _hits(hits)
      {
        listOfStraws = HitCluster::findHitNeighbours();
      }

      virtual ~HitCluster();

      //data
      std::vector<int> strawsInThisCluster;

      //methods

      virtual std::string name() const {return "HitCluster";}

      //this part written by Kutschke
      typedef std::vector<Candidate> hitNeighbours;
      hitNeighbours neighbouringStraws;

      hitNeighbours getStraws();


      void matchAndMerge    (bool& match, std::vector<HitCluster>& finalClusters);
      void cleanUpDuplicates();


    private:

      hitNeighbours listOfStraws;

      const StepPointMC* _hit;
      const Straw* _straw;
      StepPointMCCollection const* _hits;
      CLHEP::Hep3Vector clusterXYZ;

      hitNeighbours findHitNeighbours();
      void       addStraw(Candidate addThisStraw);
      int        getStraw(int ithStraw);
    };

  } //namespace HitCluster
}   //namespace mu2e

#endif /* HitCluster_HitCluster_hh */
