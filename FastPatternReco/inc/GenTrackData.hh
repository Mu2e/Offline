#ifndef GENTRACKDATA_HH
#define GENTRACKDATA_HH
//
// this is a old version of Data of the Electrons tracks that came from the targets
//
// $Id: GenTrackData.hh,v 1.1 2011/06/23 21:56:11 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/23 21:56:11 $
//
// Original author G. Tassielli
//

// C++ includes.
#include <iostream>
#include <memory>
#include <utility>
#include <map>
#include <vector>
#include <algorithm>

// Mu2e includes.
#include "MCDataProducts/inc/SimParticleCollection.hh"
//#include "FastPatternReco/inc/HitPerTrackData.hh"
//#include "FastPatternReco/inc/ITHitPerTrackData.hh"
//#include "FastPatternReco/inc/TTHitPerTrackData.hh"

namespace mu2e {

  template <typename HITTD> class GenTrackData {
  public:
          GenTrackData():trkID(0),hitTimeSorted(false),nTurns(0){}
          GenTrackData(SimParticleCollection::key_type &trkID_, std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData_):
                  trkID(trkID_),
                  hitTimeSorted(false),
                  nTurns(0)
                  {
                    genTrackData=genTrackData_;
                  }
          ~GenTrackData(){}

          friend std::ostream& operator<< ( std::ostream& ost,
                                    const GenTrackData<HITTD>& gtd ){

                  ost << "Track id " << gtd.trkID<<" vertex "<< gtd.genTrackData.first<< " Lorentz Vector "<< gtd.genTrackData.second<<std::endl;
                  ost << "N loops "<<gtd.nTurns<<std::endl;
                  for ( unsigned int j=0; j < gtd.nTurns; j++ ) {
                          ost << j+1 <<"-th loop point: "<<gtd.hitDataCollection[ gtd.turnFrstHit[j] ].hitPoint << " momentum: "<< gtd.hitDataCollection[ gtd.turnFrstHit[j] ].hitMomentum <<std::endl;
                  }
                  ost << "Track Hits:"<<std::endl;
                  for ( unsigned int i=0; i < gtd.hitDataCollection.size(); i++ ) {
                          ost << gtd.hitDataCollection[i] <<std::endl;
                  }
                  return ost;
          }

          void addHitData( HITTD hitData_ );/*{
                  hitDataCollection.push_back( hitData_ );
          }*/
          void sort();
          bool isNotPresent( size_t iHit_ );
          SimParticleCollection::key_type& getTrkID();
          CLHEP::Hep3Vector getTrkVertex();
          CLHEP::HepLorentzVector getTrkLrntzVec();
          unsigned int getNumOfHit();
          HITTD getHit(unsigned int ih);
          unsigned int& FindNTurns();
          unsigned int& getNumOfLoops();
          HITTD getithLoopHit(unsigned int ih);

  protected:
          SimParticleCollection::key_type trkID;
          std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData;
          std::vector<HITTD> hitDataCollection;
          //std::vector<HITTD>::iterator hd_it;
          typename std::vector<HITTD>::iterator hd_it;
          bool hitTimeSorted;
          unsigned int nTurns;
          //std::vector< std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > turnStartData;
          std::vector< unsigned int > turnFrstHit;

  };

  template <typename HITTD> void GenTrackData<HITTD>::addHitData( HITTD hitData_ ) {
          hitDataCollection.push_back( hitData_ );
  }

  template <typename HITTD> void GenTrackData<HITTD>::sort(){
          std::sort( hitDataCollection.begin(), hitDataCollection.end() );
          hitTimeSorted=true;
  }

  template <typename HITTD> bool GenTrackData< HITTD >::isNotPresent( size_t iHit_ ){
          //typename std::vector<HITTD>::iterator hd_it = hitDataCollection.begin();
          for ( /*std::vector<HITTD>::iterator*/ hd_it = hitDataCollection.begin(); hd_it!= hitDataCollection.end(); ++hd_it ) {
                  if(hd_it->iHit==iHit_) return false;
          }
          return  true;
  }

  template <typename HITTD> SimParticleCollection::key_type& GenTrackData<HITTD>::getTrkID() {
          return trkID;
  }

  template <typename HITTD> CLHEP::Hep3Vector GenTrackData<HITTD>::getTrkVertex() {
          return genTrackData.first;
  }

  template <typename HITTD> CLHEP::HepLorentzVector GenTrackData<HITTD>::getTrkLrntzVec() {
          return genTrackData.second;
  }

  template <typename HITTD> unsigned int GenTrackData<HITTD>::getNumOfHit() {
          return hitDataCollection.size();
  }

  template <typename HITTD> HITTD GenTrackData<HITTD>::getHit(unsigned int ih){
          if ( ih<hitDataCollection.size() ) return hitDataCollection[ih];
          else throw cet::exception("RANGE") << "Element "<<ih<<" out of size of the HitDataCollection "<<hitDataCollection.size()<<std::endl;
  }

  template <typename HITTD> unsigned int& GenTrackData<HITTD>::FindNTurns(){
          turnFrstHit.clear();
          if ( hitDataCollection.empty() ) {
                  std::cerr<<"Not point stored for track to try find loops"<<std::endl;
                  nTurns=0;
                  return nTurns;
          }
          else if ( hitDataCollection.size() < 10 ) {
                  std::cerr<<"Not enough point stored for track to try find loops"<<std::endl;
                  nTurns=1;
                  turnFrstHit.push_back(0);
                  return nTurns;
          }
          if ( !hitTimeSorted ) sort();
          nTurns=1;
          turnFrstHit.push_back(0);
          for ( unsigned int ih=1 ; ih<hitDataCollection.size(); ih++){
                  if ( (hitDataCollection[ih].mcHitTime - hitDataCollection[ih-1].mcHitTime) > 1.0 ) {
                          nTurns++;
                          turnFrstHit.push_back(ih);
                  }
          }
          return nTurns;
  }

  template <typename HITTD> unsigned int& GenTrackData<HITTD>::getNumOfLoops(){
          return nTurns;
  }

  template <typename HITTD> HITTD GenTrackData<HITTD>::getithLoopHit(unsigned int ih){
          if (ih>=nTurns) throw cet::exception("RANGE") << "Requested Loop "<<ih<<" -th not found for the track"<<std::endl;
          return hitDataCollection[ turnFrstHit[ih] ];
  }


}

#endif /* GENTRACKDATA_HH */
