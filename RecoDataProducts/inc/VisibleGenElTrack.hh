#ifndef RecoDataProducts_VisibleGenElTrack_hh
#define RecoDataProducts_VisibleGenElTrack_hh

// C++ includes.
#include <iostream>
#include <memory>
#include <utility>
#include <map>
#include <vector>
#include <algorithm>

// Mu2e includes.
#include "RecoDataProducts/inc/GenElHitData.hh"
//#include "MCDataProducts/inc/SimParticle.hh"
#include "art/Persistency/Common/Ptr.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"
//
// Data of the Electrons tracks that came from the targets
//
// $Id: VisibleGenElTrack.hh,v 1.1 2011/06/23 21:52:04 tassiell Exp $
// $Author: tassiell $
// $Date: 2011/06/23 21:52:04 $
//
// Original author G. Tassielli
//

namespace mu2e {

typedef art::Ptr<SimParticle> SimParticlePtr;
typedef std::vector<mu2e::GenElHitData> GenElHitDataCollection;

  class VisibleGenElTrack {
  public:
          VisibleGenElTrack():_trkID(0),_isConversionEl(false),_hitTimeSorted(false),_nTurns(0){}
          /*VisibleGenElTrack(SimParticlePtr &trkPtr_ ):
                  _trkPtr(trkPtr_),*/
          VisibleGenElTrack( SimParticleCollection::key_type &trkID_, std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData_ ):
                  _trkID(trkID_),
                  _isConversionEl(false),
                  _hitTimeSorted(false),
                  _nTurns(0)
                  {
                    _genTrackData=genTrackData_;
                  }
          ~VisibleGenElTrack(){}

          friend std::ostream& operator<< ( std::ostream& ost,
                                    const VisibleGenElTrack& gtd ){
                  //ost << "Track id " << gtd._trkPtr.key()<<" vertex "<< gtd._trkPtr->startPosition()<< " Lorentz Vector "<< gtd._trkPtr->startMomentum()<<std::endl;
                  ost << "Track id " << gtd._trkID<<" vertex "<< gtd._genTrackData.first<< " Lorentz Vector "<< gtd._genTrackData.second<<std::endl;
                  ost << "N loops "<<gtd._nTurns<<std::endl;
                  for ( unsigned short j=0; j < gtd._nTurns; j++ ) {
                          ost << j+1 <<"-th loop point: "<<gtd._hitDataCollection[ gtd._turnFrstHit[j] ]._hitPoint << " momentum: "<< gtd._hitDataCollection[ gtd._turnFrstHit[j] ]._hitMomentum <<std::endl;
                  }
                  ost << "Track Hits:"<<std::endl;
                  for ( size_t i=0; i < gtd._hitDataCollection.size(); i++ ) {
                          ost << gtd._hitDataCollection[i] <<std::endl;
                  }
                  return ost;
          }

          void addHitData( StrawHitPtr iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ );/*{
                  GenElHitData hitData_
                  _hitDataCollection.push_back( hitData_ );
          }*/
          void sort();
          bool isNotPresent( size_t iHit_ );
          bool& isConversionEl( );
          void setConversionEl();
          SimParticleCollection::key_type& getTrkID();
          CLHEP::Hep3Vector getTrkVertex();
          CLHEP::HepLorentzVector getTrkLrntzVec();
          size_t getNumOfHit();
          GenElHitData& getHit(size_t ih);
          unsigned short& FindNTurns();
          unsigned short& getNumOfLoops();
          GenElHitData& getithLoopHit(unsigned short ih);

  protected:
          //SimParticlePtr _trkPtr;
          SimParticleCollection::key_type _trkID;
          std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> _genTrackData;
          bool _isConversionEl;

          GenElHitDataCollection _hitDataCollection;
          //std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData;
          //std::vector hitDataCollection;
          //std::vector::iterator hd_it;
          //typename std::vector::iterator hd_it;
          bool _hitTimeSorted;
          unsigned short _nTurns;
          //std::vector< std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > turnStartData;
          std::vector< size_t > _turnFrstHit;

  };

   void VisibleGenElTrack::addHitData( StrawHitPtr iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ ) {
          _hitDataCollection.push_back(
                          GenElHitData( iHit_, mcHitTime_, isOverlapped_, isFirst_, hitPoint_, hitMomentum_ )
                          );
  }

   void VisibleGenElTrack::sort(){
          std::sort( _hitDataCollection.begin(), _hitDataCollection.end() );
          _hitTimeSorted=true;
  }

   bool VisibleGenElTrack::isNotPresent( size_t iHit_ ){
          //typename std::vector::iterator hd_it = _hitDataCollection.begin();
          for ( GenElHitDataCollection::iterator hd_it = _hitDataCollection.begin(); hd_it!= _hitDataCollection.end(); ++hd_it ) {
                  if(hd_it->_iHit.key()==iHit_) return false;
          }
          return  true;
  }

   bool& VisibleGenElTrack::isConversionEl(){
           return _isConversionEl;
   }

   void VisibleGenElTrack::setConversionEl(){
           _isConversionEl=true;
   }

   SimParticleCollection::key_type& VisibleGenElTrack::getTrkID() {
          //return _trkPtr.key();
          return _trkID;
  }

   CLHEP::Hep3Vector VisibleGenElTrack::getTrkVertex() {
          //return _trkPtr->startPosition();
          return _genTrackData.first;
  }

   CLHEP::HepLorentzVector VisibleGenElTrack::getTrkLrntzVec() {
          //return _trkPtr->startMomentum();
          return _genTrackData.second;
  }

   size_t VisibleGenElTrack::getNumOfHit() {
          return _hitDataCollection.size();
  }

   GenElHitData& VisibleGenElTrack::getHit(size_t ih){
          if ( ih<_hitDataCollection.size() ) return _hitDataCollection[ih];
          else throw cet::exception("RANGE") << "Element "<<ih<<" out of size of the HitDataCollection "<<_hitDataCollection.size()<<std::endl;
  }

   unsigned short& VisibleGenElTrack::FindNTurns(){
          _turnFrstHit.clear();
          if ( _hitDataCollection.empty() ) {
                  std::cerr<<"Not point stored for track to try find loops"<<std::endl;
                  _nTurns=0;
                  return _nTurns;
          }
          else if ( _hitDataCollection.size() < 10 ) {
                  std::cerr<<"Not enough point stored for track to try find loops"<<std::endl;
                  _nTurns=1;
                  _turnFrstHit.push_back(0);
                  return _nTurns;
          }
          if ( !_hitTimeSorted ) sort();
          _nTurns=1;
          _turnFrstHit.push_back(0);
          for ( size_t ih=1 ; ih<_hitDataCollection.size(); ih++){
                  if ( (_hitDataCollection[ih]._mcHitTime - _hitDataCollection[ih-1]._mcHitTime) > 1.0 ) {
                          _nTurns++;
                          _turnFrstHit.push_back(ih);
                  }
          }
          return _nTurns;
  }

   unsigned short& VisibleGenElTrack::getNumOfLoops(){
          return _nTurns;
  }

   GenElHitData& VisibleGenElTrack::getithLoopHit(unsigned short ih){
          if (ih>=_nTurns) throw cet::exception("RANGE") << "Requested Loop "<<ih<<" -th not found for the track"<<std::endl;
          return _hitDataCollection[ _turnFrstHit[ih] ];
  }


}

#endif /* RecoDataProducts_VisibleGenElTrack_hh */
