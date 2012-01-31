#ifndef MCDataProducts_VisibleGenElTrack_hh
#define MCDataProducts_VisibleGenElTrack_hh
//
// Data of the Electrons tracks that came from the targets
//
// $Id: VisibleGenElTrack.hh,v 1.2 2012/01/31 00:38:10 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/01/31 00:38:10 $
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
#include "MCDataProducts/inc/GenElHitData.hh"
//#include "MCDataProducts/inc/SimParticle.hh"
#include "art/Persistency/Common/Ptr.h"
#include "MCDataProducts/inc/SimParticleCollection.hh"

//using namespace art;

namespace mu2e {

typedef art::Ptr<SimParticle> SimParticlePtr;
//typedef std::vector<mu2e::GenElHitData> GenElHitDataCollection;
typedef std::map< art::Ptr<SimParticle>::key_type, mu2e::GenElHitData> GenElHitDataCollection;
typedef std::map<unsigned long, art::Ptr<SimParticle>::key_type> GenElHitDCll_OrderByTime;

  class VisibleGenElTrack {
  public:
          VisibleGenElTrack():_trkKey(0),_trkID(0),_isConversionEl(false),_hitTimeSorted(false),_nTurns(0){
                  _hitDataCollection.clear();
                  _orederedHitDataColl.clear();
          }
          /*VisibleGenElTrack(SimParticlePtr &trkPtr_ ):
                  _trkPtr(trkPtr_),*/
          VisibleGenElTrack( art::Ptr<SimParticle>::key_type &trkKey_, SimParticleCollection::key_type &trkID_, std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData_ ):
                  _trkKey(trkKey_),
                  _trkID(trkID_),
                  _isConversionEl(false),
                  _hitTimeSorted(false),
                  _nTurns(0)
                  {
                    //_trkPtr=const_cast<SimParticlePtr&>(trkPtr_);
                    //_trkID=_trkPtr->id();
                    _genTrackData=genTrackData_;
                    _hitDataCollection.clear();
                    _orederedHitDataColl.clear();
                  }
          ~VisibleGenElTrack(){}

          friend std::ostream& operator<< ( std::ostream& ost,
                                    const VisibleGenElTrack& gtd_ ){
                  VisibleGenElTrack &gtd = const_cast<VisibleGenElTrack &>(gtd_);
                  //ost << "Track id " << gtd._trkPtr.key()<<" vertex "<< gtd._trkPtr->startPosition()<< " Lorentz Vector "<< gtd._trkPtr->startMomentum()<<std::endl;
                  ost << "Track id " << gtd._trkID<<" trak key "<<gtd._trkKey<<" vertex "<< gtd._genTrackData.first<< " Lorentz Vector "<< gtd._genTrackData.second<<std::endl;
                  ost << "N loops "<<gtd._nTurns<<std::endl;
                  for ( unsigned short j=0; j < gtd._nTurns; j++ ) {
                          ost << j+1 <<"-th loop point: "<<gtd._hitDataCollection[ gtd._turnFrstHit[j] ]._hitPoint << " momentum: "<< gtd._hitDataCollection[ gtd._turnFrstHit[j] ]._hitMomentum <<std::endl;
                  }
                  ost << "Track Hits:"<<std::endl;
//                  for ( size_t i=0; i < gtd._hitDataCollection.size(); i++ ) {
//                          ost << gtd._hitDataCollection[i] <<std::endl;
//                  }
                  for ( GenElHitDataCollection::const_iterator hd_it = gtd._hitDataCollection.begin(); hd_it!= gtd._hitDataCollection.end(); ++hd_it ) {
                          ost << hd_it->second <<std::endl;
                  }
                  return ost;
          }

          void addHitData( StrawHitPtr iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isOvrlpByProton_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ );/*{
                  GenElHitData hitData_
                  _hitDataCollection.push_back( hitData_ );
          }*/
          void sort();
          bool isNotPresent( size_t iHit_ );
          bool& isConversionEl( );
          void setConversionEl();
          //SimParticlePtr const& getSimPart();
          SimParticlePtr::key_type& getTrkKey();
          SimParticleCollection::key_type& getTrkID();
          CLHEP::Hep3Vector getTrkVertex();
          CLHEP::HepLorentzVector getTrkLrntzVec();
          size_t getNumOfHit();
          GenElHitData& getHit(/*size_t*/art::Ptr<SimParticle>::key_type ih);
          GenElHitData& getHit(int i_hit);
          unsigned short& FindNTurns();
          unsigned short& getNumOfLoops();
          GenElHitData& getithLoopHit(unsigned short ih);

  protected:
          //SimParticlePtr _trkPtr;
          SimParticlePtr::key_type _trkKey;
          SimParticleCollection::key_type _trkID;
          std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> _genTrackData;
          bool _isConversionEl;

          GenElHitDataCollection _hitDataCollection;
          GenElHitDCll_OrderByTime _orederedHitDataColl;
          //std::pair<CLHEP::Hep3Vector,CLHEP::HepLorentzVector> genTrackData;
          //std::vector hitDataCollection;
          //std::vector::iterator hd_it;
          //typename std::vector::iterator hd_it;
          bool _hitTimeSorted;
          unsigned short _nTurns;
          //std::vector< std::pair<CLHEP::Hep3Vector,CLHEP::Hep3Vector> > turnStartData;
          std::vector< art::Ptr<SimParticle>::key_type/*size_t*/ > _turnFrstHit;

          unsigned long rndup(double n);

  };

   void VisibleGenElTrack::addHitData( StrawHitPtr iHit_, double &mcHitTime_, bool &isOverlapped_, bool &isOvrlpByProton_, bool &isFirst_, CLHEP::Hep3Vector hitPoint_, CLHEP::Hep3Vector hitMomentum_ ) {
//          _hitDataCollection.push_back(
//                          GenElHitData( iHit_, mcHitTime_, isOverlapped_, isOvrlpByProton_, isFirst_, hitPoint_, hitMomentum_ )
//                          );
          _hitDataCollection.insert( GenElHitDataCollection::value_type( iHit_.key(),
                          GenElHitData( iHit_, mcHitTime_, isOverlapped_, isOvrlpByProton_, isFirst_, hitPoint_, hitMomentum_ ) )
                          );
          _orederedHitDataColl.insert( GenElHitDCll_OrderByTime::value_type( rndup(mcHitTime_*100.0), iHit_.key() ) );
  }

   void VisibleGenElTrack::sort(){
          //std::sort( _hitDataCollection.begin(), _hitDataCollection.end() );
          _hitTimeSorted=true;
  }

   bool VisibleGenElTrack::isNotPresent( size_t iHit_ ){
          //typename std::vector::iterator hd_it = _hitDataCollection.begin();
          for ( GenElHitDataCollection::iterator hd_it = _hitDataCollection.begin(); hd_it!= _hitDataCollection.end(); ++hd_it ) {
                  if(hd_it->first/*_iHit.key()*/==iHit_) return false;
          }
          return  true;
  }

   bool& VisibleGenElTrack::isConversionEl(){
           return _isConversionEl;
   }

   void VisibleGenElTrack::setConversionEl(){
           _isConversionEl=true;
   }

   //SimParticlePtr const& VisibleGenElTrack::getSimPart() {
   //       SimParticlePtr const &tmpPtr = _trkPtr;
   //       return tmpPtr;
   //}

   SimParticlePtr::key_type& VisibleGenElTrack::getTrkKey(){
          return _trkKey;
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

   GenElHitData& VisibleGenElTrack::getHit(art::Ptr<SimParticle>::key_type/*size_t*/ ih){
          //if ( ih<_hitDataCollection.size() ) return _hitDataCollection[ih];
          GenElHitDataCollection::iterator hd_it=_hitDataCollection.find(ih);
          if ( hd_it!=_hitDataCollection.end() ) return hd_it->second;
          else throw cet::exception("RANGE") << "Element "<<ih<<" out of size of the HitDataCollection "<<_hitDataCollection.size()<<std::endl;
  }

   GenElHitData& VisibleGenElTrack::getHit(int i_hit){
          if ( ((unsigned int)i_hit)<_hitDataCollection.size() ) {
                  GenElHitDataCollection::iterator hd_it=_hitDataCollection.begin();
                  for (int i=0; i<i_hit; i++) ++hd_it;
                  return hd_it->second;
          }
          else throw cet::exception("RANGE") << "Element "<<i_hit<<" out of size of the HitDataCollection "<<_hitDataCollection.size()<<std::endl;
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
                  _turnFrstHit.push_back(_orederedHitDataColl.begin()->second);
                  return _nTurns;
          }
          if ( !_hitTimeSorted ) sort();
//          for ( size_t ih=1 ; ih<_hitDataCollection.size(); ih++){
//                  if ( (_hitDataCollection[ih]._mcHitTime - _hitDataCollection[ih-1]._mcHitTime) > 1.0 ) {
//                          _nTurns++;
//                          _turnFrstHit.push_back(ih);
//                  }
//          }
          _nTurns=1;
          GenElHitDCll_OrderByTime::iterator prev_hitto_it=_orederedHitDataColl.begin();
          _turnFrstHit.push_back(prev_hitto_it->second);

          GenElHitDCll_OrderByTime::iterator hitto_it=prev_hitto_it;
          hitto_it++;
          for ( ; hitto_it!=_orederedHitDataColl.end(); ++hitto_it ){
//                  if ( (hitto_it->first - prev_hitto_it->first) > 1.0 ) {
                  if ( (hitto_it->first - prev_hitto_it->first) > 100 ) {
                          _nTurns++;
                          _turnFrstHit.push_back(hitto_it->second);
                  }
                  prev_hitto_it=hitto_it;
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

   unsigned long VisibleGenElTrack::rndup(double n)//round up a float type and show one decimal place
   {
           double t;
           t=n-floor(n);
           if (t>=0.5)
           {
                   n*=10.00000;//where n is the multi-decimal float
                   n=ceil(n);
                   n/=10.00000;
           }
           else
           {
                   n*=10.00000;//where n is the multi-decimal float
                   n=floor(n);
                   n/=10.00000;
           }
           return (unsigned long)n;
   }


}

#endif /* MCDataProducts_VisibleGenElTrack_hh */
