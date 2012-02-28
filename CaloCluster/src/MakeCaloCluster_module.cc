/*
 * MakeCaloCluster_module.cc
 *
 *  Created on: Feb 10, 2012
 *      Author: gianipez
 */

//creazione del Producer per il clustering del calorimetro partendo dal file MakeStrawHit_module.cc
// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Handle.h"

// From the art tool-chain
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

//calorimeter packages
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "Mu2eUtilities/inc/sort_functors.hh"
//#include "CaloCluster/inc/CaloClusterUtilities.hh"
#include "CaloCluster/inc/CaloClusterer.hh"
//#include "CaloCluster/inc/CaloClusterFinder.hh"
//#include "CaloCluster/inc/ClosestCaloClusterFinder.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

// Other includes.
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib/exception.h"
#include "TMath.h"

using namespace std;

namespace mu2e {

//define a map which key is the index of the row "R" of the vane, and the contained object is an other map which key is the column index "Z" and also contain a vector.
//The vector contains pairs of "CaloCrystalHit" and a index which is the position of the "CaloCrystalHit" in the vector "CaloCrystalHitCollection" generated in the event
typedef std::map<unsigned int, std::map<unsigned int, std::vector<std::pair<CaloCrystalHit *, size_t > > >  > MatrixCaloHit;

//define a map which key is the vane's index and contains object of type "CaloCrystalHit". In that way we have a complete topology description of the calorimeter
typedef std::map<unsigned int, MatrixCaloHit> VanesMap;

//on the following we define a class, a vector and a map used to sort the "CaloCrystalHit" by energy
class CryID {

public:

        double _edep;//[MeV] energy deposited
        //unsigned int _vane;
        unsigned int _row;// "R" row
        unsigned int _colum;// "Z" column
        size_t _iCaloCrystalHit;//index of an element of the vector contained in "MatrixCaloHit"


        CryID(double edep,/*unsigned int van,*/ unsigned int r,  unsigned int c,  size_t index):
                _edep(edep),/*_vane(van),*/ _row(r), _colum(c),
                _iCaloCrystalHit(index) {}

        // This operator is overloaded in order to energy-sort the CrystalHits
        bool operator <(const CryID& b) const { return (_edep < b._edep); }

        bool operator ==(const CryID& b) const { return (_edep == b._edep && _row == b._row && _colum == b._colum && _iCaloCrystalHit == b._iCaloCrystalHit); }

//        inline std::ostream& operator<<( std::ostream& ost,
//                                          CryID & hit){
//           ost<<hit._row<<" , "<<hit._colum<<"  ,   "<< hit._edep<<"   ,   "<< hit._iCaloCrystalHit<<endl;
//           return ost;
//         }

};
//we use a set instead of a vector so we can use the "less" relation to fill that with energy order
typedef std::set<  CryID , less<CryID> > EnergyVec;

//the key of the following map is the vane's index
typedef std::map<unsigned int, EnergyVec> EnergyMap;


//define the object in which we store a single cluster. the key is the vane's index and the pair contains (row, column)
typedef std::multimap<unsigned int, std::pair<unsigned int, unsigned int> > ClusterData;//row, cloumn, hitId

class MakeCaloCluster : public art::EDProducer {

public:
        explicit MakeCaloCluster(fhicl::ParameterSet const& pset) :

        // Parameters
        _diagLevel(pset.get<int>("diagLevel",0)),
        _maxFullPrint(pset.get<int>("maxFullPrint",5)),
        _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)),
        _deltaTime(pset.get<double>("deltaTime", 100.)),// ns
        _nCryPerCluster(pset.get<int>("nCryPerCrystal", 0)),
        _EnoiseCut(pset.get<double>("EnoiseCut", 0.200)),//MeV 3 sigma noise
        _EclusterCut(pset.get<double>("EclusterCut", 0.200)),//MeV
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "energy")),
        _producerName("Algo"+TOUpper(_caloClusterAlgorithm)+"SeededBy"+TOUpper(_caloClusterSeeding)),
        _messageCategory("HITS"),

        // Control some information messages.
        _firstEvent(true){
                // Tell the framework what we make.
                produces<CaloClusterCollection>(_producerName);
        }
        virtual ~MakeCaloCluster() { }

        virtual void beginJob();

        void produce( art::Event& e);

private:

        // Diagnostics level.
        int _diagLevel;

        // Limit on number of events for which there will be full printout.
        int _maxFullPrint;

        // Name of the calorimeter StepPoint collection
        std::string _caloStepPoints;

        // Parameters
        double _minimumEnergy;  // minimum energy deposition of G4 step
        double _deltaTime;
        double _nCryPerCluster;
        double _EnoiseCut;
        double _EclusterCut;
        string _g4ModuleLabel;  // Name of the module that made these hits.
        string _caloReadoutModuleLabel;
        string _caloCrystalModuleLabel;
        string _caloClusterAlgorithm;
        string _caloClusterSeeding;
        const string _producerName;

        // A category for the error logger.
        const std::string _messageCategory;

        // Give some informationation messages only on the first event.
        bool _firstEvent;



};



void MakeCaloCluster::beginJob(){
        cout<<"selected clustering algorithm--> "<<_caloClusterAlgorithm <<", seeded by "<< _caloClusterSeeding<<endl;
        //        cout << "Diaglevel: "
        //                        << _diagLevel << " "
        //                        << _maxFullPrint
        //                        << endl;

}


void
MakeCaloCluster::produce(art::Event& evt) {


//        cout << "START PRODUCER MAKECALOCLUSTER..."<<endl;

        // Containers to hold the output information.
        auto_ptr<CaloClusterCollection>             caloClustersPointer(new CaloClusterCollection);

        CaloCluster Cluster;
        CaloCrystalHitPtrVector ptrtoHits;

        if ( _diagLevel > 0 ) cout << "MakeCaloCluster: produce() begin" << endl;

        static int ncalls(0);
        ++ncalls;

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;

        // Get handles to calorimeter collections
        art::Handle<CaloHitCollection> caloHits;
        evt.getByLabel(_caloReadoutModuleLabel, caloHits);

        art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
        evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

        //Get information on the calorimetr's geometry
        //int nVane = cg->nVane();       //starts from 1
        int nCryR = cg->nCrystalR();   //starts from 1
        int nCryZ = cg->nCrystalZ();   //starts from 1

        CaloClusterer clusterer(nCryR, nCryZ, _deltaTime, _nCryPerCluster, _EnoiseCut, _EclusterCut, _caloClusterAlgorithm);

        VanesMap vanesMap;

        EnergyMap energyMap;
        if (!( caloHits.isValid())) {
                return;
        }

        if (!caloCrystalHits.isValid()) {
                // cout << "NO CaloCrystalHits" << endl;
                return;
        }



        if(caloCrystalHits->size()>0 ){

                for(size_t i=0; i<caloCrystalHits->size(); ++i){
                        //cout << "start loop(caloCrystalHits->size())"<<endl;

                        CaloCrystalHit const& hit = (*caloCrystalHits).at(i);
                        //int cryId = hit.id();

                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();

                        if(ROIds.size()<1 ) continue;

                        if(hit.energyDep() < _minimumEnergy) continue;

                        CaloHit const& thehit = *ROIds.at(0);
                        int idVane = cg->getVaneByRO(thehit.id());//runs from 0 to nVane-1
                        //                        int idCry  = cg->getCrystalByRO(thehit.id());
                        int cryZ = cg->getCrystalZByRO(thehit.id());//get z-coordinate (from 0 to nCryZ-1)
                        int cryR = cg->getCrystalRByRO(thehit.id());//get r-coordinate (from 0 to nCryR-1)


                        vanesMap[idVane][cryR][cryZ].push_back(std::pair<CaloCrystalHit */*const&*/, size_t > (const_cast <CaloCrystalHit *>(&hit), i) );

                        energyMap[idVane].insert(CryID(hit.energyDep(), cryR, cryZ, vanesMap[idVane][cryR][cryZ].size() - 1  ));


                }//end caloCrystalHits loop
                //cout << "end loop(caloCrystalHits->size())"<<endl;

                VanesMap::iterator it = vanesMap.begin();

                unsigned int tmpRow;
                unsigned int tmpColum;
                unsigned int tmpHitId;

                if(_caloClusterSeeding.compare("ENERGY") == 0){
                        //cout <<"@@@@@@ ENERGY @@@@@@"<<endl;

                        while(it != vanesMap.end()){

                                //cout << "vane = "<< it->first <<", energyMap[it->first].size() = "<<energyMap[it->first].size() <<endl;
                                if(energyMap[it->first].size() == 0){
                                        ++it;
                                        continue;
                                }
                                EnergyVec::iterator   itMaxE = energyMap[it->first].end();
                                --itMaxE;

                                while(itMaxE != energyMap[it->first].end()){
                                        if(energyMap[it->first].size()!=0){
                                                itMaxE = energyMap[it->first].end();
                                                --itMaxE;
                                        }else{
                                                break;
                                        }
                                        if(itMaxE ==energyMap[it->first].end()){
                                                break;
                                        }
                                        tmpRow = itMaxE->_row;
                                        tmpColum =itMaxE->_colum;
                                        tmpHitId = itMaxE->_iCaloCrystalHit;

                                        //cout<< "it->second[tmpRow][tmpColum].size = " << it->second[tmpRow][tmpColum].size()<< endl;
                                        if(it->second[tmpRow][tmpColum].size()==0){
                                                CryID tmpCryID(itMaxE->_edep, tmpRow, tmpColum, tmpHitId);

                                                --itMaxE;
                                                energyMap[it->first].erase(tmpCryID);
                                                continue;
                                        }


                                        std::pair<CaloCrystalHit *, size_t > /*&*/itHit = it->second[tmpRow][tmpColum][tmpHitId];


                                        ClusterData cluster;
//                                        cout <<"itHit preso..."<< endl;
//                                        cout<<"itHit.first->time() = "<< itHit.first->time() <<endl;
//                                        cout<<"itHit.first->energyDep() = "<<itHit.first->energyDep()<<endl;

                                        clusterer.setFirstHitTime(itHit.first->time());

//                                        cout<<"clusterer.setFirstHitTime(itHit.first->time())------>done"<<endl;
//                                        cout<<"itHit.first->energyDep() = "<<itHit.first->energyDep()<<endl;

                                        clusterer.initializeNewSearch();

//                                        cout<<"clusterer.initializeNewSearch()---------->done"<<endl;
//                                        cout<<"itHit.first->energyDep() = "<<itHit.first->energyDep()<<endl;

                                        if( clusterer.find( cluster, it->second, tmpRow, tmpColum/*, looked*/ )){ // call to find cluster
                                               // cout << "-------- cluster found ----------"<<endl;

                                                caloClustersPointer->push_back(CaloCluster());
                                                CaloClusterCollection::iterator tmpCluster = caloClustersPointer->end();
                                                tmpCluster--;
                                                tmpCluster->_iVane = it->first;
                                                for(ClusterData::iterator itCD = cluster.begin(); itCD != cluster.end(); ++itCD){
                                                        //cout << "-------- 1 ----------"<<endl;
                                                        tmpRow = itCD->first;
                                                        tmpColum = itCD->second.first;
                                                        tmpHitId = itCD->second.second;

                                                        double tmpEnergy = it->second[tmpRow][tmpColum].at(tmpHitId).first->energyDep();
                                                        //double tmpTime = it->second[tmpRow][tmpColum].at(tmpHitId).first->time();
                                                        //tmpCluster->_energy += tmpEnergy;
                                                        //tmpCluster->_time += tmpTime;

                                                        //tmpCluster->_caloCrystalHits.push_back( CaloCrystalHitPtr( caloCrystalHits, it->second[tmpRow][tmpColum].at(tmpHitId).second ) );

                                                        CaloCrystalHitPtr tmpCrystalPtr( caloCrystalHits, it->second[tmpRow][tmpColum].at(tmpHitId).second );
                                                        tmpCluster->AddHit( tmpCrystalPtr );
                                                        //cout << "---------> AddHit done... <--------"<<endl;
                                                        CryID tmpCryID(tmpEnergy, tmpRow, tmpColum, tmpHitId);

                                                        it->second[tmpRow][tmpColum].erase( it->second[tmpRow][tmpColum].begin() + tmpHitId);

                                                        energyMap[it->first].erase(tmpCryID);

                                                        if( it->second[tmpRow][tmpColum].size() == 0 ){
                                                                //cout << "-------- it->second[tmpRow][tmpColum].size() == 0 ----------"<<endl;
                                                                it->second[tmpRow].erase(tmpColum);
                                                                //cout << "-------- 3 done erase ----------"<<endl;

                                                        }

                                                        if( it->second[tmpRow].size() == 0){
//                                                                cout << "-------- 4 it->second[tmpRow].size() == 0 ----------"<<endl;

                                                                it->second.erase(tmpRow);
//                                                                cout << "-------- 5 done erase ----------"<<endl;

                                                        }

                                                }//end for(ClusterData)
//                                                cout << "---------> end for(ClusterData)... <--------"<<endl;

                                                cog(*tmpCluster);
  //                                              cout<<"cog added..."<<endl;
                                                //                                                cout<< "cog.X = " <<tmpCluster->_impactPoint.getX() <<endl;
                                                //                                                cout<< "cog.Y = " <<tmpCluster->_impactPoint.getY() <<endl;
                                                //                                                cout<< "cog.Z = " <<tmpCluster->_impactPoint.getZ() <<endl;


                                        }else{
                                                CryID tmpCryID(itHit.first->energyDep(), tmpRow, tmpColum, tmpHitId);
//                                                cout<<"itHit.first->energyDep() = "<<itHit.first->energyDep()<<endl;
//                                                cout << "-------- cluster not found ----------"<<endl;
//                                                cout << "-------- it->second[tmpRow][tmpColum].size() = "<< it->second[tmpRow][tmpColum].size() <<endl;

                                                it->second[tmpRow][tmpColum].erase( it->second[tmpRow][tmpColum].begin() + tmpHitId);
//                                                cout << "-------- it->second[tmpRow][tmpColum].size() = "<< it->second[tmpRow][tmpColum].size() <<endl;
//
//                                                cout<< "energyMap[it->first].sixe = "<< energyMap[it->first].size() <<endl;

                                                EnergyVec::iterator j = energyMap[it->first].find(tmpCryID);
                                                if(j != energyMap[it->first].end()){
                                                        //cout<<"j_edep = "<< j->_edep<<endl;
                                                        --itMaxE;
                                                        energyMap[it->first].erase(  j  );

                                                }else{
                                                      //  cout<< "CryID not found*********************"<<endl;
                                                }




                                               // cout << "-------- 6 done erase ----------"<<endl;

                                                if( it->second[tmpRow][tmpColum].size() == 0 ){
                                                     //   cout << "-------- it->second[tmpRow][tmpColum].size() == 0 ----------"<<endl;

                                                    //    cout << "-------- it->second[tmpRow][tmpColum].size() = "<< it->second[tmpRow][tmpColum].size() <<endl;
                                                        it->second[tmpRow].erase( it->second[tmpRow].find(tmpColum) );

                                                }
                                                if( it->second[tmpRow].size() == 0){
                                                       // cout << "-------- it->second[tmpRow].size() = "<< it->second[tmpRow].size()<<endl;

                                                        //                                                it->second.erase(tmpRow);
                                                        it->second.erase(it->second.find(tmpRow));
//                                                        cout << "-------- 8 done erase ----------"<<endl;
//                                                        cout << "-------- it->second[tmpRow].size() = "<< it->second[tmpRow].size()<<endl;

                                                }

                                                if(energyMap[it->first].size() == 0) {
                                                      //  cout << "-------- energyMap[it->first].size() == 0 ----------"<<endl;

                                                        break;
                                                }
                                                //cout<<"itHit.first->energyDep() = "<<itHit.first->energyDep()<<endl;

                                                itMaxE = energyMap[it->first].end();
                                                --itMaxE;
                                                //cout<<"----------------------------------start 2 ---------------------------------------"<<endl;
                                                for(EnergyVec::iterator q = energyMap[it->first].begin(); q!= energyMap[it->first].end(); ++q){

                                                //        cout<<"(*q)._edep = "<<(*q)._edep<<", row = "<< (*q)._row<<", colum = "<<(*q)._colum<< ", caloCryHitSize = "<<(*q)._iCaloCrystalHit<<endl;
                                                }
                                               // cout<<"------------------------------------end 2 -------------------------------------"<<endl;

                                        }
                                }//end cry
                                //cout<<"end itEmax loop..."<<endl;
                                ++it;
                        }//end vanes
                      //  cout << "---------> end for(vanes)... <--------"<<endl;

                }
                else if(_caloClusterSeeding.compare("TIME") == 0 ) {
                        //cout <<"@@@@@@ TIME @@@@@@"<<endl;
                        while(it != vanesMap.end()){
                                MatrixCaloHit::iterator itrow = it->second.begin();
                                while(itrow != it->second.end()){
                                        std::map<unsigned int, std::vector < std::pair<CaloCrystalHit *, size_t > > >::iterator itcolum = itrow->second.begin();
                                        while(itcolum != itrow->second.end()){
                                                std::vector < std::pair<CaloCrystalHit *, size_t > >::iterator itHit = itcolum->second.begin();
                                                while(itHit != itcolum->second.end()){
                                                        ClusterData cluster;
                                                        //HitTimeMin = itHit->first->time();
                                                        clusterer.setFirstHitTime(itHit->first->time());
                                                        //                                                cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
                                                        //
                                                        //                                                cout << "HitTimeMin = " << itHit->first->time() <<endl;
                                                        //                                                                        for(int y=0; y<dim1; ++y){
                                                        //                                                                                for(int u=0; u<dim2; ++u){
                                                        //                                                                                        looked[y][u]=false;
                                                        //                                                                                }
                                                        //                                                                        }
                                                        clusterer.initializeNewSearch();
                                                        //                                              cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
                                                        //                                              cout<< "return bool(findClosesteCluster)"<<endl;
                                                        //                                              cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
                                                        if( clusterer.find( cluster, it->second, itrow->first, itcolum->first/*, looked*/ )){ // call to find cluster

                                                                caloClustersPointer->push_back(CaloCluster());
                                                                CaloClusterCollection::iterator tmpCluster = caloClustersPointer->end();

                                                                tmpCluster--;
                                                                tmpCluster->_iVane = it->first;

                                                                for(ClusterData::iterator itCD = cluster.begin(); itCD != cluster.end(); ++itCD){

                                                                        unsigned int tmpRow = itCD->first;
                                                                        unsigned int tmpColum = itCD->second.first;
                                                                        unsigned int tmpHitId = itCD->second.second;
                                                                        //                                                              cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
                                                                        //                                                              cout<<"itCD->first =  "<< tmpRow<<endl;
                                                                        //                                                              cout<<"itCD->secon.first =  "<< tmpColum<<endl;
                                                                        //                                                              cout<<"itCD->secon.second =  "<< tmpHitId<<endl;
                                                                        //                                                              cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
                                                                        //
                                                                        CaloCrystalHitPtr tmpCrystalPtr( caloCrystalHits, it->second[tmpRow][tmpColum].at(tmpHitId).second );
                                                                        tmpCluster->AddHit( tmpCrystalPtr );

                                                                        it->second[tmpRow][tmpColum].erase( it->second[tmpRow][tmpColum].begin() + tmpHitId);

                                                                        if( itcolum!=itrow->second.find(tmpColum) && it->second[tmpRow][tmpColum].size() == 0 ){
                                                                                it->second[tmpRow].erase(tmpColum);
                                                                                //cout<<"all colums erased too..." <<endl;
                                                                        }

                                                                }//end for(ClusterData)
                                                                cog(*tmpCluster);
                                                                //cout<<"end of FOR(CLUSTERDATA)"<<endl;
                                                        }else{

                                                                itcolum->second.erase(itHit);
                                                        }
                                                        itHit= itcolum->second.begin();
                                                }//end itHit
                                                //cout<<"itrow->second.size() =  "<< itrow->second.size() <<endl;
                                                if(itcolum->second.size()==0/*itrow->second.size() == 0*/ ){
                                                        itrow->second.erase(itcolum);
                                                        //cout<<" itrow->second.size() AFTER ERASE =  "<< itrow->second.size() <<endl;
                                                }
                                                itcolum = itrow->second.begin();
                                        }//end columns
                                        //cout<<" it->second.size() =  "<< itrow->second.size() <<endl;
                                        if(itrow->second.size() == 0 /*it->second.size() == 0*/ ){
                                                it->second.erase(itrow);
                                                //cout<<" it->second.size() AFTER ERASE =  "<< itrow->second.size() <<endl;
                                        }
                                        itrow = it->second.begin();

                                }//end rows
                                it++;
                        }//end vanes
                }else{
                        throw cet::exception("CALOCLUSTERING")
                        << "Could not find seeding procedure with the name: "
                        << _caloClusterSeeding
                        << endl;
                }

                evt.put(caloClustersPointer, _producerName);
                cout << "Event "<<evt.id().event()<<" CaloClustering done..."<<endl;
        }
} // end MakeCaloCluster::produce.


}// end namespace mu2e



using mu2e::MakeCaloCluster;
DEFINE_ART_MODULE(MakeCaloCluster);
