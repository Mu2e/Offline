/*
 * CaloClusterEff_module.cc
 *
 *  Created on: Feb 16, 2012
 *      Author: giani
 */



// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
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

//#include "Analyses/inc/MCCaloUtilities.hh"

// Mu2e includes.
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/VisibleGenElTrack.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"

//root includes
#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TGraph.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"

// From the art tool-chain
#include <cmath>
#include <deque>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include "TMath.h"

//#include "CaloCluster/inc/CaloClusterUtilities.hh"

using namespace std;

namespace mu2e {


struct elecData{
        double _cluEnergy;
        double _cluTime;
        double _impTime;
        double _impEnergy;
        CLHEP::Hep3Vector _cluCog;
        CLHEP::Hep3Vector _impPos;
        CLHEP::Hep3Vector _impPosCryFrame;
        unsigned int _vane;
        unsigned int _row;
        unsigned int _colum;

        bool operator<( const elecData other) const{
                return ( _impTime< other._impTime);
        }
        elecData & operator=(const elecData& other) {
                _cluEnergy = other._cluEnergy;
                _cluTime   = other._cluTime;
                _impTime = other._impTime;
                _impEnergy = other._impEnergy;
                _cluCog    = other._cluCog;
                _impPos = other._impPos;
                _impPosCryFrame = other._impPosCryFrame;
                _vane     = other._vane;
                _row      = other._row;
                _colum    = other._colum;
                return *this;
        }
        elecData():
                _impTime(1e10),
                _impEnergy(0.0){
        }
};

//the key is the the vane
typedef std::map<unsigned int,std::map<unsigned int, elecData > > ElecMap;


static int ncalls(0);

class CaloClusterCog : public art::EDAnalyzer {
public:
        explicit CaloClusterCog(fhicl::ParameterSet const& pset):
        _diagLevel(pset.get<int>("diagLevel",0)),
        _generatorModuleLabel(pset.get<std::string>("generatorModuleLabel", "generate")),
        _g4ModuleLabel(pset.get<std::string>("g4ModuleLabel", "g4run")),
        _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
        _caloCrystalModuleLabel(pset.get<std::string>("caloCrystalModuleLabel", "CaloCrystalHitsMaker")),
        _caloClusterModuleLabel(pset.get<std::string>("caloClusterModuleLabel", "makeCaloCluster")),
        _caloClusterAlgorithm(pset.get<std::string>("caloClusterAlgorithm", "closest")),
        _caloClusterSeeding(pset.get<std::string>("caloClusterSeeding", "time")),
        _producerName("Algo"+mu2e::TOUpper(_caloClusterAlgorithm)+"SeededBy"+mu2e::TOUpper(_caloClusterSeeding)),
        _elextractModuleLabel(pset.get<std::string>("elextractModuleLabel", "extractElData")),
        _extractElectronsData(pset.get<string>("elextractModuleLabel")),
        _rowToCanc(pset.get<int>("rowToCancel",-1)),//row index running from 1, 2, ... nCryR
        _columnToCanc(pset.get<int>("columnToCancel",-1)),//column index running from 1, 2, ... nCryZ
        _nAnalyzed(0),
        _hTHistGlobalEffNorm(0),
        _hTHistDeltaEnergy(0),
        _hTHistDeltaEnergyRec(0),
        _hTHistGlobalEffRec(0),
        _hTHistDeltaXquality(0),
        _hTHistDeltaXRec(0),
        _hTHistDeltaYquality(0),
        _hTHistDeltaYRec(0),
        _hTHistDeltaZquality(0),
        _hTHistDeltaZRec(0),
        _CogOffSet(pset.get<double>("cogOffset",3.)),
        _application(0),
        _directory(0)
        {
        }
        virtual ~CaloClusterCog() {
        }
        virtual void beginJob();
        virtual void endJob();

        void analyze(art::Event const& e );

private:

        void doCalorimeter(art::Event const& evt, bool skip);

        // Diagnostic level
        int _diagLevel;

        // Label of the generator.
        std::string _generatorModuleLabel;

        // Label of the G4 module
        std::string _g4ModuleLabel;

        // Label of the calo readout hits maker
        std::string _caloReadoutModuleLabel;

        // Label of the calo crystal hists maker
        std::string _caloCrystalModuleLabel;

        // Label of the calo clusters  maker
        std::string _caloClusterModuleLabel;

        string _caloClusterAlgorithm;
        string _caloClusterSeeding;
        string _producerName;
        string _elextractModuleLabel;
        string _extractElectronsData;

        double _minimumEnergy; //minimum energy deposition of hits
        int _rowToCanc, _columnToCanc;

        //number of analyzed events
        int _nAnalyzed;

        TH1D* _hTHistGlobalEffNorm;
        TH2D* _hTHistDeltaEnergy;
        TH2D* _hTHistDeltaEnergyRec;
        TH1D* _hTHistGlobalEffRec;

        TH2D* _hTHistDeltaXquality;
        TH2D* _hTHistDeltaXRec;
        TH2D* _hTHistDeltaYquality;
        TH2D* _hTHistDeltaYRec;
        TH2D* _hTHistDeltaZquality;
        TH2D* _hTHistDeltaZRec;

        double _CogOffSet;

        double globalNtrkCut;
        double globalNtrkTot;
        double *globalCaloCut;
        double *globalRecCaloCut;

        std::auto_ptr<MCCaloUtilities> CaloManager;

        bool _skipEvent;

        // The job needs exactly one instance of TApplication.  See note 1.
        auto_ptr<TApplication> _application;

        // Save directory from beginJob so that we can go there in endJob. See note 3.
        TDirectory* _directory;


};

bool findTrkId(std::vector<unsigned int> vec, unsigned int t){
        bool res = false;

        unsigned int size = vec.size();
        if(size!=0){
                unsigned int cont = 0;
                while(/*!res ||*/ cont!=size){
                        if(vec[cont] == t) {
                                res = true;
                        }
                        ++cont;

                }
        }
        return res;
}


void CaloClusterCog::beginJob( ) {

        cout << "start CaloClusterCog..."<<endl;

        CaloManager = auto_ptr<MCCaloUtilities>(new MCCaloUtilities());

        // If needed, create the ROOT interactive environment. See note 1.
        if ( !gApplication ){
                int    tmp_argc(0);
                char** tmp_argv(0);
                _application = auto_ptr<TApplication>(new TApplication( "noapplication", &tmp_argc, tmp_argv ));
        }

        gStyle->SetPalette(1);
        gROOT->SetStyle("Plain");

        _directory = gDirectory;


}

void CaloClusterCog::analyze(art::Event const & evt ) {

        ++_nAnalyzed;
        ++ncalls;

        CaloClusterer c;

        art::ServiceHandle<GeometryService> geom;
        if (ncalls == 1) {

                // cout << "This should be done only in the first event" << endl;

                art::ServiceHandle<art::TFileService> tfs;

                _hTHistGlobalEffNorm = tfs->make<TH1D>( "CaloGlobalEffTrkNorm", "CaloGlobalEffTrkNorm;EnergyLowCut [MeV];#epsilon", 100,  0.1, _CogOffSet);

                _hTHistDeltaEnergy   = tfs->make<TH2D>( "DeltaEnergy", "DeltaEnergy;EnergyLowCut [MeV];Eseed-Eclu [MeV]", 100,  0.1, _CogOffSet, 100, -20., 110.);
                _hTHistDeltaEnergyRec= tfs->make<TH2D>( "DeltaEnergyRec", "DeltaEnergyRec;EnergyLowCut [MeV];Eseed-Eclu [MeV]", 100,  0.1, _CogOffSet, 100, -20., 110.);
                _hTHistGlobalEffRec  = tfs->make<TH1D>( "CaloGlobalEffRec", "CaloGlobalEffRec;EnergyLowCut [MeV];#epsilon", 100,  0.1, _CogOffSet);

                _hTHistDeltaXquality= tfs->make<TH2D>( "DeltaXquality", "DeltaXquality;OffSet;Xseed-Xclu [mm]",100,  0.1, _CogOffSet, 200, -100., 100.);
                _hTHistDeltaXRec= tfs->make<TH2D>( "DeltaXRec", "DeltaXRec;OffSet;Xseed-Xclu [mm]",100, 0.1, _CogOffSet, 200, -100., 100.);
                _hTHistDeltaYquality= tfs->make<TH2D>( "DeltaYquality", "DeltaYquality;OffSet;Yseed-Yclu [mm]", 100, 0.1, _CogOffSet, 200, -100., 100.);
                _hTHistDeltaYRec= tfs->make<TH2D>( "DeltaYRec", "DeltaYRec;OffSet;Yseed-Yclu [mm]", 100, 0.1, _CogOffSet,200, -100., 100.);
                _hTHistDeltaZquality= tfs->make<TH2D>( "DeltaZquality", "DeltaZquality;OffSet;Zseed-Zclu [mm]", 100, 0.1, _CogOffSet,200, -100., 100.);
                _hTHistDeltaZRec= tfs->make<TH2D>( "DeltaZRec", "DeltaZRec;OffSet;Zseed-Zclu [mm]", 100, 0.1, _CogOffSet,200, -100., 100.);

                globalCaloCut = new double [_hTHistDeltaXRec->GetNbinsX()];
                globalRecCaloCut= new double [_hTHistDeltaXRec->GetNbinsX()];
                for(int j= 0; j< _hTHistDeltaXRec->GetNbinsX(); ++j ){

                        globalCaloCut[j] = 0.0;
                        globalRecCaloCut[j] = 0.0;
                }

        }

        doCalorimeter(evt, _skipEvent);

} // end of analyze

void CaloClusterCog::endJob() {

        for(unsigned int b=1; b <=(unsigned int) _hTHistGlobalEffNorm->GetNbinsX(); ++b){

                double tmpEff2 = globalCaloCut[b-1] / globalNtrkCut;
                _hTHistGlobalEffNorm->SetBinContent(b,  tmpEff2);
                _hTHistGlobalEffNorm->SetBinError(b, TMath::Sqrt(tmpEff2*fabs(1.0-tmpEff2) / globalNtrkCut ));

                double tmpEff4 = globalRecCaloCut[b-1] / globalNtrkTot;
                _hTHistGlobalEffRec->SetBinContent(b,  tmpEff4);
                _hTHistGlobalEffRec->SetBinError(b, TMath::Sqrt(tmpEff4*fabs(1.0-tmpEff4) / globalNtrkTot ));

        }

}



void CaloClusterCog::doCalorimeter(art::Event const& evt, bool skip){


        if ( _diagLevel > 0 ) cout << "MakeCaloCluster: produce() begin" << endl;

        //Get handle to calorimeter
        art::ServiceHandle<GeometryService> geom;
        if(! geom->hasElement<Calorimeter>() ) return;
        GeomHandle<Calorimeter> cg;

        // Get handles to calorimeter collections
        art::Handle<CaloHitCollection> caloHits;
        evt.getByLabel(_caloReadoutModuleLabel, caloHits);

        art::Handle<CaloCrystalHitCollection>  caloCrystalHits;
        evt.getByLabel(_caloCrystalModuleLabel, caloCrystalHits);

        // Get the persistent data about pointers to StepPointMCs
        art::Handle<PtrStepPointMCVectorCollection> mcptrHandle;
        evt.getByLabel(_caloReadoutModuleLabel,"CaloHitMCCrystalPtr",mcptrHandle);

        PtrStepPointMCVectorCollection const* hits_mcptr = mcptrHandle.product();
        if (!( caloHits.isValid())) {
                return;
        }

        if (!caloCrystalHits.isValid()) {
                cout << "NO CaloCrystalHits" << endl;
                return;
        }

        // Get handles to the generated and simulated particles.
        art::Handle<GenParticleCollection> genParticles;
        evt.getByLabel(_generatorModuleLabel, genParticles);

        art::Handle<SimParticleCollection> simParticles;
        evt.getByLabel(_g4ModuleLabel, simParticles);

        art::Handle<CaloClusterCollection> caloClusters;
        evt.getByLabel(_caloClusterModuleLabel,_producerName,caloClusters );

        auto_ptr<CaloClusterCollection>             caloClustersPointer(new CaloClusterCollection);

        art::Handle<VisibleGenElTrackCollection> genEltrksHandle;
        evt.getByLabel(_extractElectronsData,genEltrksHandle);

        VisibleGenElTrackCollection const* genEltrks = genEltrksHandle.product();
        std::vector<mu2e::VisibleGenElTrack>::const_iterator genEltrk_it;

        double trkMomCut = 100.0;//MeV
        std::vector<unsigned int> tmpV, tmpVTot;
        int NtrkCut =0;
        int NtrkTot = 0;

        //mapping &counting the electrons with quality cuts in the TRK
        for ( genEltrk_it = genEltrks->begin(); genEltrk_it!= genEltrks->end(); ++genEltrk_it ){
                ++NtrkTot;
                VisibleGenElTrack &iEltrk = const_cast<VisibleGenElTrack &>(*genEltrk_it);
                GenElHitData& hdil = iEltrk.getithLoopHit(0);

                if(!findTrkId(tmpVTot, iEltrk.getTrkID().asUint() ) ){

                        //cout<<"faccio il puschback in tmpV..."<<endl;
                        tmpVTot.push_back( iEltrk.getTrkID().asUint() );
                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                }

                if (iEltrk.getNumOfHit()>=20){



                        if(hdil._hitMomentum.mag() >= trkMomCut){
                                NtrkCut++;// # of electrons at the entrance of the TRK which have trkMomCut MeV of kinetic energy and which will light at least 20 straws
                                if(!findTrkId(tmpV, iEltrk.getTrkID().asUint() ) ){

                                        //cout<<"faccio il puschback in tmpV..."<<endl;
                                        tmpV.push_back( iEltrk.getTrkID().asUint() );
                                        //cout<<"--------> puschback in tmpV eseguito..."<<endl;

                                }

                                //                                for(unsigned int b=1; b <=(unsigned int) _hTHistEff->GetNbinsX(); ++b){
                                //
                                //                                        if (hdil._hitMomentum.mag() >= (_hTHistEff->GetXaxis()->GetBinCenter(b) -  0.5*_hTHistEff->GetXaxis()->GetBinWidth(b) )){
                                //                                               // NGoodElec[b-1]++;
                                //                                                if ( iEltrk.isConversionEl() ) {
                                //                                                     //   NGoodConvEl[b-1]++;
                                //                                                }
                                //                                        }
                                //
                                //                                }
                        }

                }

        }//end TRK mapping
        //cout <<"NtrkCut = " <<NtrkCut<<endl;

        if (NtrkCut==0) return;

        globalNtrkCut += NtrkCut;
        globalNtrkTot += NtrkTot;

        for(int bin =1; bin <= _hTHistDeltaXRec->GetNbinsX(); ++bin){
                //for(int bin =1; bin <= _hTHistEff->GetNbinsX(); ++bin ){
                std::vector<unsigned int> trkVec, trkVecTot;
                for(unsigned int f = 0; f != tmpV.size(); ++f){
                        trkVec.push_back(tmpV[f]);
                }
                for(unsigned int f = 0; f != tmpVTot.size(); ++f){
                        trkVecTot.push_back(tmpVTot[f]);
                }
                ElecMap elecMap;

                if(caloClusters->size()>0 ){
                        int iVane;
                        for(size_t icl=0; icl<caloClusters->size(); ++icl){
                                CaloCluster const& clu = (*caloClusters).at(icl);

                                caloClustersPointer->push_back(clu);
                                CaloClusterCollection::iterator tmpCluster = caloClustersPointer->end();
                                tmpCluster--;


                                CLHEP::Hep3Vector logCog = LOGcog(*tmpCluster,_hTHistDeltaXRec->GetBinCenter(bin) );

                                iVane = clu._iVane;

                                CaloCrystalHitPtrVector caloClusterHits = clu._caloCrystalHits;

                                for(size_t i=0; i<caloClusterHits.size(); ++i){
                                        CaloCrystalHit const& hit = *(caloClusterHits.at(i));

                                        std::vector<art::Ptr<CaloHit> > const& ROIds = hit.readouts();
                                        if(ROIds.size()<1 ) continue;

                                        CaloHit const& thehit = *ROIds.at(0);
                                        size_t collectionPosition = ROIds.at(0).key();

                                        PtrStepPointMCVector const & mcptr(hits_mcptr->at(collectionPosition));
                                        size_t nHitsPerCrystal = mcptr.size();

                                        for (size_t j2=0; j2<nHitsPerCrystal; ++j2) {

                                                //cout<< "Start loop..."<< "j2 = "<< j2<<endl;

                                                StepPointMC const& mchit = *mcptr[j2];

                                                // The simulated particle that made this hit.
                                                SimParticleCollection::key_type trackId(mchit.trackId());
                                                CLHEP::Hep3Vector cryFrame = cg->toCrystalFrame(thehit.id(), mchit.position());

                                                if( (cg->getCrystalRByRO(thehit.id()) != ( _rowToCanc - 1 ) && cg->getCrystalZByRO(thehit.id()) != ( _columnToCanc - 1 ) ) ){

                                                        if(elecMap[iVane][trackId.asUint()]._impTime > mchit.time() ){
                                                                elecMap[iVane][trackId.asUint()]._cluEnergy = clu._energy;
                                                                elecMap[iVane][trackId.asUint()]._cluTime = clu._time;
                                                                elecMap[iVane][trackId.asUint()]._impTime = mchit.time();
                                                                elecMap[iVane][trackId.asUint()]._impEnergy = mchit.momentum().mag();
                                                                elecMap[iVane][trackId.asUint()]._cluCog = logCog;//clu._impactPoint;
                                                                elecMap[iVane][trackId.asUint()]._impPos = mchit.position();
                                                                elecMap[iVane][trackId.asUint()]._impPosCryFrame = cryFrame;
                                                                elecMap[iVane][trackId.asUint()]._vane = iVane;
                                                                elecMap[iVane][trackId.asUint()]._row = cg->getCrystalRByRO(thehit.id());
                                                                elecMap[iVane][trackId.asUint()]._colum = cg->getCrystalZByRO(thehit.id());

                                                                //                                                                cout<< "###################"<<endl;
                                                                //                                                                cout<< "idVande = "<< iVane<<endl;
                                                                //                                                                cout << "cluX = "<<elecMap[iVane][trackId.asUint()]._impPos.getX()<<endl;
                                                                //                                                                cout << "cluY = "<<elecMap[iVane][trackId.asUint()]._impPos.getY()<<endl;
                                                                //                                                                cout << "cluZ = "<<elecMap[iVane][trackId.asUint()]._impPos.getZ()<<endl;
                                                                //                                                                cout<< "###################"<<endl;


                                                        }
                                                }
                                        }

                                }


                        }//end for(caloClsuters.size)
                }//end if(caloClsuters.size() > 0)

                //counting how many reconstructed clusters match with the generated particles
                unsigned int canc2 = 0;

                //cout <<"-------------------> numero di elettroni generati = "<< trkVecTot.size() <<endl;
                unsigned int size2 = trkVecTot.size();
                unsigned int it2=0;
                while( it2 < size2){
                        ElecMap::iterator ite = elecMap.begin();
                        bool trovato = false;
                        while(!trovato && ite!=elecMap.end() ){
                                if(ite->second.find(trkVecTot[it2]) != ite->second.end()){
                                        //                                        cout<< "$$$$$$$$$$$$----> vane = "<< ite->first<<endl;
                                        //                                        cout<<"------------------------------------------------------"<<endl;
                                        //                                        cout << "impX = "<<ite->second[trkVec[it2]]._impPos.getX()<<endl;
                                        //                                        cout << "impY = "<<ite->second[trkVec[it2]]._impPos.getY()<<endl;
                                        //                                        cout << "impZ = "<<ite->second[trkVec[it2]]._impPos.getZ()<<endl;
                                        //                                        cout << "impXcryFRame = "<<ite->second[trkVec[it2]]._impPosCryFrame.getX()<<endl;
                                        //                                        cout << "impYcryFRame = "<<ite->second[trkVec[it2]]._impPosCryFrame.getY()<<endl;
                                        //                                        cout << "impZcryFRame = "<<ite->second[trkVec[it2]]._impPosCryFrame.getZ()<<endl;
                                        //                                        cout<< "@@@@@@@@@@@@@@@@@@@"<<endl;
                                        //cout << "cluX = "<<ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                        //                                        cout << "cluY = "<<ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                        //                                        cout << "cluZ = "<<ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
                                        //                                        cout<< "@@@@@@@@@@@@@@@@@@@"<<endl;
                                        //                                        cout << "delta X = " << ite->second[trkVec[it2]]._impPos.getX() - ite->second[trkVec[it2]]._cluCog.getX()<<endl;
                                        //                                        cout << "delta Y = " << ite->second[trkVec[it2]]._impPos.getY() - ite->second[trkVec[it2]]._cluCog.getY()<<endl;
                                        //                                        cout << "delta Z = " << ite->second[trkVec[it2]]._impPos.getZ() - ite->second[trkVec[it2]]._cluCog.getZ()<<endl;
                                        globalRecCaloCut[bin - 1] += 1.0;

                                        _hTHistDeltaEnergyRec->Fill(_hTHistDeltaEnergyRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impEnergy - ite->second[trkVecTot[it2]]._cluEnergy);

                                        if(ite->first == 0 || ite->first == 2){

                                                _hTHistDeltaXRec->Fill(_hTHistDeltaXRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getX() - ite->second[trkVecTot[it2]]._cluCog.getX() );
                                        }else{
                                                _hTHistDeltaYRec->Fill(_hTHistDeltaYRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getY() - ite->second[trkVecTot[it2]]._cluCog.getY() );
                                        }
                                        _hTHistDeltaZRec->Fill(_hTHistDeltaZRec->GetXaxis()->GetBinCenter(bin), ite->second[trkVecTot[it2]]._impPos.getZ() - ite->second[trkVecTot[it2]]._cluCog.getZ() );

                                        trovato = true;
                                        canc2 = it2;
                                        //cout<<"$$$ yes! the calo works! $$$"<<endl;
                                }
                                ++ite;
                        }

                        if(trovato){
                                std::vector<unsigned int>::iterator er = trkVecTot.begin();
                                er +=canc2;
                                trkVecTot.erase(er);
                                size2 = trkVecTot.size();
                        }else{
                                ++it2;
                        }

                }
                //                cout <<"after.... trkVecTot.size = "<< trkVecTot.size() <<endl;
                //-----------------------------------

                //counting how many reconstructed clusters match with the generated particles which pass the quality cuts of the tracker
                unsigned int canc = 0;
                //cout <<"-------------------> numero di elettroni di qualita' generati = "<< trkVec.size() <<endl;
                unsigned int size = trkVec.size();
                unsigned int it=0;
                while( it < size){
                        ElecMap::iterator ite = elecMap.begin();
                        bool trovato = false;
                        while(!trovato && ite!=elecMap.end() ){
                                if(ite->second.find(trkVec[it]) != ite->second.end()){
                                        //                                        cout << "impX = "<<ite->second[trkVec[it]]._impPosCryFrame.getX()<<endl;
                                        //                                        cout << "impY = "<<ite->second[trkVec[it]]._impPosCryFrame.getY()<<endl;
                                        //                                        cout << "impZ = "<<ite->second[trkVec[it]]._impPosCryFrame.getZ()<<endl;

                                        globalCaloCut[bin - 1] += 1.0;

                                        _hTHistDeltaEnergy->Fill(_hTHistDeltaEnergy->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._impEnergy - ite->second[trkVec[it]]._cluEnergy);
                                        //cout<< "cervavo "<< trkVec[it]<<", trovato "<< jo->first <<endl;
                                        trovato = true;
                                        canc = it;
                                        //cout<<"$$$ yes! the calo works! $$$"<<endl;
                                        if(ite->first == 0 || ite->first == 2){

                                                _hTHistDeltaXquality->Fill(_hTHistDeltaXquality->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._impPos.getX() - ite->second[trkVec[it]]._cluCog.getX() );
                                        }else{
                                                _hTHistDeltaYquality->Fill(_hTHistDeltaYquality->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._impPos.getY() - ite->second[trkVec[it]]._cluCog.getY() );
                                        }
                                        _hTHistDeltaZquality->Fill(_hTHistDeltaZquality->GetXaxis()->GetBinCenter(bin), ite->second[trkVec[it]]._impPos.getZ() - ite->second[trkVec[it]]._cluCog.getZ() );
                                        // }
                                }
                                ++ite;
                        }

                        if(trovato){
                                std::vector<unsigned int>::iterator er = trkVec.begin();
                                er +=canc;
                                trkVec.erase(er);
                                size = trkVec.size();
                        }else{
                                ++it;
                        }

                }

                //                cout <<"after.... trkVec.size = "<< trkVec.size() <<endl;
                //                cout<<"...\n...\n "<<endl;

        }//end for(_hTHistEff->GetNbinsX())
        cout << "Event "<<evt.id().event()<<" CaloClusterCog done..."<<endl;

}

}


using mu2e::CaloClusterCog;
DEFINE_ART_MODULE(CaloClusterCog);

