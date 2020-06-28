
// To run multithreader:  g++ -std=c++14 BasicMinimizer.cc -o a.out -pthread
#include <iostream>
#include <math.h>
#include <random>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <ctime>
#include <algorithm>
#include<thread>
#include<mutex>
#include <future>
#include "ThreaderPool.hh"
using namespace std;

bool use_multi = false;
bool mixed = true;
bool MDCprimaryNofilter = false;
bool MDCprimaryPrescale = false;
//System Details - from /proc/cpuinfo or lscpu can remain hardcoded once we know system (?)
unsigned N_cores_per_socket = 14;
unsigned N_sockets = 2;
unsigned N_threads_per_core = 2;
unsigned N_CPUs = 56;
unsigned N_Management_Threads =1; //used to oversee threading
unsigned MAX_THREADS = N_cores_per_socket*N_threads_per_core*N_sockets;
unsigned THREAD_COUNT = 1500;//MAX_THREADS;

bool diag = true;
unsigned int 	N_EVENTS;
unsigned int N_CONVERGED = 0;
unsigned int N_CRYSTALS = 674;
double Loss = 0;
double step_size = 0.0001;
double error = 1; //0.2671; //sigma of response.
double MaxIterations = 100;
double MaxFunction = 10;
double dcmin = 0.0001;
double dcmax = 0.1;
double max_c = 1.8;
double dFmax = 5;
double dSumMax = 0.1;

std::vector<double> CalibrationConstants;
std::vector<double> RawCalibrationResults;
std::vector<double> TrueConstants;
//std::vector<double> Vel;

struct CrystalList{
	std::vector<double> crystal_energy;
	std::vector<unsigned int> crystal_number;
	std::vector<double> crystal_energy_error;
  std::vector<double> crystal_corrected_energy;
  CrystalList();
  CrystalList(std::vector<double> _crystal_energy, std::vector<unsigned int> _crystal_number)
	: crystal_energy(_crystal_energy), crystal_number(_crystal_number) {};

	CrystalList(std::vector<double> _crystal_energy, std::vector<unsigned int> _crystal_number,
		std::vector<double> _crystal_energy_error)
	: crystal_energy(_crystal_energy), crystal_number(_crystal_number),
	crystal_energy_error(_crystal_energy_error){};

  CrystalList(std::vector<double> _crystal_energy, std::vector<unsigned int> _crystal_number,
		std::vector<double> _crystal_energy_error, std::vector<double> predicted)
	: crystal_energy(_crystal_energy), crystal_number(_crystal_number),
	crystal_energy_error(_crystal_energy_error), crystal_corrected_energy(predicted){};
};

struct Event{
  unsigned int EventNumber;
  double calo_energy;
  double track_energy;
  double track_mom;
	//double track_mom_err;
	unsigned int cluster_size;
	CrystalList crystal_list;
	Event();
	Event(unsigned int _n, double _energy, double _size) : EventNumber(_n),
	  track_energy(_energy), cluster_size(_size){};

	Event(unsigned int _n, double _energy, double _size, CrystalList _list )
	   : EventNumber(_n), track_energy(_energy), cluster_size(_size), crystal_list(_list){};

 Event(unsigned int _n, double _caloenergy, double _trackenergy,
	 		double _mom, double _size, CrystalList _list)
	   : EventNumber(_n), calo_energy(_caloenergy), track_energy(_trackenergy),
		 track_mom(_mom),  cluster_size(_size), crystal_list(_list){};
};

std::vector<Event> FakeDateMaker(std::vector<double> RawCalibrationResults, std::vector<double> offset_vector){

    std::vector<Event> event_list;

    double offset;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::normal_distribution<double> te(46.,3.);
    std::normal_distribution<double> cs(4.,1);
    std::uniform_real_distribution<double> cn(0, N_CRYSTALS);
    std::uniform_real_distribution<double> randoff(0.1, 1.0);

    for(unsigned int n=0;n<N_EVENTS;n++){

        //if(diag) std::cout<<"[In FakeDataMaker()] Event Loop ..."<<n<<std::endl;
				auto const track_energy = te(mt);
				auto const size = cs(mt);
        unsigned int cluster_size = round(size);

        std::vector<double> crystal_energy;
				std::vector<unsigned int> crystal_number;

				for(unsigned int m=0;m<cluster_size; m++){
            unsigned int C_number = round(cn(mt));
						crystal_number.push_back(C_number);
            offset = offset_vector[C_number];
            crystal_energy.push_back((1/offset)*track_energy/(cluster_size));
				}
				CrystalList crystal_list(crystal_energy,crystal_number);
				Event event(n, track_energy,cluster_size,crystal_list);
        event_list.push_back(event);
     }

    return event_list;

}

void SetOffsetVector(std::vector<double> &RawCalibrationResults, std::vector<double> &offset_vector){
    if(diag) std::cout<<"[In OffsetMaker()] Building Raw Values ..."<<std::endl;
    std::vector<Event> event_list;
    std::random_device rd;
    std::mt19937 mt(rd());
    //std::uniform_real_distribution<double> randoff(0.1, 1.0);
    std::normal_distribution<double> randoff(1.0,0.2);
    ofstream rawfile;
    rawfile.open("raw.csv");
    for(unsigned int c=0;c<N_CRYSTALS;c++){
       
        auto const off = randoff(mt);
        offset_vector.push_back(off);
        TrueConstants.push_back(off);
				RawCalibrationResults.push_back(off*0.8);
        rawfile<<off*0.8<<std::endl;
	 }
}

CrystalList FillCrystals(const char *CrystalsFile, unsigned int &n, unsigned int eventN){
		//if(diag) std::cout<<"[In ReadInputData()] : opening "<<CrystalsFile<<std::endl;

    FILE *fC = fopen(CrystalsFile, "r");
		if ( fC == NULL) {
			std::cout<<"[In ReadDataInput()] : Error: Cannot open Cry files "<<std::endl;
			exit(1);
		}

    float CryE, CryErr;
    unsigned int eventC, runC, CryId;

    std::vector<double> crystal_energy;
		std::vector<double> crystal_energy_error;
		std::vector<unsigned int> crystal_number;
    n=0;
		while(fscanf(fC, "%i,%i,%i,%f,%f\n", &eventC, &runC, &CryId, &CryE, &CryErr)!=EOF){
            if (eventC == eventN){
                crystal_energy.push_back(CryE*(1/TrueConstants[CryId]));
                crystal_number.push_back(CryId);
								crystal_energy_error.push_back(TrueConstants[CryId]*CryErr);
                n++;
            }
	}
    CrystalList crystal_list(crystal_energy, crystal_number, crystal_energy_error);
		fclose(fC);
    return crystal_list;
}

std::vector<Event> BuildEventsFromDataNew(const char *TrackFile, const char *predicteddE){
  
	std::vector<Event> event;
	FILE *fT = fopen(TrackFile, "r");
  FILE *fP = fopen(predicteddE, "r");
	if (fT == NULL or fP==NULL ) {
		std::cout<<"[In ReadDataInput()] : Error: Cannot open Track files "<<std::endl;
		exit(1);
	}

  float CaloE, TrackerE, P, PErr;
  unsigned int eventT, runT;

  float CryE, CryErr, dE;
  unsigned int CryId, clSize;

  std::vector<double> crystal_energy;
	std::vector<double> crystal_energy_error;
	std::vector<unsigned int> crystal_number;
  std::vector<double> corrected_crystal_energy;
  unsigned int currentEv=0, currentRun =0, count=0, EvCount=0;
  bool endEvent = false;
  //double correctiondE = 1;
  ofstream output;
  output.open("EnergySum.csv");
  std::vector<double> predicted;
  while(fscanf(fP, "%f\n", &dE)!=EOF){
    predicted.push_back(dE);
  }
  float Calop;
	while(fscanf(fT, "%i,%i,%i,%i,%f,%f,%f,%f,%f,%f\n", &eventT, &runT, &clSize, &CryId, &CryE, &CryErr,&CaloE, &TrackerE, &P, &PErr)!=EOF){
    
    if(currentEv == eventT and currentRun==runT) {
        crystal_energy.push_back(CryE*(1/TrueConstants[CryId]));
        crystal_number.push_back(CryId);
		    crystal_energy_error.push_back(CryErr);
        count++;
        Calop = predicted[EvCount];
        
         if(count==clSize){ 
            endEvent =true;
            EvCount +=1;  
         } 
        
      }

      if(currentEv != eventT or currentRun !=runT) {
        currentEv = eventT;
        currentRun = runT;
        crystal_energy.push_back(CryE*(1/TrueConstants[CryId]));
        crystal_number.push_back(CryId);
		    crystal_energy_error.push_back(CryErr);
        count ++;
        Calop = predicted[EvCount];
     
        if(count==clSize){ 
            endEvent =true;
            EvCount +=1;
         } 
        
      }
      if (endEvent) {
        CrystalList crystal_list(crystal_energy, crystal_number, crystal_energy_error);
	      Event e(eventT, Calop, TrackerE, P, clSize ,crystal_list);
        
        output<<TrackerE-CaloE<<endl;
        crystal_energy.clear();
        crystal_number.clear();
        crystal_energy_error.clear();
        event.push_back(e);
        count = 0;
        endEvent = false;
      }

	}
	
  fclose(fT);
	return event;
}

std::vector<Event> BuildEventsFromData(const char *CrystalsFile, const char *TrackFile){
	//if(diag) std::cout<<"[In ReadInputData()] : opening "<<TrackFile<<std::endl;
	std::vector<Event> event;
	FILE *fT = fopen(TrackFile, "r");

	if (fT == NULL ) {
		std::cout<<"[In ReadDataInput()] : Error: Cannot open Track files "<<std::endl;
		exit(1);
	}

  float CaloE, TrackerE, P, PErr;
  unsigned int eventT, runT;

	while(fscanf(fT, "%i,%i,%f,%f,%f,%f\n", &eventT, &runT, &CaloE, &TrackerE, &P, &PErr)!=EOF){
    unsigned int clSize;
    CrystalList crystal_list = FillCrystals(CrystalsFile, clSize, eventT);
	  Event e(eventT, CaloE, TrackerE, P, clSize ,crystal_list);
    event.push_back(e);

	}
	
  fclose(fT);
	return event;
}

unsigned int GetLines (const char *filename){
    ifstream file(filename);
    unsigned int lines=0;
    std::string line;
    while (getline(file, line)){
	   lines++;
    }
    if(diag) std::cout<<"[In GetLines()] File has "<<lines<<" Lines"<<std::endl;
    return lines;
}

std::vector<double> SGD(Event event, unsigned int j, std::vector<double> constants){

	double Loss = 0;
  double InitLoss = 0;
	double old_c ;
	double new_c ;
	double dc;
	double dFdCm;

	bool converged = false;
  
	if(diag) std::cout<<"[In SGD()] Initial Loss is "<<InitLoss<<std::endl;

	std::vector<double> previous_constants = constants;
	unsigned int k = 0;
	double Etrk = event.track_energy;
  double Ecalo = event.calo_energy;
  

  double Csum = 0; double Csum_km1 = 0; double dCsum_km1 = 0; double F_km1 = 0;
	while(converged == false and k < MaxIterations){
    Csum_km1 = Csum;
    Csum = 0;
  	for(unsigned int m=0; m < event.cluster_size; m++){
      
			unsigned int Cm = event.crystal_list.crystal_number[m];
			old_c = constants[Cm];
			
      double prediction = 0;

			double sigma = 1;//sqrt(1.68*1.68 +(old_c*event.crystal_list.crystal_energy_error[m]*event.crystal_list.crystal_energy_error[m]));
      double Esum=0;
			for(unsigned int i=0;i<event.cluster_size;i++){
          unsigned int Ci = event.crystal_list.crystal_number[i];
          prediction +=constants[Ci]*event.crystal_list.crystal_energy[i];
          Esum+=event.crystal_list.crystal_energy[i];
        
       }
      if(diag) cout<<"Calo E"<<Esum<<" tracker e "<<Etrk<<" calo e "<<Ecalo<<endl;
      double Vm = (event.crystal_list.crystal_energy[m]);
      Loss = pow((prediction - Etrk + Ecalo )/sigma ,2);
      if(k==0) {
        InitLoss = Loss;
        F_km1 = InitLoss;
      }
      dFdCm = 2*Vm*(1/(sigma))*(prediction -Etrk + Ecalo); //TODO - ignoring the sigma
      if(diag) cout<<"Grad Cm "<<dFdCm<<endl;
     
      new_c = (old_c - (1/sigma*sigma)*step_size*constants[Cm]*dFdCm);
      if (diag) cout<<"Old Cm "<<old_c<<" New Cm "<<new_c<<"True Cm "<<TrueConstants[Cm]<<endl;
      Csum +=new_c;
      dc = abs(new_c - old_c);
      if(!isnan(new_c ) and dc > dcmin and dc < dcmax and abs(new_c) < max_c){
        constants[Cm] = new_c;
        if(diag) cout<<"constants updated"<<endl;
      }
    }
    
    double dF = Loss-F_km1;
    double dCsum = Csum - Csum_km1;
    if (diag) cout<<"init  loss "<<InitLoss<<" new loss"<<Loss<<" previous "<<F_km1<<"change "<<dF<<endl;
    if (diag) cout<<"dc Csum "<<Csum<<" - "<<Csum_km1<<" = "<<dCsum<<" ddCsum "<<dCsum-dCsum_km1<<endl;
		if(k>0 and dCsum < dSumMax and Loss<MaxFunction and abs(dF)<dFmax){
      if (diag) cout<<"Converged with "<<k<<" iterations "<<endl;
      converged =true;
      N_CONVERGED +=1;
    }
    F_km1 = Loss;
    dCsum_km1 = dCsum;
		k++;
   
  }

	if(converged ==true){
   
		  for(unsigned int m=0; m<event.cluster_size;m++){
          int Ci = event.crystal_list.crystal_number[m];
          if (diag) cout<<"[In SGD() ] Updated "<<Ci<<" from "<<previous_constants[Ci]<<" to "<<constants[Ci]<<" with k iterations "<<k<<" True "<<TrueConstants[Ci]<<endl;
          
      }
	    return constants;
  } else{
      return previous_constants;
  }
}

void Randomize(std::vector<Event> &EventList){
    std::random_shuffle(EventList.begin(), EventList.end());
}

int main(int argc, char* argv[]){
    
     std::string thread_arg = argv[1];
     bool fake = false;
     std::cout<<"[In Main()] Beginning ..."<<std::endl;
     ofstream outputfile, TrackFile, CrystalsFile;
     outputfile.open("SGDv1.csv");
     outputfile<<"cryId,reco,true,Residual"<<std::endl;
     std::vector<double> offset_vector;


     for(unsigned int c=0;c<N_CRYSTALS;c++){
			CalibrationConstants.push_back(0);
      
		}

	  SetOffsetVector(RawCalibrationResults, offset_vector);
    std::vector<Event> event_list;
    if(fake) {
        N_EVENTS = 10000;
        event_list = FakeDateMaker(RawCalibrationResults, offset_vector);
    }
    if(!fake){
	    if(mixed) event_list = BuildEventsFromDataNew("/mu2e/data/users/sophie/mixed-IPA/GridOutcomes/Filtered200K/Combined.csv", "/mu2e/data/users/sophie/mixed-IPA/GridOutcomes/Filtered200K/predicteddEXboost.csv");
      /*if(MDCprimaryNofilter) event_list = BuildEventsFromDataNew("/mu2e/data/users/sophie/Simulations/IPAFromMDC/Combined.csv", "/mu2e/data/users/sophie/Simulations/IPAFromMDC/PredicteddE.csv");*/
      if(MDCprimaryNofilter) event_list = BuildEventsFromDataNew("/mu2e/data/users/sophie/Simulations/Primary/RecoResults/BigCombined.csv", "/mu2e/data/users/sophie/Simulations/Primary/RecoResults/predicteddEXBoost.csv");
      if(MDCprimaryPrescale) event_list = BuildEventsFromDataNew("/mu2e/data/users/sophie/primary-IPA/new_format/Prescaler/ForGrid/GridOutcomes/Combined.csv", "/mu2e/data/users/sophie/primary-IPA/new_format/Prescaler/ForGrid/GridOutcomes/predicteddE.csv");
			if(diag) std::cout<<"Found "<<event_list.size()<<" Events "<<std::endl;
     }
      N_EVENTS =  event_list.size();//GetLines("PreScaleTracks.csv");
     auto start = chrono::high_resolution_clock::now();
     CalibrationConstants = RawCalibrationResults;

    /*
      counting_barrier barrier(event_list.size());
      thread_pool Pool(THREAD_COUNT);
     */
    
     for(unsigned int l = 0; l < 1 ; l++){
          
         N_CONVERGED = 0;  
         if(diag) cout<<"Randomizing ... "<<endl;
	       Randomize(event_list);
         for(auto const& event : event_list){
            
            if(!use_multi or thread_arg == "--single"){
              CalibrationConstants = SGD(event, event.EventNumber, CalibrationConstants);
            }
         /*  if(use_multi or thread_arg == "--all"){
            auto done = Pool.add_task([&barrier, event=event, n= event.EventNumber, constants=CalibrationConstants]{

	            CalibrationConstants = SGD(event, n, constants);
               --barrier;
            });}*/
         }
     }
     /* if(use_multi || thread_arg == "--all"){
       barrier.wait();
     }*/
		 for(unsigned int i =0 ;i<N_CRYSTALS;i++){

        std::cout<<"End Offset "<<i<<" is "<<CalibrationConstants[i]
				<<" Raw "<<RawCalibrationResults[i] <<" True "<<TrueConstants[i]<<" Residuals "
				<<CalibrationConstants[i]-TrueConstants[i]<<" change "<<CalibrationConstants[i]-RawCalibrationResults[i]<<std::endl;
        outputfile<<i<<","<<CalibrationConstants[i]<<","
				<<offset_vector[i]<<","<<CalibrationConstants[i]-offset_vector[i]
				<<std::endl;

    }

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

    cout<<"NEvents Processed "<<N_EVENTS<<" NEVents converged "<<N_CONVERGED<<"Time  taken to converge "<<duration.count()<<endl;

	return 0;
}

