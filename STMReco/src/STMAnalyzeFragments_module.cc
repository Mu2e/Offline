// ======================================================================
//
// STMAnalyzeFragments_module:  Get and fit timing from STM Fragments
//
// ======================================================================

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"

#include "art_root_io/TFileService.h"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <list>
#include <numeric>
#include <random>
#include <vector>
 
// ROOT includes
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TBufferFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TString.h>
#include <TBufferJSON.h>
#include <TLine.h>
#include <TGraph.h>

namespace art
{
  class STMAnalyzeFragments;
}

using art::STMAnalyzeFragments;

// ======================================================================

class art::STMAnalyzeFragments : public EDAnalyzer
{
  public:
  using Name=fhicl::Name;
  using Comment=fhicl::Comment;
  struct Config {
    fhicl::Atom<double> width_guess{ Name("width_guess"), Comment("Estimate of pulse width for fit (default = 150e-9)")};
    fhicl::Atom<double> rise_time_guess{ Name("rise_time_guess"), Comment("Estimate of pulse rise time for fit (default = 6e-9)")};
    fhicl::Atom<double> threshold{ Name("threshold"), Comment("Threshold to define the peak (default = 20000)")};
    fhicl::Atom<double> samp_freq{ Name("samp_freq"), Comment("Sampling frequency of ADC (default = 300e6)")};
    fhicl::Atom<int> pulse_region{ Name("pulse_region"), Comment("Number of ADC values to fit around pulse (default 51)")};
    fhicl::Atom<int> pulses_per_event{Name("pulses_per_event"), Comment("Number of pulses per event to find Delta t (default = 2)")};
    fhicl::Atom<std::string> analysis_dir{Name("analysis_dir"), Comment("Directory to store output files (default = /home/mu2estm/sweetmor/pulses/)")};
    fhicl::Atom<std::string> data_file{Name("data_file"), Comment("Output file name (default = STM)")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  // --- C'tor/d'tor:
  explicit STMAnalyzeFragments(const Parameters& config);

  // --- Production:
  virtual void analyze(const Event&);

  private:

  int event_counter = 0;
  // Define vector for pulse separations 
  std::vector<double> pulse_seps;
  
  // ROOT Service
  art::ServiceHandle<art::TFileService> tfs;

  // ROOT products
  TCanvas* c_stm;
  TH1F* DeltaT;

  // Fhicl params
  double _width_guess;
  double _rise_time_guess;
  double _threshold;
  double _samp_freq;
  int _pulse_region;
  int _pulses_per_event;
  std::string _analysis_dir;
  std::string _data_file;

  // Define functions
  void book_histograms(art::ServiceHandle<art::TFileService> tfs);
  static double rising_edge(double*t, double*p);
  double fit_rising_edge(std::vector<double> x, std::vector<double> y, double pulse_time);
  std::vector<double> divide_element(std::vector<double> x);
  //void update_canvas();

}; // STMAnalyzeFragments

// ======================================================================

STMAnalyzeFragments::STMAnalyzeFragments(const Parameters& config) :
    art::EDAnalyzer{config}
    ,_width_guess(config().width_guess())
    ,_rise_time_guess(config().rise_time_guess())
    ,_threshold(config().threshold())
    ,_samp_freq(config().samp_freq())
    ,_pulse_region(config().pulse_region())
    ,_pulses_per_event(config().pulses_per_event())
    ,_analysis_dir(config().analysis_dir())
    ,_data_file(config().data_file())
{
  // Nothing here for now
}

// ----------------------------------------------------------------------

void STMAnalyzeFragments::book_histograms(art::ServiceHandle<art::TFileService> tfs) {
    
  DeltaT = tfs->make<TH1F>("delta_t","Time between pulses",30,10009,10011);
  c_stm = tfs->makeAndRegister<TCanvas>("c_stm", "c_stm");

}

// Rising edge function
double STMAnalyzeFragments::rising_edge(double* t, double* p){

  double offset = p[0];
  double amplitude = p[1];
  double rise_time = p[2];
  double t0 = p[3];
  return offset + amplitude * (1 - exp(-(t[0] - t0) / rise_time)) * (t[0] >= t0);

}

// Fitting function
double STMAnalyzeFragments::fit_rising_edge(std::vector<double> x, std::vector<double> y, double pulse_time) {

  // Function guess values
  double offset_guess = *std::min_element(y.begin(),y.end()); // offset
  double amp_guess = *std::max_element(y.begin(),y.end()); // amplitude
  double center_guess = pulse_time; // pulse centre
  double t0_guess = center_guess - _width_guess / 2; // rising edge t0

  // Find the middle of the pulse (same as middle of data)
  int middle = (x.size() / 2);
  // x values to middle 
  std::vector<double> fit_x(x.begin(), x.begin()+middle);
  // y values to middle 
  std::vector<double> fit_y(x.begin(), x.begin()+middle);

  // Create function
  TF1* fn = new TF1("fn", rising_edge, 0, middle, 4);
  // Initial guesses for the parameters
  fn->SetParameters(offset_guess, amp_guess, _rise_time_guess, t0_guess);

  TGraph *g = new TGraph(fit_x.size(), &fit_x[0], &fit_y[0]);
  g->Fit("fn");
  double rise_time = fn->GetParameter(2);
  double t0 = fn->GetParameter(3);
  // Calculate rising edge
  double rising_edge_fit = t0 + rise_time;
  
  // Return rising edge
  return rising_edge_fit;
}

std::vector<double> STMAnalyzeFragments::divide_element(std::vector<double> x){
  std::vector<double> vec;
  for(int i=0; i < x.size(); i++){
    double val = x[i] / _samp_freq;
    vec.push_back(val);
  }
  return vec;
}

void STMAnalyzeFragments::analyze(const Event& event)
{
  art::EventNumber_t eventNumber = event.event();

  // Pulse time average num
  double avg_num = 0;
  // Pulse time average den
  double avg_den = 0;
  // Pulse separation
  //double pulse_sep = 0;
  // Old average value
  //double prev_pulse_time = 0;
  // Number of pulses
  int pulse_num = 0;
  // Number of events
  int event_count = 0;
  
  auto stmH = event.getValidHandle<artdaq::Fragments>("daq:STM");
  const std::vector<artdaq::Fragment> *stmFrags = stmH.product(); 

  std::cout << std::dec << "Analyzer: Run " << event.run() << ", subrun " << event.subRun()
	    << ", event " << eventNumber << " has " << std::endl;
  std::cout << stmFrags->size() << " STM fragments." << std::endl;

  std::vector<double> pulse_locs;
  std::vector<double> fit_vals;
  
  int frag_counter = 0;
  for (auto& frag : *stmFrags) {

    auto stm_frag = static_cast<mu2e::STMFragment>(frag);
    //std::cout << "Frag #" << frag_counter << ": EvNum: " << *(stm_frag.EvNum()) << std::endl;
    //std::cout << "Frag #" << frag_counter << ": DataType: " << *(stm_frag.DataType()) << std::endl;
    //std::cout << "Frag #" << frag_counter << ": EvLen: " << *(stm_frag.EvLen()) << std::endl;
    //int16_t EWT = *(stm_frag.EvNum());
    for(int i=0; i < *(stm_frag.EvLen()); i++){
      // If the data is above the pulse threshold
      if (*(stm_frag.DataBegin() + i) > _threshold) {
	// Average pulse time numerator
	avg_num += i;
	// Average pulse time denonimator
	avg_den += 1;
      } else {
	// Else if not in a pulse region
	if (avg_num != 0) {
          // Calculate the pulse time average
          int avg_time = static_cast<int>(avg_num / avg_den);
          // Store current pulse location
          double pulse_loc = avg_time;
          // Minimum pulse fit region
          int pulse_min = static_cast<int>(pulse_loc - (_pulse_region - 1) / 2);
          // Maximum pulse fit region
          int pulse_max = static_cast<int>(pulse_loc + (_pulse_region - 1) / 2);
	  // Get diff
	  int pulse_diff = (pulse_max+1) - pulse_min;
          // x values to fit
	  //std::cout << "Pulse info 1: " << pulse_loc << "  " << pulse_min << "  " << pulse_max << std::endl;

	  std::vector<double> xa(_pulse_region);
          std::iota(xa.begin(), xa.end(), pulse_min / _samp_freq);
	  std::vector<double> x = divide_element(xa);
          // y values to fit
	  std::vector<double> y;
	  for(int k=0; k < pulse_diff; k++){
	    y.push_back(*(stm_frag.DataBegin() + pulse_min + k));
	  }
          // Fit rising edge
          double rising_edge_fit = pulse_loc;//fit_rising_edge(x, y, pulse_loc / _samp_freq);	  
	  std::cout << "Rising edge fit: " << rising_edge_fit << "\n";
	  // Save pulse time in this event
	  pulse_locs.push_back(pulse_loc);
	  //fit_vals.push_back(pulse_loc);
	  pulse_num++;
        }
	// Only reset average variables if a pulse was detected
        if (avg_den > 0) {
            avg_num = 0;
            avg_den = 0;
        }
      }
    }
    event_count++; 
    ++frag_counter;
  }

  std::cout << "Pulse locs: " << pulse_locs.size() << "\n";
  //std::cout << "Entering fit loop: " << fit_vals.size() << "\n";
  double pulse_0 = 0;
  if(pulse_locs.size()>0){
    // Loop over number of events
    for (int i = 0; i < pulse_locs.size(); i++) {
      // Get the EWT
      double EWT = eventNumber;
      if(i==0) pulse_0 = pulse_locs[i];
      if(i==1){
	double pulse_1 = pulse_locs[i];
	double pulse_sep = pulse_1 - pulse_0;
	std::cout << "Pulse sep: " << pulse_sep*1e9/300e6 << "  |  EWT: " << EWT << "\n";
	//pulse_seps.push_back(pulse_sep * 1e9);
      }
    }
  }
  
  std::cout << "Entering seps loop: " << pulse_seps.size() << "\n";
  if(pulse_seps.size() > 0){
    for(int i=0; i < pulse_seps.size(); i++){
      DeltaT->Fill(pulse_seps[i]);
      std::cout << "Index: " << i << "  |   Val: " << pulse_seps[i] << "\n";  
    }
    c_stm->cd(1);
    DeltaT->Draw();
    DeltaT->Fit("gaus");
    DeltaT->GetXaxis()->SetTitle("#Delta t [s]");
  }
  event_counter++;
      
} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMAnalyzeFragments)
