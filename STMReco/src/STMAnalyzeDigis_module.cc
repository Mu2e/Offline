// ======================================================================
//
// STMAnalyzeFragments_module:  Get and fit timing from STM Digis
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
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"

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
#include <TSpline.h>

namespace art
{
  class STMAnalyzeDigis;
}

using art::STMAnalyzeDigis;

// ======================================================================

class art::STMAnalyzeDigis : public EDAnalyzer
{
  public:
  using Name=fhicl::Name;
  using Comment=fhicl::Comment;
  struct Config {
    fhicl::Atom<art::InputTag> stmDigiCollection  { Name("stmDigiCollection"),  Comment("STM Digi module label") };
    fhicl::Atom<double> width_guess{ Name("width_guess"), Comment("Estimate of pulse width for fit (default = 150e-9)")};
    fhicl::Atom<double> rise_time_guess{ Name("rise_time_guess"), Comment("Estimate of pulse rise time for fit (default = 6e-9)")};
    fhicl::Atom<double> threshold{ Name("threshold"), Comment("Threshold to define the peak (default = 20000)")};
    fhicl::Atom<double> samp_freq{ Name("samp_freq"), Comment("Sampling frequency of ADC (default = 300e6)")};
    fhicl::Atom<int> pulse_region{ Name("pulse_region"), Comment("Number of ADC values to fit around pulse (default 51)")};
    fhicl::Atom<int> pulses_per_event{Name("pulses_per_event"), Comment("Number of pulses per event to find Delta t (default = 2)")};
  };
  using Parameters = art::EDAnalyzer::Table<Config>;

  // --- C'tor/d'tor:
  explicit STMAnalyzeDigis(const Parameters& config);

  // --- Production:
  virtual void analyze(const Event&);

  private:

  int event_counter = 0;
  
  // ROOT Service
  art::ServiceHandle<art::TFileService> tfs;

  // ROOT products
  TCanvas* c_stm;
  TH1F* DeltaT;
  TGraph* fEvent;
  TGraph* fitFirstPulse;
  TGraph* splinePulse;
  
  const  art::ProductToken<mu2e::STMWaveformDigiCollection> _stmDigisToken;
  // Fhicl params
  double _width_guess;
  double _rise_time_guess;
  double _threshold;
  double _samp_freq;
  int _pulse_region;
  int _pulses_per_event;

  // Averaging params
  // Pulse time average num
  double avg_num = 0;
  // Pulse time average den
  double avg_den = 0;
  
  // Boolean to plot the first event waveform
  bool firstEvent = true;
  // Boolean to plot the first pulse
  bool firstPulse = true;
  // Params to store end of last event
  int lastEventLength = 0;
  std::vector<double> prevEvent;
  
  // Define functions
  void book_histograms(art::ServiceHandle<art::TFileService> tfs);
  static double rising_edge(double*t, double*p);
  static double osc_rising_edge(double*t, double*p);
  double fit_rising_edge(std::vector<double> x, std::vector<double> y, double pulse_time);
  std::vector<double> linspace(double start_in, double end_in, int num_in);

}; // STMAnalyzeDigis

// ======================================================================

STMAnalyzeDigis::STMAnalyzeDigis(const Parameters& config) :
    art::EDAnalyzer{config}
    ,_stmDigisToken {consumes<mu2e::STMWaveformDigiCollection>(config().stmDigiCollection())}
    ,_width_guess(config().width_guess())
    ,_rise_time_guess(config().rise_time_guess())
    ,_threshold(config().threshold())
    ,_samp_freq(config().samp_freq())
    ,_pulse_region(config().pulse_region())
    ,_pulses_per_event(config().pulses_per_event())
{
  // Book the histograms
  book_histograms(tfs);
  // Set the size of the vector
  prevEvent.resize(_pulse_region);
}

// ----------------------------------------------------------------------

void STMAnalyzeDigis::book_histograms(art::ServiceHandle<art::TFileService> tfs) {
    
  DeltaT = tfs->make<TH1F>("delta_t","Time between 1st and 2nd pulse per event",100,10009,10011);
  fEvent = tfs->makeAndRegister<TGraph>("fEvent","First event; time [us]; ADC");
  fitFirstPulse = tfs->makeAndRegister<TGraph>("fitFirstPulse","First fitted pulse; time [us]; ADC");
  splinePulse = tfs->makeAndRegister<TGraph>("splinePulse","TSpline3 of first fitted pulse; time [us]; ADC");
  c_stm = tfs->makeAndRegister<TCanvas>("c_stm", "c_stm");

}

// Rising edge function
double STMAnalyzeDigis::rising_edge(double* t, double* p){

  double offset = p[0];
  double amplitude = p[1];
  double rise_time = p[2];
  double t0 = p[3];
  return offset + amplitude * (1 - exp(-(t[0] - t0) / rise_time)) * (t[0] >= t0);

}

// Oscillating rising edge function
double STMAnalyzeDigis::osc_rising_edge(double* t, double* p){

  double offset = p[0];
  double amplitude = p[1];
  double rise_time = p[2];
  double t0 = p[3];
  double w = p[4];
  double phi = p[5];
  return offset + amplitude * (exp(-(t[0] - t0) / rise_time))*cos(w*t[0] - phi) * (t[0] >= t0);

}

// Fitting function
double STMAnalyzeDigis::fit_rising_edge(std::vector<double> x, std::vector<double> y, double pulse_time) {

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
  std::vector<double> fit_y(y.begin(), y.begin()+middle);
  // Get start point
  double start = fit_x.front();
  // Get end point
  double end = fit_x.back();

  // Fit first pulse
  int vsize = fit_x.size();
  fitFirstPulse->Set(fit_x.size());
  
  if(firstPulse){
    std::vector<double> xi;	  
    for(int i=0; i < vsize; i++){
      if(fit_y[i] > 15000) xi.push_back(fit_x[i]);
      fitFirstPulse->SetPoint(i, fit_x[i], fit_y[i]);
    }    
    // Get reduced start point
    double startx = xi.front();
    // Get reduced end point
    double endx = xi.back();
    TSpline3* interpPulse = new TSpline3("interpPulse", fitFirstPulse);
    int nPoints = xi.size() * 100;
    std::vector<double> xx(nPoints);
    std::vector<double> yy;	  
    // Find linspace in pulse region
    xx = linspace(startx, endx, nPoints); 
    splinePulse->Set(nPoints);
    for(int ii=0; ii < nPoints; ii++){
      double v = interpPulse->Eval(xx[ii]);
      yy.push_back(v);
      splinePulse->SetPoint(ii, xx[ii], v);
      std::cout << xx[ii] << ", " << v << "\n";
    }

    // Create function
    TF1* fnx = new TF1("fnx", osc_rising_edge, startx, endx, 6);
    double offset_guessx = *std::min_element(yy.begin(),yy.end()); // offset
    double amp_guessx = *std::max_element(yy.begin(),yy.end()); // amplitude
    double center_guessx = pulse_time; // pulse centre
    double t0_guessx = center_guessx - _width_guess / 2; // rising edge t0
    double w_guessx = 1;
    double phi_guessx = 1;    
    // Initial guesses for the parameters
    fnx->SetParameters(offset_guessx, amp_guessx, _rise_time_guess, t0_guessx, w_guessx, phi_guessx);
    fnx->SetNpx(1E5);
    fnx->SetParNames("Offset","A","#tau","t_{0}", "#omega", "#phi");
    fnx->SetLineColor(kRed);
    fnx->SetLineWidth(3);

    fitFirstPulse->Draw("A");
    splinePulse->Draw("SAME");
    splinePulse->SetLineWidth(3);
    splinePulse->Fit("fnx","EM");
    firstPulse = false;
  }
  
  // Fit all pulses
  TF1* fn = new TF1("fn", rising_edge, start, end, 4);
  // Initial guesses for the parameters
  fn->SetParameters(offset_guess, amp_guess, _rise_time_guess, t0_guess);
  fn->SetParNames("Offset","A","#tau","t_{0}");
  fn->SetLineWidth(3);
  
  // Draw and fit graph
  TGraph *g = new TGraph(fit_x.size(), &fit_x[0], &fit_y[0]);
  g->Draw("A*");
  g->Fit("fn","Q");
  double rise_time = fn->GetParameter(2);
  double t0 = fn->GetParameter(3);
  // Calculate rising edge
  double rising_edge_fit = t0 + rise_time;
  
  // Return rising edge
  return rising_edge_fit;
}

std::vector<double> STMAnalyzeDigis::linspace(double start_in, double end_in, int num_in)
{

  std::vector<double> linspaced;

  double start = start_in;
  double end = end_in;
  double num = num_in;

  if (num == 0) { return linspaced; }
  if (num == 1) 
    {
      linspaced.push_back(start);
      return linspaced;
    }

  double delta = (end - start) / (num - 1);

  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end); // I want to ensure that start and end
                            // are exactly the same as the input
  return linspaced;
}

void STMAnalyzeDigis::analyze(const Event& event)
{
  //art::EventNumber_t eventNumber = event.event();
  
  // Pulse time average num
  //double avg_num = 0;
  // Pulse time average den
  //double avg_den = 0;

  // Number of pulses
  int pulse_num = 0;
  // Number of events
  int event_count = 0;

  auto stmH = event.getValidHandle<mu2e::STMWaveformDigiCollection>(_stmDigisToken);
  const mu2e::STMWaveformDigiCollection *stmDigis = stmH.product(); 

  if(stmDigis->size()!=0){
    
    /*std::cout << std::dec << "Analyzer: Run " << event.run() << ", subrun " << event.subRun()
	      << ", event " << eventNumber << " has " << std::endl;
	      std::cout << stmDigis->size() << " STM Digis" << std::endl;*/

    // Initialise vectors for pulse locations and fit vals
    std::vector<double> pulse_locs;
    std::vector<double> fit_vals;

    // Initialise first event vectors
    std::vector<double> t_vec;
    std::vector<double> adc_vec;
    
    int Digi_counter = 0;
    // Loop over the Digis
    for (auto& Digi : *stmDigis) {
      auto stm_Digi = static_cast<mu2e::STMWaveformDigi>(Digi);
      std::vector<int16_t> adcs = stm_Digi.adcs();
      // Plot the first waveform
      if(firstEvent){
	fEvent->Set(adcs.size());
	for(int ii=0; ii < adcs.size(); ii++){
	  double t = ii / _samp_freq;
	  fEvent->SetPoint(ii, t, adcs[ii]);
	}
	firstEvent=false;
      }
      // Get and fit the pulses
      for(int i=0; i < adcs.size(); i++){
	// If the data is above the pulse threshold
	if (adcs[i] > _threshold) {
	  //std::cout << "Index: " << i << "  |  ADCs: " << adcs[i] << "  |  ADCs[i+1]: " << adcs[i+1] << "\n"; 
	  // If the pulse started in the previous event
	  if (avg_den > 0 && i == 0) {
	    // Reset the average number...
	    avg_num = 0;
	    //std::cout << "Pulse in last event!: " << i << "   " << avg_den << "\n";
	    // Sum negative value from the previous event
	    for (int k = -avg_den; k < 0; ++k) {
	      avg_num += k;
	    }
	  } // End if pulse in previous event
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
	    double p_min = pulse_min / _samp_freq;
	    double p_max = pulse_max / _samp_freq;
	    //std::cout << "Pulse info :  pulse_num: " << pulse_num << "  | index: " << i << "  | loc:  " << pulse_loc << "  | min: " << pulse_min << "  | max: " << pulse_max << std::endl;
	    std::vector<double> x(_pulse_region);	  
	    // Find linspace in pulse region
	    x = linspace(p_min, p_max, _pulse_region); 
	    // y values to fit
	    std::vector<double> y;
	    for(int k=0; k < pulse_diff; k++){
	      y.push_back(adcs[pulse_min + k]);
	    }
	    // If the minimum pulse fit region is in the previous event
	    if (pulse_min < 0) {
	      //std::cout << "Pulse min in last event: " << pulse_min << "\n";
	      // get start value in last event to copy
	      int del = pulse_min + _pulse_region; 
	      // Store the y data that is in the previous event
	      std::vector<double> prev_event(prevEvent.begin() + del, prevEvent.end());
	      // Store the y data that is in this event
	      std::vector<double> this_event(adcs.begin(), adcs.begin() + pulse_max + 1);
	      // New y values to fit
	      y.insert(y.begin(), prev_event.begin(), prev_event.end());
	      y.insert(y.end(), this_event.begin(), this_event.end());
	    }
	    // Fit rising edge
	    double rising_edge_fit = fit_rising_edge(x, y, pulse_loc / _samp_freq);	  
	    //std::cout << "Rising edge fit: " << rising_edge_fit << "\n";
	    if(rising_edge_fit < 0){
	      // rising edge is in last event
	      double rise = rising_edge_fit;
	      rising_edge_fit = lastEventLength - rise;
	      // Save pulse time and fit vals in this event
	      pulse_locs.push_back(pulse_loc);
	      fit_vals.push_back(rising_edge_fit);	      
	    } else {
	      // Save pulse time and fit vals in this event
	      pulse_locs.push_back(pulse_loc);
	      fit_vals.push_back(rising_edge_fit);
	    }
	    // Save the length of the event	    
	    lastEventLength = adcs.size();
	    // Save last pulse region of event
	    int copyloc = adcs.size() - _pulse_region;
	    prevEvent.insert(prevEvent.begin(), adcs.begin() + copyloc, adcs.end());	  
	    // Increment number of pulses
	    pulse_num++;
	  }
	  // Reset average variables
	  avg_num = 0;
	  avg_den = 0;
	}
      }
      event_count++; 
      ++Digi_counter;
    }    
    //std::cout << "Pulse locs: " << pulse_locs.size() << "\n";
    // Find separation of first two consecutive pulses
    if(pulse_locs.size()>2){
      //double EWT = eventNumber;
      double pulse_0 = fit_vals[0];
      double pulse_1 = fit_vals[1];
      double pulse_sep = (pulse_1 - pulse_0)*1e9;
      DeltaT->Fill(pulse_sep);
      //std::cout << "Pulse sep: " << pulse_sep << "  |  EWT: " << EWT << "\n";
    }
    
    event_counter++;
    
  }      
} // analyze()

// ======================================================================

DEFINE_ART_MODULE(STMAnalyzeDigis)
