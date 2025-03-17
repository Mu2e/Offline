// ======================================================================
//
// STMWaveformDigisFromFragments: create STMWaveformDigis from STMFragments
//
// ======================================================================

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "art/Framework/Principal/Handle.h"
#include "artdaq-core-mu2e/Overlays/STMFragment.hh"
#include <artdaq-core/Data/Fragment.hh>

#include <string>
#include <memory>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <list>
#include <numeric>
#include <random>
#include <vector>
 
// ROOT includes
#include <TROOT.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>
#include <TLine.h>
#include <TGraph.h>

namespace art
{
  class STMWaveformDigisFromFragments;
}

using art::STMWaveformDigisFromFragments;

// ======================================================================

class art::STMWaveformDigisFromFragments : public EDProducer
{
  public:
  struct Config
  {
    fhicl::Atom<art::InputTag> stmTag {fhicl::Name("stmTag"), fhicl::Comment("Input module")};
    fhicl::Atom<double> width_guess{ fhicl::Name("width_guess"), fhicl::Comment("Estimate of pulse width for fit (default = 150e-9)")};
    fhicl::Atom<double> rise_time_guess{ fhicl::Name("rise_time_guess"), fhicl::Comment("Estimate of pulse rise time for fit (default = 6e-9)")};
    fhicl::Atom<double> threshold{ fhicl::Name("threshold"), fhicl::Comment("Threshold to define the peak (default = 20000)")};
    fhicl::Atom<double> samp_freq{ fhicl::Name("samp_freq"), fhicl::Comment("Sampling frequency of ADC (default = 300e6)")};
    fhicl::Atom<int> pulse_region{ fhicl::Name("pulse_region"), fhicl::Comment("Number of ADC values to fit around pulse (default 51)")};
    fhicl::Atom<int> pulses_per_event{fhicl::Name("pulses_per_event"), fhicl::Comment("Number of pulses per event to find Delta t (default = 2)")};
  };

  // --- C'tor/d'tor:
  explicit STMWaveformDigisFromFragments(const art::EDProducer::Table<Config>& config);

  // --- Production:
  virtual void produce(Event&);

  private:

  art::InputTag _stmFragmentsTag;
  
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
  
  // Params to store end of last event
  int lastEventLength = 0;
  uint64_t lastEWT = 0;
  std::vector<double> prevEvent;
  
  int eventCounter = 0;
  
  static double rising_edge(double*t, double*p);
  static double osc_rising_edge(double*t, double*p);
  double fit_rising_edge(std::vector<double> x, std::vector<double> y, double pulse_time);
  std::vector<double> linspace(double start_in, double end_in, int num_in);
  
}; // STMWaveformDigisFromFragments

// ======================================================================

// Rising edge function
double STMWaveformDigisFromFragments::rising_edge(double* t, double* p){

  double offset = p[0];
  double amplitude = p[1];
  double rise_time = p[2];
  double t0 = p[3];
  return offset + amplitude * (1 - exp(-(t[0] - t0) / rise_time)) * (t[0] >= t0);

}

// Fitting function
double STMWaveformDigisFromFragments::fit_rising_edge(std::vector<double> x, std::vector<double> y, double pulse_time) {

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

  int vsize = fit_x.size();

  // Fit all pulses
  TF1* fn = new TF1("fn", rising_edge, start, end, 4);
  // Initial guesses for the parameters
  fn->SetParameters(offset_guess, amp_guess, _rise_time_guess, t0_guess);
  fn->SetParNames("Offset","A","#tau","t_{0}");
  fn->SetLineWidth(3);
  
  // Draw and fit graph
  TGraph *g = new TGraph(vsize, &fit_x[0], &fit_y[0]);
  g->Draw("A*");
  g->Fit("fn","Q");
  double rise_time = fn->GetParameter(2);
  double t0 = fn->GetParameter(3);
  // Calculate rising edge
  double rising_edge_fit = t0 + rise_time;
  // Return rising edge
  return rising_edge_fit;
}

std::vector<double> STMWaveformDigisFromFragments::linspace(double start_in, double end_in, int num_in)
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

STMWaveformDigisFromFragments::STMWaveformDigisFromFragments(const art::EDProducer::Table<Config>& config) :
  art::EDProducer{config}
  ,_stmFragmentsTag(config().stmTag())
  ,_width_guess(config().width_guess())
  ,_rise_time_guess(config().rise_time_guess())
  ,_threshold(config().threshold())
  ,_samp_freq(config().samp_freq())
  ,_pulse_region(config().pulse_region())
  ,_pulses_per_event(config().pulses_per_event())
  
{
  // Set the size of the vector
  prevEvent.resize(_pulse_region);
  produces<mu2e::STMWaveformDigiCollection>();
}

// ----------------------------------------------------------------------


void STMWaveformDigisFromFragments::produce(Event& event)
{
  eventCounter++;
  std::unique_ptr<mu2e::STMWaveformDigiCollection> stm_waveform_digis(new mu2e::STMWaveformDigiCollection);
  
  //auto STMFragments = event.getValidHandle<artdaq::Fragments>("daq:STM");
  art::Handle<artdaq::Fragments> STMFragmentsH;
  if(!event.getByLabel(_stmFragmentsTag, STMFragmentsH)) {
    event.put(std::move(stm_waveform_digis));
    return;
  }
  const auto STMFragments = STMFragmentsH.product();

  int16_t detID = 0;
  for (const auto& frag : *STMFragments) {
    const auto& stm_frag = static_cast<mu2e::STMFragment>(frag);

    uint64_t EWT = *(stm_frag.EvNum());
    int16_t event_len = *(stm_frag.EvLen());    
    detID = *(stm_frag.detID()) & 0xFF;
    //std::cout << "DetID: " << *(stm_frag.detID()) << "  |   " << detID << "\n";
    //std::cout << "In event: " << event_number << " | Length of event: " << event_len << "\n";

    // Get 40MHz DTC clock time
    uint16_t DTC0 =  static_cast<uint16_t>((*(stm_frag.GetTHdr() + 28) >> 8) & 0xF) & 0xFFFF;
    uint16_t DTC1 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 29)) & 0xFFFF;
    uint16_t DTC2 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 30)) & 0xFFFF;
    uint16_t DTC3 =  static_cast<uint16_t>(*(stm_frag.GetTHdr() + 31)) & 0xFFFF;
    uint64_t DTC = static_cast<uint64_t>(DTC3) << 40 | static_cast<uint64_t>(DTC2) << 24 | static_cast<uint64_t>(DTC1) << 8 | static_cast<uint64_t>(DTC0);
    //std::cout << "DTCs: " << std::hex << DTC0 << " " << DTC1 << " " << DTC2 << " " << DTC3 << " --> " << DTC << std::dec << std::endl;
    //std::cout << "DTC time = " << DTC << " = " << (DTC/40e6) << "\n";

    // Get 75MHz ADC clock time
    int16_t adc[4] = {*(stm_frag.GetTHdr() + 4),*(stm_frag.GetTHdr() + 5),*(stm_frag.GetTHdr() + 6),*(stm_frag.GetTHdr() + 7)};
    uint64_t ADC;
    memcpy(&ADC, adc, sizeof(ADC));
    //std::cout << "ADC time = " << ADC << " = " << (ADC/75e6) << "\n";
    
    std::vector<int16_t> adcs;
    adcs.reserve(event_len);
    for (short int i_adc_sample = 0; i_adc_sample < event_len; ++i_adc_sample) {
      adcs.emplace_back(*(stm_frag.DataBegin()+i_adc_sample));
    }

    // Number of pulses
    int pulse_num = 0;

    // Initialise vectors for pulse locations and fit vals
    std::vector<double> pulse_locs;
    std::vector<double> fit_vals;
    
    // Get and fit the pulses
    for(long unsigned int i=0; i < adcs.size(); i++){
      // If the data is above the pulse threshold
      if (adcs[i] > _threshold) {
	// If the pulse started in the previous event
	if (avg_den > 0 && i == 0) {
	  // Reset the average number...
	  avg_num = 0;
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
	  if(rising_edge_fit > 0){
	    // rising edge is not in last event
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
    // Make sure variables are reset if events are not consecutive
    if(EWT != lastEWT+1){
      avg_num = 0;
      avg_den = 0;
    }

    double peak_fitTime1 = 0;
    double peak_fitTime2 = 0;
    double peak_sep = 0;
    // Find separation of first two consecutive pulses    
    if(pulse_locs.size()>2){
      peak_fitTime1 = fit_vals[0]*1e9;
      peak_fitTime2 = fit_vals[1]*1e9;
      peak_sep = (peak_fitTime2 - peak_fitTime1);
    }
    std::cout << "pulse 0 = " << peak_fitTime1 << "  |  pulse 1 = " << peak_fitTime2 << "  |  pulse sep = " << peak_sep << "\n";
  
    // Store EWT
    lastEWT = EWT;
    
    // Create the STMWaveformDigi and put it in the event
    uint32_t trig_time_offset = 0;
    mu2e::STMWaveformDigi stm_waveform(detID, EWT, DTC, ADC, trig_time_offset, peak_fitTime1, peak_fitTime2, peak_sep, adcs);
    stm_waveform_digis->push_back(stm_waveform);
  }

  event.put(std::move(stm_waveform_digis));

} // produce()

// ======================================================================

DEFINE_ART_MODULE(STMWaveformDigisFromFragments)
