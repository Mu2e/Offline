#ifndef TrackerConditions_StrawElectronics_hh
#define TrackerConditions_StrawElectronics_hh

//
// StrawElectronics collects the electronics response behavior 
// of a Mu2e straw in several functions and parameters
//

// C++ includes
#include <iostream>
#include <vector>
#include <array>

// Mu2e includes
#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/StrawEnd.hh"
#include "DataProducts/inc/TrkTypes.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  class StrawElectronics : virtual public ProditionsEntity {
  public:
    // separately describe the 2 analog paths
    enum Path{thresh=0,adc,npaths};

    struct WireDistancePoint{
      double _distance;
      double _mean;
      double _normalization;
      double _sigma;
      double _t0;
      double _tmax[StrawElectronics::npaths]; // time at which value is maximum
      double _linmax[StrawElectronics::npaths]; // linear response to unit charge at max
      std::vector<double> _currentPulse;
      std::vector<double> _preampResponse;
      std::vector<double> _adcResponse;
      std::vector<double> _preampToAdc1Response;
      WireDistancePoint() {};
      WireDistancePoint(double distance, double mean, double normalization, 
			double sigma, double t0) : 
	_distance(distance), _mean(mean), _normalization(normalization), 
	_sigma(sigma), _t0(t0) {};
    };

    typedef std::shared_ptr<StrawElectronics> ptr_t;
    typedef std::shared_ptr<const StrawElectronics> cptr_t;

    StrawElectronics():_name("StrawElectronics") {}

    // construct with constants, then some values are computed and filled below
    StrawElectronics( double tdeadAnalog, double tdeadDigital, double vsat, 
		      double snoise, double ADCLSB, int maxADC, unsigned nADCpre, 
		      double ADCPeriod, double ADCOffset, 
		      unsigned maxtsep, unsigned tCoince, 
		      double TDCLSB, unsigned maxTDC, double TOTLSB, unsigned maxTOT, 
		      double tdcResolution, double electronicsTimeDelay, 
		      double ewMarkerROCJitter, double flashStart, 
		      double flashEnd, double flashClockSpeed, 
		      int responseBins, double sampleRate, int saturationSampleFactor, 
		      std::vector<double> preampPoles, std::vector<double> preampZeros, 
		      std::vector<double> adcPoles, std::vector<double> adcZeros, 
		      std::vector<double> preampToAdc1Poles, 
		      std::vector<double> preampToAdc1Zeros, 
		      std::vector<double> preampToAdc2Poles, 
		      std::vector<double> preampToAdc2Zeros, 
		      std::vector<double> wireDistances, 
		      std::vector<double> currentMeans, 
		      std::vector<double> currentNormalizations, 
		      std::vector<double> currentSigmas, 
		      std::vector<double> currentT0s, 
		      double clusterLookbackTime, 
		      std::vector<double> timeOffsetPanel, 
		      std::vector<double> timeOffsetStrawHV, 
		      std::vector<double> timeOffsetStrawCal) :
      _name("StrawElectronics"),
      _tdeadAnalog(tdeadAnalog), _tdeadDigital(tdeadDigital), _vsat(vsat), 
      _snoise(snoise), _ADCLSB(ADCLSB), _maxADC(maxADC), 
      _nADCpre(nADCpre), _ADCPeriod(ADCPeriod), _ADCOffset(ADCOffset), 
      _maxtsep(maxtsep), _TCoince(tCoince), 
      _TDCLSB(TDCLSB), _maxTDC(maxTDC), _TOTLSB(TOTLSB), _maxTOT(maxTOT), 
      _tdcResolution(tdcResolution), _electronicsTimeDelay(electronicsTimeDelay), 
      _ewMarkerROCJitter(ewMarkerROCJitter), 
      _flashStart(flashStart), _flashEnd(flashEnd), 
      _flashClockSpeed(flashClockSpeed),
      _responseBins(responseBins), 
      _sampleRate(sampleRate), _saturationSampleFactor(saturationSampleFactor), 
      _preampPoles(preampPoles), _preampZeros(preampZeros), _adcPoles(adcPoles), 
      _adcZeros(adcZeros), _preampToAdc1Poles(preampToAdc1Poles), 
      _preampToAdc1Zeros(preampToAdc1Zeros), _preampToAdc2Poles(preampToAdc2Poles), 
      _preampToAdc2Zeros(preampToAdc2Zeros), _wireDistances(wireDistances), 
      _currentMeans(currentMeans), _currentNormalizations(currentNormalizations), 
      _currentSigmas(currentSigmas), _currentT0s(currentT0s), 
      _clusterLookbackTime(clusterLookbackTime), 
      _timeOffsetPanel(timeOffsetPanel), 
      _timeOffsetStrawHV(timeOffsetStrawHV), 
      _timeOffsetStrawCal(timeOffsetStrawCal) {}

    virtual ~StrawElectronics() {}
    
    // linear response to a charge pulse.  This does NOT include saturation effects,
    // since those are cumulative and cannot be computed for individual charges
    double linearResponse(StrawId sid, Path ipath, double time,double charge,double distance,bool forsaturation=false) const; // mvolts per pCoulomb
    double adcImpulseResponse(StrawId sid, double time, double charge) const;
    // Given a (linear) total voltage, compute the saturated voltage
    double saturatedResponse(double lineearresponse) const;
    // relative time when linear response is maximal
    double maxResponseTime(Path ipath,double distance) const;
    // digization
    TrkTypes::ADCValue adcResponse(StrawId id, double mvolts) const; // ADC response to analog inputs
    TrkTypes::TDCValue tdcResponse(double time) const; // TDC response to a signal input to electronics at a given time (in ns since eventWindowMarker)
    void digitizeWaveform(StrawId id, TrkTypes::ADCVoltages const& wf,TrkTypes::ADCWaveform& adc) const; // digitize an array of voltages at the ADC
    bool digitizeTimes(TrkTypes::TDCTimes const& times,TrkTypes::TDCValues& tdc) const; // times in ns since eventWindowMarker
    bool digitizeAllTimes(TrkTypes::TDCTimes const& times,double mbtime, TrkTypes::TDCValues& tdcs) const; // for straws which are being read regardless of flash blanking
    void uncalibrateTimes(TrkTypes::TDCTimes &times, const StrawId &id) const; // convert time from beam t0 to tracker channel t0
    bool combineEnds(double t1, double t2) const; // are times from 2 ends combined into a single digi?
    // interpretation of digital data
    void tdcTimes(TrkTypes::TDCValues const& tdc, TrkTypes::TDCTimes& times) const;
    double adcVoltage(StrawId sid, uint16_t adcval) const; // mVolts
    double adcCurrent(StrawId sid, uint16_t adcval) const; // microAmps
    // accessors
    double adcLSB() const { return _ADCLSB; } //LSB in mvolts
    double tdcLSB() const { return _TDCLSB; } //LSB in nseconds
    double totLSB() const { return _TOTLSB; } //LSB in nseconds
    uint16_t maxADC() const { return _maxADC; }
    uint16_t maxTDC() const { return _maxTDC; }
    uint16_t maxTOT() const { return _maxTOT; }
    uint16_t ADCPedestal(StrawId sid) const { return _ADCped[sid.getStraw()]; };
    size_t nADCSamples() const { return TrkTypes::NADC; }
    size_t nADCPreSamples() const { return _nADCpre; }
    double adcPeriod() const { return _ADCPeriod; } // period of ADC clock in nsec
    double adcOffset() const { return _ADCOffset; } // offset WRT clock edge for digitization
    double flashStart() const { return _flashStart; } // time flash blanking starts
    double flashEnd() const { return _flashEnd; } // time flash blanking ends
    void adcTimes(double time, TrkTypes::ADCTimes& adctimes) const; // given crossing time, fill sampling times of ADC CHECK THIS IS CORRECT IN DRAC FIXME!
    double saturationVoltage() const { return _vsat; }
    double threshold(StrawId const &sid, StrawEnd::End iend) const { return _vthresh[sid.getStraw()*2 + iend]; }
    double analogNoise(Path ipath) const { return _analognoise[ipath]; }  // incoherent noise
    double strawNoise() const { return _snoise;} // coherent part of threshold circuit noise
    double deadTimeAnalog() const { return _tdeadAnalog; }
    double deadTimeDigital() const { return _tdeadDigital; }
    double TDCResolution() const { return _tdcResolution; }
    double electronicsTimeDelay() const { return _electronicsTimeDelay; }
    double eventWindowMarkerROCJitter() const { return _ewMarkerROCJitter; }
    
    double currentToVoltage(StrawId sid, Path ipath) const { return _dVdI[ipath][sid.getStraw()]; }
    double maxLinearResponse(StrawId sid, Path ipath,double distance,double charge=1.0) const;
    double normalization(Path ipath) const { return 1.;} //FIXME
    double fallTime(Path ipath) const { return 22.;} //FIXME
    double clusterLookbackTime() const { return _clusterLookbackTime;}
    
    double truncationTime(Path ipath) const { return _ttrunc[ipath];}
    double saturationTimeStep() const { return _saturationSampleFactor/_sampleRate;}
    void calculateResponse(std::vector<double> const& poles, 
		  std::vector<double> const& zeros, std::vector<double> const& input, 
                  std::vector<double> &response);

    std::string const& name() const { return _name; }
    void print(std::ostream& os) const;
    void printVector(std::ostream& os, std::string const& name, 
		    std::vector<double> const& a) const;


    // all of these must be called to fill this object
    void setdVdI(std::array<std::vector<double>,npaths> dVdI) { _dVdI = dVdI; }
    void setAnalogNoise(std::array<double,npaths> analognoise) 
              { _analognoise = analognoise; }
    void setvthresh(std::vector<double> vthresh) { _vthresh = vthresh;}
    void setFlashTDC( double flashStartTDC, double flashEndTDC) {
      _flashStartTDC = flashStartTDC;
      _flashEndTDC = flashEndTDC;
    }
    void setADCPed( std::vector<uint16_t> ADCped) { _ADCped = ADCped; }
    void setttrunc(std::array<double,npaths> ttrunc) { _ttrunc = ttrunc; }
    void setCurrentImpulse(std::vector<double> currentImpulse) {
      _currentImpulse = currentImpulse;
    }
    void setPreampToAdc2Response(std::vector<double> preampToAdc2Response) {
      _preampToAdc2Response = preampToAdc2Response;
    }
    void setwPoints(std::vector<WireDistancePoint> wPoints) {
      _wPoints = wPoints;
    }

    // this is used to update values from the database
    void setOffsets( std::vector<double> timeOffsetPanel,
		     std::vector<double> timeOffsetStrawHV,
		     std::vector<double> timeOffsetStrawCal ) {
      _timeOffsetPanel = timeOffsetPanel;
      _timeOffsetStrawHV = timeOffsetStrawHV;
      _timeOffsetStrawCal = timeOffsetStrawCal;
    }

  private:
    
    std::string _name;

    // generic waveform parameters
    std::array<std::vector<double>,npaths> _dVdI; // scale factor between charge and voltage (milliVolts/picoCoulombs)
    std::array<double,npaths> _ttrunc; // time to truncate signal to 0
    // threshold path parameters
    double _tdeadAnalog; // electronics dead time
    double _tdeadDigital; // electronics readout dead time
    // scale factor between current and voltage (milliVolts per microAmps)
    double _vsat; // saturation parameters.  _vmax is maximum output, _vsat is where saturation starts
    std::vector<double> _vthresh; // threshold voltage for electronics discriminator (mVolt)
    double _snoise; // straw noise at threshold
    std::array<double,npaths> _analognoise; //noise (mVolt) from the straw itself
    double _ADCLSB; // least-significant bit of ADC (mVolts)
    TrkTypes::ADCValue _maxADC; // maximum ADC value
    std::vector<uint16_t> _ADCped; // ADC pedestal (reading for 0 volts)
    size_t _nADCpre; // Number of ADC presamples
    double _ADCPeriod; // ADC period in nsec
    double _ADCOffset; // Offset of 1st ADC sample WRT threshold crossing (nsec)
    unsigned _maxtsep; // maximum # of ADC clock ticks between straw end threshold crossings to form a digi
    unsigned _TCoince; // maxing threshold xing pair time separation to create a digi, in number of ADC clock cycles
    double _TDCLSB; // least-significant bit of TDC (nsecs)
    TrkTypes::TDCValue _maxTDC; // maximum TDC value
    double _TOTLSB; // least-significant bit of TOT (nsecs)
    TrkTypes::TOTValue _maxTOT; // maximum TOT value
    double _tdcResolution; // tdc resolution (electronics effects only) (nsecs)
    double _electronicsTimeDelay; // Absolute time delay in electronics due to firmware signal propagation etc (ns)
    double _ewMarkerROCJitter; // jitter of ewMarker per ROC (ns)
    // electronicsTimeDelay is the time offset between a hit arriving at the electronics and the time that is digitized
    double _flashStart, _flashEnd, _flashClockSpeed; // flash blanking period (no digitizations during this time!!!) (ns from eventWindowMarker arrival, what will actually be set)
    TrkTypes::TDCValue _flashStartTDC, _flashEndTDC; // TDC values corresponding to the above. Note ignores electronicsTimeDelay since this is not a digitized signal
    // but an actual TDC value that will be compared against
    // helper functions
    static inline double mypow(double,unsigned);
    
    int _responseBins;
    double _sampleRate;
    int _saturationSampleFactor;
    std::vector<double> _preampPoles;
    std::vector<double> _preampZeros;
    std::vector<double> _adcPoles;
    std::vector<double> _adcZeros;
    std::vector<double> _preampToAdc1Poles;
    std::vector<double> _preampToAdc1Zeros;
    std::vector<double> _preampToAdc2Poles;
    std::vector<double> _preampToAdc2Zeros;
    
    std::vector<double> _wireDistances;
    std::vector<double> _currentMeans;
    std::vector<double> _currentNormalizations;
    std::vector<double> _currentSigmas;
    std::vector<double> _currentT0s;
    
    std::vector<double> _currentImpulse;
    std::vector<double> _preampToAdc2Response;
    
    std::vector<WireDistancePoint> _wPoints;
    
    double _clusterLookbackTime;
    
    std::vector<double> _timeOffsetPanel;
    std::vector<double> _timeOffsetStrawHV;
    std::vector<double> _timeOffsetStrawCal;
    
  };
  
}

#endif

