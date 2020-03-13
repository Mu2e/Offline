#ifndef ParticleTrack_h
#define ParticleTrack_h
#include "TROOT.h"
#include "DataProducts/inc/XYZVec.hh"

namespace mu2e{
	class ParticleTrack{
		ParticleTrack();
	        ParticleTrack(const ParticleTrack &);
	        ParticleTrack& operator=(const ParticleTrack &);
		//XYZVec const& _direction;
		int _pid;
		double _mom;
		int _nHits;
		double _startT;
		double _endT;

		bool _showElectrons;
		bool _showMuons;
		bool _showPions;
		bool _showGammas;

		public:
			ParticleTrack(double x0, double y0, double z0, double t0, double xN, double yN, double zN, double tN, int pid, double mom, bool showE, bool showMuon, bool showPion, bool showGamma) {
			_mom = mom;
			_pid = pid;
			_startT = t0;
			_endT = tN;
			_showElectrons = showE;
			_showMuons = showMuon;
			_showPions = showPion;
			_showGammas = showGamma;
			
		}
		
  		virtual ~ParticleTrack(){};

		//Some accessors:
		int GetPID(){return _pid;}
		double GetStartTime(){return _startT;}
		double GetEndTime(){return _endT;}
		int GetNHits(){ return _nHits;}

		//Simple Filters for track selection:
		bool FilterTrackHits(int minHits, int maxHits) {
			if(_nHits < minHits || _nHits > maxHits) { 
				return false; 
			}
			else { 
				return true;
			}
		}

		bool FilterTrackTime(double timeMin, double timeMax) {
			if(_startT < timeMin || _endT > timeMax) { 
				return false; 
			}
			else { 
				return true;
			}
		}

		bool FilterTrackPID(){ return true;}

		ClassDef(ParticleTrack,0);

	};

}//end namespace
#endif /*Particle_Track_h*/
