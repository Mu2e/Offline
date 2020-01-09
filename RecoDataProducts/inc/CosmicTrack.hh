#ifndef RecoDataProducts_CosmicTrack_hh
#define RecoDataProducts_CosmicTrack_hh
////S. Middleton, Feb 2019 - Cosmic track class, main purpose id to store diagnostics.
#include "TMath.h"
#include "TMatrixD.h"
#include "DataProducts/inc/XYZVec.hh"
#include "Mu2eUtilities/inc/PointLinePCA_XYZ.hh"
#include<vector>
#include<bitset>

using CLHEP::Hep3Vector;
using namespace std;
   //Struct To Hold Covarience Info:
   struct TrackCov{
   	std::vector<double> Covariance;
	double sigA0;
	double sigA0A1;
        double sigA1A0;
  	double sigA1;
  	double sigB0;
	double sigB0B1;
        double sigB1B0;
  	double sigB1;
  	TrackCov();
  	TrackCov(double siga0, double siga0a1,double siga1a0, double siga1, double sigb0,double sigb0b1, double sigb1b0, double sigb1) : sigA0(siga0), sigA0A1(siga0a1), sigA1A0(siga1a0), sigA1(siga1), sigB0(sigb0), sigB0B1(sigb0b1), sigB1B0(sigb1b0), sigB1(sigb1) {};
   
   };
   //Struct To Hold Track Parameters
   struct TrackParams{
  	double A0;
  	double A1;
  	double B0;
  	double B1;
  	double T0;
	double deltaA0;
  	double deltaA1;
  	double deltaB0;
  	double deltaB1;
  	double deltaT0;
  	TrackParams();
  	TrackParams(double a0, double a1, double  b0, double b1) : A0(a0), A1(a1), B0(b0) , B1(b1){};
	TrackCov Covarience;
	
   };
   //Struct to hold Coordinate System
   struct TrackAxes{
  	XYZVec _XDoublePrime;
  	XYZVec _YDoublePrime;
  	XYZVec _ZPrime;
  	TrackAxes();
  	TrackAxes(XYZVec X, XYZVec Y, XYZVec Z) : _XDoublePrime(X),_YDoublePrime(Y),_ZPrime(Z){};
  	
   };
   //Struct to store a Track Equation (optional)
   struct TrackEquation{
   	XYZVec Pos;
        XYZVec Dir;
        TrackEquation();
        TrackEquation(XYZVec P, XYZVec D) : Pos(P), Dir(D){};
   };
   //Struct to store Diagnostics associated with seed fit:
   struct TrackSeedDiag{
   	double FinalChiX;
   	double FinalChiY;
   	double FinalChiTot;
   	
   	double InitialChiX;
   	double InitialChiY;
   	double InitialChiTot;
   	
	TrackSeedDiag();

   };

namespace mu2e {
  
  class CosmicTrack{
  public:
  
    
    // Constructors
    public:
	    CosmicTrack();
	    CosmicTrack(TrackParams params) : FitParams(params) {};
	    ~CosmicTrack();
	    
	 //---------------Accessors:--------------//
	  
	    double get_fit_phi() const { return FitPhi;}
	    double get_fit_theta() const { return FitTheta;}
	    int get_N() const { return _Nhits;}
	    int get_iter(){return _niters;}
	    
	    TrackParams GetFitParams()const{ 
	    	return FitParams;
	    }
	    TrackParams GetInitParams()const{ 
	    	return InitParams;
	    }

	    XYZVec GetTrackDirection() const{ //TODO use private direction 
	    	XYZVec Direction(FitParams.A1, FitParams.B1, 1);
	    	return Direction;
	    }
	    
	    XYZVec GetInitTrackDirection() const{
	    	XYZVec Direction(InitParams.A1, InitParams.B1, 1);
	    	return Direction.Unit();
	    }
	    
	    XYZVec GetTrackPosition() const{ //TODO use private positon?
	    	XYZVec Position(FitParams.A0, FitParams.B0, 0);
	    	return Position;
	    }

 	     XYZVec GetPOCA() const{
		return POCA;
             }

	     double GetDOCA() const{
		return DOCA;
             }

	     void clear_errors(); //clears track info and diagnostics
	     void clear_parameters(); //clears just track info
	     void clear_diag();
	     void clear_cov();
	//-------------Modiffiers of Track Parameters/Features---------------//
	    void SetFitParams(TrackParams par){ 
	    	this->FitParams = par;
	    }
	    void SetInitParams(TrackParams par){ 
	    	this->InitParams = par;
	    }

	    void SetMinuitParams(double par_a0, double par_a1, double par_b0, double par_b1, double par_t0 ){
	    	this->MinuitFitParams.A0 = par_a0;
	 	this->MinuitFitParams.A1 = par_a1;
	 	this->MinuitFitParams.B0 = par_b0;
	 	this->MinuitFitParams.B1 = par_b1;
	  	this->MinuitFitParams.T0 = par_t0;
	    
	    }
	    void SetFitTrackCoOrdSystem(TrackAxes coordsys){
	    	this->TrackCoordSystem = coordsys;
	    }

	    void SetInitTrackCoordSystem(TrackAxes coordsys){
	    	this->InitTrackCoordSystem = coordsys;
    	    }
    	    
    	    void SetCovarience(double siga0, double siga0a1, double siga1a0, double siga1, double sigb0, double sigb0b1, double sigb1b0, double sigb1){ 
    	    	TrackCov Cov(siga0, siga0a1, siga1a0, siga1, sigb0, sigb0b1, sigb1b0, sigb1);
	    	this->FitParams.Covarience = Cov;
	    }
	     void SetMinuitCoordSystem(TrackAxes coordsys){
	    	this->MinuitCoordSystem = coordsys;
	    }
	    void Set_N(unsigned N){_Nhits = N;}
	    
	    void SetTrackEquation(TrackEquation Track){ 
	        this->FitEquation = Track;
	     }

	    void SetTrackEquationXYZ(TrackEquation Track){ 
	        this->FitEquationXYZ = Track;
	     }

	    void SetMinuitTrackEquation(TrackEquation Track){ 
	        this->MinuitFitEquation = Track;
	     }
	    
	    void SetTrackDirection(XYZVec dir){
	    	this->Direction = dir;
	    }

	    void SetTrackPosition(XYZVec pos){
	    	this->Position = pos;
	    }
	   void SetFirstHitVec(double x, double y, double z) {
		this->FirstHitVec.SetXYZ(x,y,z);

	    }
  	    void SetLastHitVec(double x, double y, double z) {
		this->LastHitVec.SetXYZ(x,y,z);
	    }

	    void SetTrackPOCA() {
		XYZVec TrackerCenter(0,0,0);//TODO 🙄
        	PointLinePCA_XYZ PCA = PointLinePCA_XYZ(TrackerCenter, this->FirstHitVec, this->LastHitVec);
		this->POCA = PCA.pca();
		this->DOCA = PCA.dca();
	    }
  
	    void set_fit_theta(double track_angle){ FitTheta = track_angle;}
	    void set_fit_phi(double track_angle){ FitPhi = track_angle;}

	    //-----------Fill Diag info this is legacy it will probably be removed soon----------//
	    void set_finalchisq_dof(double finalchisq_dof)   { Diag.FinalChiTot = finalchisq_dof; }
	    void set_finalchisq_dofX(double finalchisq_dofX) { Diag.FinalChiX = finalchisq_dofX; }
	    void set_finalchisq_dofY(double finalchisq_dofY) { Diag.FinalChiY = finalchisq_dofY; }
	    
	    void set_initchisq_dof(double initchisq_dof)   { Diag.InitialChiTot = initchisq_dof; }
	    void set_initchisq_dofX(double initchisq_dofX) { Diag.InitialChiX = initchisq_dofX; }
	    void set_initchisq_dofY(double initchisq_dofY) { Diag.InitialChiY = initchisq_dofY; }
	    void set_niter(int iter){ _niters= (iter); }
            //----------------------------End Diag Fill---------------------------------//
    	     TrackParams FitParams; //From Iterative Seed Fit (not drift)
             TrackParams InitParams; //Initial fit (first estimate)
	     TrackParams MinuitFitParams; // Minuit Params (seed drift final)

	     TrackAxes MinuitCoordSystem;//Result from Minuit Fit
	     TrackAxes TrackCoordSystem;//Seed Fit Result Axes
	     TrackAxes InitTrackCoordSystem;//Initial Axes (start->end line)
	     
	     TrackEquation FitEquation;
	     TrackEquation FitEquationXYZ;
	     TrackEquation MinuitFitEquation;

	     TrackSeedDiag Diag;
	     
	     XYZVec sigmaPos;
	     XYZVec sigmaDir;
	     
	     bool converged = false;
	     bool minuit_converged = false;
	     unsigned n_outliers = 0;

	     double FitPhi;
	     double FitTheta;
	     XYZVec POCA;
	     double DOCA;

  private:
	    
	    unsigned _Nhits;
	   
	    int _niters;
	    XYZVec Direction;
	    XYZVec Position;
	    XYZVec FirstHitVec;
	    XYZVec LastHitVec; 
	    
  };
  
}

#endif

