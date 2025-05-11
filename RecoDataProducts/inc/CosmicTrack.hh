#ifndef RecoDataProducts_CosmicTrack_hh
#define RecoDataProducts_CosmicTrack_hh
////S. Middleton, Feb 2019 - Cosmic track class, main purpose id to store diagnostics.
#include "KinKal/General/Vectors.hh"
#include <vector>

#include "Offline/DataProducts/inc/GenVector.hh"


using CLHEP::Hep3Vector;
namespace mu2e {
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
          TrackCov Covarience; //FIXME backwards compatibility
    std::vector<double> cov;

          XYZVectorF Direction() const { return XYZVectorF(A1, B1, 1).unit();};
          XYZVectorF Position() const { return XYZVectorF(A0, B0, 0);};

   };

   //Struct to hold Coordinate System
   struct TrackAxes{
          XYZVectorF _XDoublePrime;
          XYZVectorF _YDoublePrime;
          XYZVectorF _ZPrime;
          TrackAxes();
          TrackAxes(XYZVectorF const& X, XYZVectorF const& Y, XYZVectorF const& Z) : _XDoublePrime(X),_YDoublePrime(Y),_ZPrime(Z){};
   };

   //Struct to store a Track Equation (optional)
   struct TrackEquation{
           XYZVectorF Pos;
        XYZVectorF Dir;
        TrackEquation();
        TrackEquation(XYZVectorF const& P, XYZVectorF const& D) : Pos(P), Dir(D){};
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

  class CosmicTrack{

    // Constructors
    public:
      CosmicTrack();
      CosmicTrack(TrackParams params) : FitParams(params) {};

      //---------------Accessors:--------------//

      int get_iter(){return niters;}

      TrackParams GetFitParams()const{
        return FitParams;
      }
      TrackParams GetInitParams()const{
        return InitParams;
      }

      //-------------Modiffiers of Track Parameters/Features---------------//
      void SetFitParams(TrackParams par){
        this->FitParams = par;
      }
      void SetInitParams(TrackParams par){
        this->InitParams = par;
      }

      void SetMinuitParams(double par_a0, double par_a1, double par_b0, double par_b1, double par_t0_ ){
        this->MinuitParams.A0 = par_a0;
        this->MinuitParams.A1 = par_a1;
        this->MinuitParams.B0 = par_b0;
        this->MinuitParams.B1 = par_b1;
        this->MinuitParams.T0 = par_t0_;
      }

      void SetFitCoordSystem(TrackAxes coordsys){
        this->FitCoordSystem = coordsys;
      }

      void SetInitCoordSystem(TrackAxes coordsys){
        this->InitCoordSystem = coordsys;
      }

      void SetCovarience(double siga0, double siga0a1, double siga1a0, double siga1, double sigb0, double sigb0b1, double sigb1b0, double sigb1){
      TrackCov Cov(siga0, siga0a1, siga1a0, siga1, sigb0, sigb0b1, sigb1b0, sigb1);
        this->FitParams.Covarience = Cov;
      }

      void SetMinuitCoordSystem(TrackAxes coordsys){
        this->MinuitCoordSystem = coordsys;
      }

      void SetFitEquation(TrackEquation Track){
        this->FitEquation = Track;
      }

      void SetMinuitEquation(TrackEquation Track){
        this->MinuitEquation = Track;
      }

      XYZVectorF Pos0() const {
        XYZVectorF intercept(this->MinuitParams.A0, 0, this->MinuitParams.B0);
        return intercept;
      }

      XYZVectorF Dir() const {
        XYZVectorF direction(this->MinuitParams.A1, -1, this->MinuitParams.B1);
        return direction.Unit();
      }

      //----------------------------End Diag Fill---------------------------------//
      TrackParams InitParams; //Initial fit (first estimate)
      TrackParams FitParams; //From Iterative Seed Fit (not drift)
      TrackParams MinuitParams; // Minuit Params (seed drift final)

      TrackAxes InitCoordSystem; //Initial Axes (start->end line)
      TrackAxes FitCoordSystem; //Seed Fit Result Axes
      TrackAxes MinuitCoordSystem; //Result from Minuit Fit

      TrackEquation InitEquation;
      TrackEquation FitEquation; // pos and dir in detector coordinates
      TrackEquation MinuitEquation; // pos and dir in detector coordinates

      TrackSeedDiag Diag;

      bool converged = false;
      bool minuit_converged = false;

      //Fill Diag info this is legacy it will probably be removed soon: - TODO can we remove this?
      void set_finalchisq_dof(double finalchisq_dof)   { Diag.FinalChiTot = finalchisq_dof; }
      void set_finalchisq_dofX(double finalchisq_dofX) { Diag.FinalChiX = finalchisq_dofX; }
      void set_finalchisq_dofY(double finalchisq_dofY) { Diag.FinalChiY = finalchisq_dofY; }

      void set_initchisq_dof(double initchisq_dof)   { Diag.InitialChiTot = initchisq_dof; }
      void set_initchisq_dofX(double initchisq_dofX) { Diag.InitialChiX = initchisq_dofX; }
      void set_initchisq_dofY(double initchisq_dofY) { Diag.InitialChiY = initchisq_dofY; }
      void set_niter(int iter){ niters= (iter);}

      void SetKinKalParams(double d0, double phi0,double z0, double cost, double t0, double mom){
        d0_  = d0;
        phi0_ = phi0;
        z0_ = z0;
        cost_ = cost;
        t0_ = t0;
        mom_ = mom;
      }

      //Kinkal params:
      double d0() const { return d0_; }
      double z0() const { return z0_; }
      double phi0() const { return phi0_; }
      double cost() const { return cost_; }
      double t0() const { return t0_; }
      double mom() const { return mom_; }

    private:
      int niters;
      //KinematicLine parameters:
      double d0_;
      double z0_;
      double t0_;
      double cost_;
      double  phi0_;
      double mom_;

  };
}

#endif

