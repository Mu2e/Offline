#ifndef CalorimeterGeom_CaloInfo_hh
#define CalorimeterGeom_CaloInfo_hh
//
// Contains data for the disk calorimeter class
//
// Original author B. Echenard
//


#include <vector>

namespace mu2e {

    class CaloInfo {

       public:

           CaloInfo(): 
	      crystalShift_(0),crystalHalfTrans_(0),crystalHalfLength_(0),crystalVolume_(0),wrapperThickness_(0),
	      nROPerCrystal_(0),roHalfTrans_(0),roHalfThickness_(0),roElecHalfX_(0),roElecHalfY_(0),roElecHalfZ_(0),
	      crateRadiusIn_(0),crateRadiusOut_(0),crateHalfLength_(0),
	      caseThickness_(0),envelopeInRadius_(0),envelopeOutRadius_(0),envelopeZ0_(0),envelopeZ1_(0),
	      nPipes_(0),pipeRadius_(0),pipeThickness_(0),pipeTorRadius_()   	      
	   {}
	     
           ~CaloInfo() {}

           void crystalShift(bool value)           {crystalShift_ = value;}
           void crystalHalfLength(double value)    {crystalHalfLength_ = value;}
           void crystalHalfTrans(double value)     {crystalHalfTrans_ = value;}
           void crystalVolume(double value)        {crystalVolume_ = value;}
	   void wrapperThickness(double value)     {wrapperThickness_ = value;}
           
	   void nROPerCrystal(int value)           {nROPerCrystal_ = value;}
           void roHalfTrans(double value)          {roHalfTrans_ = value;}
           void roHalfThickness(double value)      {roHalfThickness_ = value;}
           void roElecHalfX(double value)          {roElecHalfX_ = value;}
           void roElecHalfY(double value)          {roElecHalfY_ = value;}
           void roElecHalfZ(double value)          {roElecHalfZ_ = value;}
           void crateRadiusIn(double value)        {crateRadiusIn_ = value;}
           void crateRadiusOut(double value)       {crateRadiusOut_ = value;}
           void crateHalfLength(double value)      {crateHalfLength_ = value;}
          
	   void caseThickness(double value)        {caseThickness_ = value;}
           void envelopeInRadius(double value)     {envelopeInRadius_ = value;}
           void envelopeOutRadius(double value)    {envelopeOutRadius_ = value;}
           void envelopeZ0(double value)           {envelopeZ0_ = value;}
           void envelopeZ1(double value)           {envelopeZ1_ = value;} 
	   void refractiveIndex(double value)      {refractiveIndex_ = value;}
	   void crystalDecayTime(double value)     {crystalDecayTime_ = value;}

           
	   
           bool   crystalShift()        const      {return crystalShift_;}
           double crystalHalfLength()   const      {return crystalHalfLength_;}
           double crystalHalfTrans()    const      {return crystalHalfTrans_;}
           double crystalVolume()       const      {return crystalVolume_;}
           double wrapperThickness()    const      {return wrapperThickness_;}
           
	   int    nROPerCrystal()       const      {return nROPerCrystal_;}
           double roHalfTrans()         const      {return roHalfTrans_;}
           double roHalfThickness()     const      {return roHalfThickness_;}
           double roElecHalfX()         const      {return roElecHalfX_;}
           double roElecHalfY()         const      {return roElecHalfY_;}
           double roElecHalfZ()         const      {return roElecHalfZ_;}
           double crateRadiusIn()       const      {return crateRadiusIn_;}
           double crateRadiusOut()      const      {return crateRadiusOut_;}
           double crateHalfLength()     const      {return crateHalfLength_;}
           
	   double caseThickness()       const      {return caseThickness_;}
           double envelopeInRadius()    const      {return envelopeInRadius_;}
           double envelopeOutRadius()   const      {return envelopeOutRadius_;}
           double envelopeZ0()          const      {return envelopeZ0_;}
           double envelopeZ1()          const      {return envelopeZ1_;}
	   double refractiveIndex()     const      {return refractiveIndex_; }
	   double crystalDecayTime()    const      {return crystalDecayTime_; }




           void nPipes(int value)                         {nPipes_ = value;}
           void pipeRadius(double value)                  {pipeRadius_ = value;}
           void pipeThickness(double value)               {pipeThickness_ = value;}           
	   void pipeTorRadius(std::vector<double>& value) {pipeTorRadius_ = value;}

           int    nPipes()                            const {return nPipes_;}
           double pipeRadius()                        const {return pipeRadius_;}
           double pipeThickness()                     const {return pipeThickness_;}
           const std::vector<double>& pipeTorRadius() const {return pipeTorRadius_;}
           
           
           
       private:

          bool   crystalShift_;
	  double crystalHalfTrans_;
	  double crystalHalfLength_;
	  double crystalVolume_;
          double wrapperThickness_;

          int    nROPerCrystal_;
          double roHalfTrans_;
          double roHalfThickness_;
          double roElecHalfX_;
          double roElecHalfY_;
          double roElecHalfZ_;
	  
	  double crateRadiusIn_;
	  double crateRadiusOut_;
          double crateHalfLength_;

          double caseThickness_;
          double envelopeInRadius_;
          double envelopeOutRadius_;
          double envelopeZ0_;
          double envelopeZ1_;
	   
          double refractiveIndex_;
          double crystalDecayTime_;

	  unsigned int         nPipes_;
	  double               pipeRadius_;
	  double               pipeThickness_;
	  std::vector<double>  pipeTorRadius_;
     };

}    

#endif 
