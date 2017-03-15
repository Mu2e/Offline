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
	      crateRadiusIn_(0),crateRadiusOut_(0),crateHalfLength_(0),crateHalfX_(0),crateHalfY_(0),crateHalfZ_(0),
	      crateTopHalfY_(0),crateBottomHalfY_(0),crateSideHalfY_(0),crateSideHalfX_(0),shieldHalfY_(0),shieldDispZ_(0),
	      shieldDZFront_(0),stepsRadiusIn_(0),stepsRadiusOut_(0),outerRingEdgeDepth_(0),
	      outerRingEdgeThickness_(0),caseThickness_(0),caseThicknessIn_(0),
	      caseThicknessOut_(0),envelopeInRadius_(0),envelopeOutRadius_(0),
	      envelopeZ0_(0),envelopeZ1_(0),boardHalfY_(0),radiatorHalfY_(0),activeStripHalfY_(0),passiveStripHalfY_(0),nPipes_(0),pipeRadius_(0),
	      pipeThickness_(0),pipeTorRadius_()
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
           void crateToDiskDeltaZ(double value)    {crateToDiskDeltaZ_ = value;}
           void crateHalfLength(double value)      {crateHalfLength_ = value;}
           void crateHalfX(double value)           {crateHalfX_ = value;}
           void crateHalfY(double value)           {crateHalfY_ = value;}
           void crateHalfZ(double value)           {crateHalfZ_ = value;}
           void crateTopHalfY(double value)        {crateTopHalfY_ = value;}
           void crateBottomHalfY(double value)     {crateBottomHalfY_ = value;}
           void crateSideHalfY(double value)       {crateSideHalfY_ = value;}
           void crateSideHalfX(double value)       {crateSideHalfX_ = value;}
           void shieldHalfY(double value)          {shieldHalfY_    = value;}
           void shieldDispZ  (double value)        {shieldDispZ_    = value;}
           void shieldDZFront(double value)        {shieldDZFront_    = value;}
            
           void stepsRadiusIn(double value)          {stepsRadiusIn_ = value;}
           void stepsRadiusOut(double value)         {stepsRadiusOut_ = value;}
           void outerRingEdgeDepth(double value)     {outerRingEdgeDepth_ = value;}
           void outerRingEdgeThickness(double value) {outerRingEdgeThickness_ = value;}
	   void caseThickness(double value)          {caseThickness_ = value;}
           void caseThicknessIn(double value)        {caseThicknessIn_ = value;}
           void caseThicknessOut(double value)       {caseThicknessOut_ = value;}
           void envelopeInRadius(double value)       {envelopeInRadius_ = value;}
           void envelopeOutRadius(double value)      {envelopeOutRadius_ = value;}
           void envelopeZ0(double value)             {envelopeZ0_ = value;}
           void envelopeZ1(double value)             {envelopeZ1_ = value;} 
           void boardHalfY(double value)             {boardHalfY_ = value;}
           void radiatorHalfY(double value)          {radiatorHalfY_ = value;}
           void activeStripHalfY(double value)       {activeStripHalfY_ = value;}
           void passiveStripHalfY(double value)       {passiveStripHalfY_ = value;}
	   void refractiveIndex(double value)        {refractiveIndex_ = value;}
	   void crystalDecayTime(double value)       {crystalDecayTime_ = value;}

           
	   
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
           double crateToDiskDeltaZ()   const      {return crateToDiskDeltaZ_;}
           double crateHalfX()          const      {return crateHalfX_;}
	   double crateHalfY()          const      {return crateHalfY_;}
	   double crateHalfZ()          const      {return crateHalfZ_;}
	   double crateTopHalfY()       const      {return crateTopHalfY_;}
           double crateBottomHalfY()    const      {return crateBottomHalfY_;}
           double crateSideHalfY()      const      {return crateSideHalfY_;}
           double crateSideHalfX()      const      {return crateSideHalfX_;}         
           double shieldHalfY()         const      {return shieldHalfY_;}         
           double shieldDispZ  ()         const      {return shieldDispZ_;}         
           double shieldDZFront()         const      {return shieldDZFront_;}         
      
         
           double stepsRadiusIn()          const      {return stepsRadiusIn_;}
           double stepsRadiusOut()         const      {return stepsRadiusOut_;}
           double outerRingEdgeDepth()     const      {return outerRingEdgeDepth_;}
           double outerRingEdgeThickness() const      {return outerRingEdgeThickness_;}
	   double caseThickness()          const      {return caseThickness_;}
           double caseThicknessIn()        const      {return caseThicknessIn_;}
           double caseThicknessOut()       const      {return caseThicknessOut_;}
           double envelopeInRadius()       const      {return envelopeInRadius_;}
           double envelopeOutRadius()      const      {return envelopeOutRadius_;}
           double envelopeZ0()             const      {return envelopeZ0_;}
           double envelopeZ1()             const      {return envelopeZ1_;}
           double boardHalfY()             const      {return boardHalfY_;}
           double radiatorHalfY()          const      {return radiatorHalfY_;}
           double activeStripHalfY()       const      {return activeStripHalfY_;}
           double passiveStripHalfY()       const      {return passiveStripHalfY_;}
	   double refractiveIndex()        const      {return refractiveIndex_; }
	   double crystalDecayTime()       const      {return crystalDecayTime_; }

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
          double crateToDiskDeltaZ_;
          double crateHalfX_;
          double crateHalfY_;
          double crateHalfZ_;
          double crateTopHalfY_;
          double crateBottomHalfY_;
          double crateSideHalfY_;
          double crateSideHalfX_;
          double shieldHalfY_;
          double shieldDispZ_;  
          double shieldDZFront_;

          double stepsRadiusIn_;
	  double stepsRadiusOut_;
          double outerRingEdgeDepth_;
          double outerRingEdgeThickness_;
          double caseThickness_;
          double caseThicknessIn_;
          double caseThicknessOut_;
          double envelopeInRadius_;
          double envelopeOutRadius_;
          double envelopeZ0_;
          double envelopeZ1_;
          double boardHalfY_;
          double radiatorHalfY_;
          double activeStripHalfY_;
          double passiveStripHalfY_;
	   
          double refractiveIndex_;
          double crystalDecayTime_;

	  unsigned int         nPipes_;
	  double               pipeRadius_;
	  double               pipeThickness_;
	  std::vector<double>  pipeTorRadius_;
     };

}    

#endif 
