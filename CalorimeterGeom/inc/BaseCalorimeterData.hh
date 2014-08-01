#ifndef CalorimeterGeom_BaseCalorimeterData_hh
#define CalorimeterGeom_BaseCalorimeterData_hh
//
// $Id: BaseCalorimeterData.hh,v 1.1 2014/08/01 21:04:58 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 21:04:58 $
//
// Contains data for the BaseCalorimeter class
// Should go away in the future when there will be the condition database
//
// Might consider to have a class for the calibration system and pipes later 
//
// Original author B. Echenard
//

//C++ includes
#include <vector>



namespace mu2e {

    

    class BaseCalorimeterData {


       public:

           BaseCalorimeterData()  {}
           virtual ~BaseCalorimeterData() {}

           void nROPerCrystal(int value)           {_nROPerCrystal = value;}
           void crystalHalfLength(double value)    {_crystalHalfLength = value;}
           void crystalHalfTrans(double value)     {_crystalHalfTrans = value;}
           void crystalVolume(double value)        {_crystalVolume = value;}
	   void wrapperThickness(double value)     {_wrapperThickness = value;}
           void shellThickness(double value)       {_shellThickness = value;}
           void caseThickness(double value)        {_caseThickness = value;}
           void roHalfTrans(double value)          {_roHalfTrans = value;}
           void roHalfThickness(double value)      {_roHalfThickness = value;}
           void enveloppeInRadius(double value)    {_enveloppeInRadius = value;}
           void enveloppeOutRadius(double value)   {_enveloppeOutRadius = value;}
           void enveloppeZ0(double value)          {_enveloppeZ0 = value;}
           void enveloppeZ1(double value)          {_enveloppeZ1 = value;} 

           int    nROPerCrystal()       const      {return _nROPerCrystal;}
           double crystalHalfLength()   const      {return _crystalHalfLength;}
           double crystalHalfTrans()    const      {return _crystalHalfTrans;}
           double crystalVolume()       const      {return _crystalVolume;}
           double wrapperThickness()    const      {return _wrapperThickness;}
           double shellThickness()      const      {return _shellThickness;}
           double caseThickness()       const      {return _caseThickness;}
           double roHalfTrans()         const      {return _roHalfTrans;}
           double roHalfThickness()     const      {return _roHalfThickness;}
           double enveloppeInRadius()   const      {return _enveloppeInRadius;}
           double enveloppeOutRadius()  const      {return _enveloppeOutRadius;}
           double enveloppeZ0()         const      {return _enveloppeZ0;}
           double enveloppeZ1()         const      {return _enveloppeZ1;}
           
	   
	   
	   
	   
	   
	   void nonUniformity(double value)        {_nonUniformity = value;}
	   void timeGap(double value)              {_timeGap = value;}
	   void electronEdep(double value)         {_electronEdep = value;}
	   void electronEmin(double value)         {_electronEmin = value;}
	   void apdMeanNoise(double value)         {_apdMeanNoise = value;}
	   void apdSigmaNoise(double value)        {_apdSigmaNoise = value;}
	   void lysoLightYield(double value)       {_lysoLightYield = value;}
	   void apdQuantumEff(double value)        {_apdQuantumEff = value;}
	   void apdCollectEff(double value)        {_apdCollectEff = value;}
           
	   double nonUniformity()       const      {return _nonUniformity; }
	   double timeGap()             const      {return _timeGap; }
	   double electronEdep()        const      {return _electronEdep; }
	   double electronEmin()        const      {return _electronEmin; }
	   double apdMeanNoise()        const      {return _apdMeanNoise;}
	   double apdSigmaNoise()       const      {return _apdSigmaNoise;}
	   double lysoLightYield()      const      {return _lysoLightYield;}
	   double apdQuantumEff()       const      {return _apdQuantumEff;}
	   double apdCollectEff()       const      {return _apdCollectEff;}



           void nPipes(int value)                  {_nPipes = value;}
           void pipeRadius(double value)           {_pipeRadius = value;}
           void pipeThickness(double value)        {_pipeThickness = value;}           
	   void pipeTorRadius(std::vector<double>& value) { _pipeTorRadius = value;}

           int    nPipes()             const       {return _nPipes;}
           double pipeRadius()         const       {return _pipeRadius;}
           double pipeThickness()      const       {return _pipeThickness;}
           std::vector<double> const& pipeTorRadius() const  {return _pipeTorRadius;}



       private:

          int    _nROPerCrystal;
	  double _crystalHalfTrans;
	  double _crystalHalfLength;
	  double _crystalVolume;
          double _wrapperThickness;
          double _roHalfTrans;
          double _roHalfThickness;
          double _shellThickness;
          double _caseThickness;

          double _enveloppeInRadius;
          double _enveloppeOutRadius;
          double _enveloppeZ0;
          double _enveloppeZ1;
	   
          double _nonUniformity;
          double _timeGap;
          double _electronEdep; // energy deposition of charged particle crossing APD
          double _electronEmin; // minimum energy deposition of charged particle crossing APD

          double _apdMeanNoise; //MeV
          double _apdSigmaNoise;//MeV

          double _lysoLightYield;
          double _apdQuantumEff;//quantum efficiency for Hamamatsu S8664-1010 for a radiation wavelenght of 402nm (typical of lyso)
          double _apdCollectEff;//light collection efficiency for  30 mm2 area of the crystal efficiency for Hamamatsu S8664-1010

	  unsigned int         _nPipes;
	  double               _pipeRadius;
	  double               _pipeThickness;
	  std::vector<double>  _pipeTorRadius;


     };

}    

#endif /* CalorimeterGeom_BaseCalorimeterData_hh*/
