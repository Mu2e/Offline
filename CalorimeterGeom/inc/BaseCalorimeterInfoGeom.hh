#ifndef CalorimeterGeom_BaseCalorimeterInfoGeom_hh
#define CalorimeterGeom_BaseCalorimeterInfoGeom_hh
//
// $Id: BaseCalorimeterInfoGeom.hh,v 1.1 2014/08/01 21:49:38 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 21:49:38 $
//
// Contains data for the BaseCalorimeter class, much of it will be moved into the calibration 
// in the future
//
// Might consider to have a class for the calibration system and pipes later (or might not...)
//
// Original author B. Echenard
//


#include <vector>


namespace mu2e {

    

    class BaseCalorimeterInfoGeom {


       public:

           BaseCalorimeterInfoGeom()  {}
           virtual ~BaseCalorimeterInfoGeom() {}

           void crystalNedges(int value)           {_crystalNumEdges = value;}
           void crystalShift(bool value)           {_crystalShift = value;}
           void nROPerCrystal(int value)           {_nROPerCrystal = value;}
           void crystalHalfLength(double value)    {_crystalHalfLength = value;}
           void crystalHalfTrans(double value)     {_crystalHalfTrans = value;}
           void crystalVolume(double value)        {_crystalVolume = value;}
	   void wrapperThickness(double value)     {_wrapperThickness = value;}
           void caseThickness(double value)        {_caseThickness = value;}
           void roHalfTrans(double value)          {_roHalfTrans = value;}
           void roHalfThickness(double value)      {_roHalfThickness = value;}
           void envelopeInRadius(double value)    {_envelopeInRadius = value;}
           void envelopeOutRadius(double value)   {_envelopeOutRadius = value;}
           void envelopeZ0(double value)          {_envelopeZ0 = value;}
           void envelopeZ1(double value)          {_envelopeZ1 = value;} 

           int    crystalNedges()       const      {return _crystalNumEdges;}
           bool   crystalShift()        const      {return _crystalShift;}
           int    nROPerCrystal()       const      {return _nROPerCrystal;}
           double crystalHalfLength()   const      {return _crystalHalfLength;}
           double crystalHalfTrans()    const      {return _crystalHalfTrans;}
           double crystalVolume()       const      {return _crystalVolume;}
           double wrapperThickness()    const      {return _wrapperThickness;}
           double caseThickness()       const      {return _caseThickness;}
           double roHalfTrans()         const      {return _roHalfTrans;}
           double roHalfThickness()     const      {return _roHalfThickness;}
           double envelopeInRadius()   const      {return _envelopeInRadius;}
           double envelopeOutRadius()  const      {return _envelopeOutRadius;}
           double envelopeZ0()         const      {return _envelopeZ0;}
           double envelopeZ1()         const      {return _envelopeZ1;}
           
	   
	   
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

          int    _crystalNumEdges;
          bool   _crystalShift;
          int    _nROPerCrystal;
	  double _crystalHalfTrans;
	  double _crystalHalfLength;
	  double _crystalVolume;
          double _wrapperThickness;
          double _roHalfTrans;
          double _roHalfThickness;
          double _caseThickness;

          double _envelopeInRadius;
          double _envelopeOutRadius;
          double _envelopeZ0;
          double _envelopeZ1;
	   
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

#endif 
