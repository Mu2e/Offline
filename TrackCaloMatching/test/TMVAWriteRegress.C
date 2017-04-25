#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;


   
void TMVAWriteRegress(int choice, TString trainDS, int mode=1, bool noSpectator=0) 
{
   
   if (choice!=1 && choice !=2&& choice !=3) {cout<<"Choice must be 1 or 2 or 3"<<endl; return;}
   if (mode!=1 && mode !=2 && mode!=3 ) {cout<<"Mose must be 1 or 2 or 3"<<endl; return;}
     

   TMVA::Tools::Instance();

   // Create a new root output file
   TString outfileName( "TMVAReg.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   TMVA::Factory *factory = new TMVA::Factory("TMVARegression", outputFile,"!V:!Silent:Color:DrawProgressBar:analysisType=Regression" );
   
   TString dlname="weightsX";
   if (choice==2) dlname="weightsY";
   if (choice==3) dlname="weightsA";
   TMVA::DataLoader *dataloader=new TMVA::DataLoader(dlname);

   // Define the input variables that shall be used for the MVA training
   if (mode==3) dataloader->AddVariable( "e0",  "Energy cell 0",  "MeV", 'F' );
   if (mode==3) dataloader->AddVariable( "e1",  "Energy cell 1",  "MeV", 'F' );   
   dataloader->AddVariable( "e2",  "Energy cell 2",  "MeV", 'F' );
   dataloader->AddVariable( "e3",  "Energy cell 3",  "MeV", 'F' );
   dataloader->AddVariable( "e4",  "Energy cell 4",  "MeV", 'F' );
   dataloader->AddVariable( "e5",  "Energy cell 5",  "MeV", 'F' );
   dataloader->AddVariable( "e6",  "Energy cell 6",  "MeV", 'F' );
   dataloader->AddVariable( "e7",  "Energy cell 7",  "MeV", 'F' );
   dataloader->AddVariable( "e8",  "Energy cell 8",  "MeV", 'F' );
   dataloader->AddVariable( "e9",  "Energy cell 9",  "MeV", 'F' );
   dataloader->AddVariable( "e10", "Energy cell 10", "MeV", 'F' );
   dataloader->AddVariable( "e11", "Energy cell 11", "MeV", 'F' );
   dataloader->AddVariable( "e12", "Energy cell 12", "MeV", 'F' );
   dataloader->AddVariable( "e13", "Energy cell 13", "MeV", 'F' );
   dataloader->AddVariable( "e14", "Energy cell 14", "MeV", 'F' );
   dataloader->AddVariable( "e15", "Energy cell 15", "MeV", 'F' );
   dataloader->AddVariable( "e16", "Energy cell 16", "MeV", 'F' );
   dataloader->AddVariable( "e17", "Energy cell 17", "MeV", 'F' );
   dataloader->AddVariable( "e18", "Energy cell 18", "MeV", 'F' );
   dataloader->AddVariable( "e19", "Energy cell 19", "MeV", 'F' );
   if (mode!=3)dataloader->AddVariable( "e20", "Energy cell 20", "MeV", 'F' );

   dataloader->AddVariable( "t0","vdPhi", "", 'F' );
   dataloader->AddVariable( "t1","vdTheta", "", 'F' );
   dataloader->AddVariable( "r0","r0", "", 'F' );
   dataloader->AddVariable( "z0","z0", "", 'F' );
   //dataloader->AddVariable( "ip0","ip0", "", 'F' );
   //dataloader->AddVariable( "ip1","ip1", "", 'F' );


   if (!noSpectator)
   {
       // You can add so-called "Spectator variables", which are not used in the MVA training, 
       dataloader->AddSpectator( "c0",  "c0", "mm", 'F' );
       dataloader->AddSpectator( "c1",  "c1", "mm", 'F' );
       dataloader->AddSpectator( "c2",  "c2", "mm", 'F' );
       dataloader->AddSpectator( "c3",  "c3", "mm", 'F' );
       dataloader->AddSpectator( "c4",  "c4", "mm", 'F' );
       dataloader->AddSpectator( "c5",  "c5", "mm", 'F' );
       dataloader->AddSpectator( "ffx", "ffx", "mm", 'F' );
       dataloader->AddSpectator( "ffy", "ffy", "mm", 'F' );
       dataloader->AddSpectator( "ffz", "ffz", "mm", 'F' );
       dataloader->AddSpectator( "vdx", "vdx", "mm", 'F' );
       dataloader->AddSpectator( "vdy", "vdy", "mm", 'F' );
       dataloader->AddSpectator( "vdz", "vdz", "mm", 'F' );
   }

   // Add the variable carrying the regression target
   if (choice==1) dataloader->AddTarget( "targetX" ); 
   if (choice==2) dataloader->AddTarget( "targetY" ); 
   if (choice==3) dataloader->AddTarget( "targetA" ); 



   // --- Register the regression tree
   TFile input(trainDS);
   TTree *regTree = (TTree*)input.Get("TreeReg");


   // You can add an arbitrary number of regression trees
   dataloader->AddRegressionTree( regTree, 1.0 );


   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycut = "";   //"target>0&&target<2"; 


   // tell the factory to use all remaining events in the trees after training for testing:
   dataloader->PrepareTrainingAndTestTree( mycut,  "nTrain_Regression=0:nTest_Regression=0:SplitMode=Random:NormMode=NumEvents:!V" );


   //linear regression
   //factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "!H:!V:VarTransform=None" );
 
   //MLP
   //factory->BookMethod( dataloader,  TMVA::Types::kMLP, "MLP","!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=100:HiddenLayers=N+10:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator" );

   //svm
   //factory->BookMethod( dataloader,  TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

   
   // Boosted Decision Trees
   //factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT","!H:!V:NTrees=100:MinNodeSize=1.0%:BoostType=AdaBoostR2:SeparationType=RegressionVariance:nCuts=20:PruneMethod=CostComplexity:PruneStrength=30" );

   //decorrelated BDT
   factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG","!H:!V:NTrees=100::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.8:nCuts=40:MaxDepth=4");
   // --------------------------------------------------------------------------------------------------


   factory->TrainAllMethods();
   factory->TestAllMethods();
   factory->EvaluateAllMethods();    
   
   outputFile->Close();

   delete factory;
   delete dataloader;
}

