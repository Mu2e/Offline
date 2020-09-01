# Mu2e TEve Event Display Development Branch
This branch contains the prototype for the TEve based Mu2e Event Display.
Please read the wiki tab for details.
Here are just a few basic details.

## The Module
TEveEventDisplay/src/TEveEventDisplay_module.cc is the Analyzer mdoule which currently controls the TEveEventDisplay. This is your main function. The BeginJob function sets up your Frame (a TEveMu2eMainwindow Object) and it opens an Application. This is needed for the browser to appear.
The BeginRun function calls the Frame's (TEveMu2eMainWindow Object) SetRunGeometry. This is where the GDML is accessed and the Geom_Interface used to descend nodes and desiplay the DS.
The Analyze function fills the DataCollections ( a list of Mu2e Data Products are called Collections). The Filler is a Collection_Filler object where the DataCollection is filled.
The Analyze function calls the the Frame's SetEvent function. In that function the various AddProduct (e.g. AdComboHit) are called.

## The fcl file
The prolog.fcl file resides in TEveEventDisplay/fcl and contains whats known as module instances for the TEveEventDisplay. Currently these are for : helix tracks, cosmic tracks or calo only events. We can add more.

## Running the code
You should access some of the Art examples in TEveEventDisplay/ArtExamples. You should access fcl files in TEveEventDisplay/CallerFcls
to run: ```$ mu2e -c PATH_TO_CALLER_FCL/File.fcl PATH_TO_ART/art.art --nevts 100 (for 100 events)```
The TEve Browser will appear. The first event takes a little longer as the GUI must be created.

## The TEve Event Display Infrastructure
Currently there are a number of Interfaces defined in ``dict classes`` in the ``src`` directory:

### gdml
The GDML file used here can be regenerated using: ```mu2e -c mu2eG4/fcl/gdmldump.fcl```. It contains the entire Mu2e World. We use fix.gdml as a bug in the mu2e.gdml was found in the early stages of this development.

### Geom Interface
Contains callers for access to Tracker and Calo geometry. This class also contains functions to set visability of different elements based on their names within the gdml.

### TEveMu2e basis
Contains base classes which inherit from TEve objects. This is the interface between TEve objects and mu2e products. 

### Collection Filler and Data Collections
The DataCollection class has a list all the possible Mu2e data collections we might want to access. The full list is found in ```Offline/RecoDataProduct/inc```. The collections in DataCollections are set to 0 unless they are filled. The filling is done by a function ```FillRecoCollections``` in the Collection_Filler class. This is called in the module Analyze function.

## classes.h and class_def

Any src directory in mu2e which wants to use classes needs to list them in a classes.h and class_def.xml file. If you make a new class you must add it here. To make a new class use an existing class as a template. That way you wont run into errors.

### Main Window

This class sets up the Gui and imports the geometry. Here the plotting data functions are currently called. you should base any work on AddComboHits.

### Adding Data Products

In order to add Data Products to the Event Display you may need to add an additional plotting function to the code. Here are the steps you should take to add the function.

1. Check if the data product you want to add is in the Collection Filler and Data Collections classes. If not then add it to these classes.

2. Create a function to add your data product in the TEveMu2eDataInterface or TEveMu2eMCInterface based on whether it is a Reco or MC Data Product. Add two lists in the header file.

3. In your function call the DataLists template function. This will take care of handling your lists so that they work with the other features in the Event Display.

4. Now add a condition to make sure your Data Product Collection is not empty.

5. If you would like to add energies call the energies template. This currently works for hit and cluster type data products.

6. Then loop through your Data Product Collection and call the draw function to the corresponding data type. For example, DrawHit3D/2D, DrawCluster3D/2D

7. To utilize the 3D only feature, create a condition and pass the show2D parameter and control the 2D drawing functions with that parameter.

8. Call the function in TEveMu2eMainWindow 

9. Modify the Collection Filler in TEveMu2eModule to get your data product from the art file. Add the product to the fcl file and change the prolog.fcl as needed.

## Authors

This code is built upon the Mu2e Offline code however TEveEventDisplay is package developed within Mu2e Offline by Sophie Middleton with help from Aditi Venkatesh both from Caltech. If you have any questions/comments about TEveEventDisplay place contact: smidd@caltech.edu.
