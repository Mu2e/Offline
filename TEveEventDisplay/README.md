# Mu2e TEve Event Display Development Branch
The directory contains all the code for the TEve based mu2e Event display. This code allows 2D and 3D rendering of several of the Mu2e Data Products. It is still being developed so if we are missing some thing you need contact us (see below for details).

## The GDML

Our 3D geometry comes from the GDML stored in the src directory.

## The Module
TEveEventDisplay/src/TEveEventDisplay_module.cc is the Analyzer mdoule which currently controls the TEveEventDisplay. This is your main function. The BeginJob function sets up your Frame (a TEveMu2eMainwindow Object) and it opens an Application. This is needed for the browser to appear.
The BeginRun function calls the Frame's (TEveMu2eMainWindow Object) SetRunGeometry. This is where the GDML is accessed and the Geom_Interface used to descend nodes and desiplay the DS.
The Analyze function fills the DataCollections (a list of Mu2e Data Products are called Collections). The Filler is a Collection_Filler object where the DataCollection is filled.
The Analyze function calls the the Frame's SetEvent function. In that function the various AddProduct (e.g. AdComboHit) are called. These Add functions reside in the Data or MC Interfaces.

## Event Filter Module

The EventFilter module allows the user to call a specific starting event. This can also be done using the text boxes on the GUI.

## The fcl file
The prolog.fcl file resides in TEveEventDisplay/fcl and contains module instances for the TEveEventDisplay.

## Command Line uses
We are in the process of making this code more sophisicated. Currently you can use the ```binTEveMu2e.sh``` script for command line input (to by-pass the need to continuously change the prolog):

```TEveMu2e [fcl] [art] [number of evts] - [options]```

* ```[fcl``` the Caller fcl file which calls upon the TEveMu2e module
* ```art``` the data you would like to see represented int he display
* ```number of evts``` the nominal number of event loops that the module could reach (it works like any other). You can select a "starting" event using the EventFilter module.

```options``` are:

* Changing the Geometry: -2DOnly, -2Dand3D, -DSOnly, -CRVOnly, -DSandCRV. Note that there are default settings of 2Dand3D, DSOnly - this is displayed unless altered.
* Adding Data Products: -hits, -clusters, -tracks, -crvhits, -cosmictracks, -mctraj
* Accumulate Products i.e. collect multiple events for calibration assessments (turned off by default as only for specialist runs): -accumulate

## The Standard way to run the code
We include some example CallerFcl files in: TEveEventDisplay/CallerFcls
to run like any other ART module: ```$ mu2e -c PATH_TO_CALLER_FCL/File.fcl PATH_TO_ART/art.art --nevts 100 (for 100 events)```
The TEve Browser will appear. The first event takes a little longer as the GUI must be created.

The TEve code can be used like any other Analyzer and added to your Reco/End path as such. There is no need to use the Callers they are just guides and examples.

## The TEve Event Display Infrastructure
Current notable features of the code:

### gdml
The GDML file used here can be regenerated using: ```mu2e -c mu2eG4/fcl/gdmldump.fcl```. It contains the entire Mu2e World. We use fix.gdml as a bug in the mu2e_common.gdml was found in the early stages of this development.

### Geom Interface
Contains callers for access to Tracker and Calo geometry. This class also contains functions to set visability of different elements based on their names within the gdml.

### GeomUtils

Contains geometry transforms. The GDML was in cm so we convert all our coordinates to cm here.

### TEveMu2e basis
Contains base classes which inherit from TEve objects. This is the interface between TEve objects and mu2e products.

### Collection Filler and Data Collections
The DataCollection class has a list all the possible Mu2e data collections we might want to access. The full list is found in ```Offline/RecoDataProduct/inc```. The collections in DataCollections are set to 0 unless they are filled. The filling is done by a function ```FillRecoCollections``` in the Collection_Filler class. This is called in the module Analyze function.

### classes.h and class_def

Any src directory in mu2e which wants to use classes needs to list them in a classes.h and class_def.xml file. If you make a new class you must add it here. To make a new class use an existing class as a template. That way you wont run into errors.

### Main Window

This class sets up the GUI and imports the geometry.

## Data and MC Interfaces

These contain ```Add``` functions which add the specific data product to the event display. There is some templating but if you wish to add a new product you must do so explicitly by following the instructions in the following section (or contacting the developers).

## Adding Data Products

In order to add Data Products to the Event Display you may need to add an additional plotting function to the code. Here are the steps you should take to add the function.

1. Check if the data product you want to add is in the Collection Filler and Data Collections classes. If not then add it to these classes.

2. Create a function to add your data product in the TEveMu2eDataInterface or TEveMu2eMCInterface based on whether it is a Reco or MC Data Product. Add two lists in the header file.

3. In your function call the DataLists template function. This will take care of handling your lists so that they work with the other features in the Event Display.

4. Now add a condition to make sure your Data Product Collection is not empty.

5. If you would like to add energies call the energies template. This currently works for hit and cluster i.e. point like type data products.

6. Then loop through your Data Product Collection and call the draw function to the corresponding data type. For example, DrawHit3D/2D, DrawCluster3D/2D

7. To utilize the 3D only feature, create a condition and pass the show2D parameter and control the 2D drawing functions with that parameter.

8. Call the function in TEveMu2eMainWindow

9. Modify the Collection Filler in TEveMu2eModule to get your data product from the art file. Add the product to the fcl file and change the prolog.fcl as needed.

## Authors

This code is built upon the Mu2e Offline code however TEveEventDisplay is package developed within Mu2e Offline by Sophie Middleton (Caltech) with help from: Aditi Venkatesh (Caltech) and Namitha Chithirasee from (Pisa).

If you have any questions/comments about TEveEventDisplay please contact: smidd@caltech.edu.
