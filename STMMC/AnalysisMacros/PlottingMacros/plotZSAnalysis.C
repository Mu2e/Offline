// Plots the stages of the ZS analysis
// See plotZSGradients.sh for usage examples
// Original author - Pawel Plesniak

#include <TError.h>
#include <iostream>
#include <limits>
#include <math.h>

void customErrorHandler(int level, Bool_t abort, const char* location, const char* message) {
    /*
        Description
            Define a custom error handler that won't print the stack trace but will print an error message and exit.
    */
    std::cerr << message << std::endl;
    if (level > kInfo)
        exit(1);
};

void collectAnalysisData(const std::string fileName, const std::string treeName, std::vector<int16_t> &ADCs, std::vector<int16_t> &gradients, std::vector<double> &averagedGradients, std::vector<unsigned int> &eventIds) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotZSAnalysis"
            treeName - as documented in function "plotZSAnalysis"
            ADCs - as documented in function "plotZSAnalysis"
            gradients - as documented in function "plotZSAnalysis"
            averagedGradients - as documented in function "plotZSAnalysis"
            eventIds - as documented in function "plotZSAnalysis"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataADC - ADC from file
            ADCMin - minimum ADC value allowed
            dataGradient - gradient from file
            dataAveragedGradient - averaged gradient from file
            dataEventId - event ID from file
            entries - number of entries in the TTree
            i - iterator variable
    */
    std::cout << "Processing file " << fileName << std::endl;

    // Get the branch
    TFile *file = new TFile(fileName.c_str());
    if (!file || file->IsZombie()) {
        Fatal("collectData", "Failed to open the file.");
    };
    TTree *tree = (TTree*)file->Get(treeName.c_str());
    if (!tree)
        Fatal("collectData", "Requested tree does not exist in the file.");

    // Get the list of branches to check if they exist
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<std::string> branchNames;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch)
            branchNames.push_back(branch->GetName());
    };

    // If the branches exist, assign them to the appropriate variables
    const int16_t ADCMin = (-1 * std::pow(2, 15) + 1), overlapThreshold = 100;
    int16_t dataADC, dataGradient;
    double dataAveragedGradient;
    unsigned int dataEventId;
    if (std::find(branchNames.begin(), branchNames.end(), "ADC") != branchNames.end())
        tree->SetBranchAddress("ADC", &dataADC);
    else
        Fatal("collectData", "Requested branch 'ADC' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "gradient") != branchNames.end())
        tree->SetBranchAddress("gradient", &dataGradient);
    else
        Fatal("collectData", "Requested branch 'gradient' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "averagedGradient") != branchNames.end())
        tree->SetBranchAddress("averagedGradient", &dataAveragedGradient);
    else
        Fatal("collectData", "Requested branch 'averagedGradient' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "eventId") != branchNames.end())
        tree->SetBranchAddress("eventId", &dataEventId);
    else
        Fatal("collectData", "Requested branch 'eventId' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        ADCs.push_back(dataADC);
        gradients.push_back(dataGradient > overlapThreshold ? ADCMin : dataGradient);
        averagedGradients.push_back(dataAveragedGradient > overlapThreshold ? ADCMin : dataAveragedGradient);
        eventIds.push_back(dataEventId);
    };

    // Check if data has been collected
    if (eventIds.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void collectResultData(std::string fileName, std::string treeName, std::vector<int16_t> &ADCs, std::vector<double> &times, std::vector<unsigned int> &eventIds) {
    /*
        Description
            Collects all the required data from virtual detector TTrees

        Arguments
            fileName - as documented in function "plotZSAnalysis"
            treeName - as documented in function "plotZSAnalysis"
            ADCs - as documented in function "plotZSAnalysis"
            times - as documented in function "plotZSAnalysis"
            eventIds - as documented in function "plotZSAnalysis"

        Variables
            file - ROOT TFile interface
            branches - ROOT TFile interface to TTree branches
            branchNames - vector of branch names
            dataADC - ADC from file
            dataTime - time from file
            dataEventId - event ID from file
            entries - number of entries in the TTree
            tADC - ADC clock tick length [ns]
            i - iterator variable
    */
    std::cout << "Processing file " << fileName << std::endl;

    // Get the branch
    TFile *file = new TFile(fileName.c_str());
    if (!file || file->IsZombie()) {
        Fatal("collectData", "Failed to open the file.");
    };
    TTree *tree = (TTree*)file->Get(treeName.c_str());
    if (!tree)
        Fatal("collectData", "Requested tree does not exist in the file.");

    // Get the list of branches to check if they exist
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<std::string> branchNames;
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *branch = dynamic_cast<TBranch*>(branches->At(i));
        if (branch)
            branchNames.push_back(branch->GetName());
    };

    // If the branches exist, assign them to the appropriate variables
    int16_t dataADC;
    uint32_t dataTime;
    unsigned int dataEventId;
    if (std::find(branchNames.begin(), branchNames.end(), "ADC") != branchNames.end())
        tree->SetBranchAddress("ADC", &dataADC);
    else
        Fatal("collectData", "Requested branch 'ADC' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "time") != branchNames.end())
        tree->SetBranchAddress("time", &dataTime);
    else
        Fatal("collectData", "Requested branch 'time' does not exist in the file.");
    if (std::find(branchNames.begin(), branchNames.end(), "eventId") != branchNames.end())
        tree->SetBranchAddress("eventId", &dataEventId);
    else
        Fatal("collectData", "Requested branch 'eventId' does not exist in the file.");

    // Get the number of entries
    int entries = tree->GetEntries();

    // Choose the data entries corresponding to the chosen time
    const double tADC = 3.125;

    // Collect the data
    for (int i = 0; i < entries; i++) {
        // Get the corresponding entry
        tree->GetEntry(i);

        // Collect the data
        ADCs.push_back(dataADC);
        times.push_back(dataTime * tADC);
        eventIds.push_back(dataEventId);
    };

    // Check if data has been collected
    if (eventIds.size() == 0)
        Fatal("collectData", "No data was collected from this file");

    // Clear
    file->Close();
    delete file;
    std::cout << "Finished processing file " << fileName << std::endl;
    return;
};

void plot(std::vector<int16_t> &inputADCs, std::vector<int16_t> &gradients, std::vector<double> &averagedGradients, std::vector<int16_t> &outputADCs, std::vector<double> outputTimes, unsigned int eventId, const double tBefore, const double tAfter, const double threshold, const int window, const int nAverage, const double tMin, const double tMax) {
    /*
        Description
            Plots the ZS waveforms, their analysis steps, and the results

        Arguments
            inputADCs - as documented in function "plotZSAnalysis"
            gradients - as documented in function "plotZSAnalysis"
            averagedGradients - as documented in function "plotZSAnalysis"
            outputADCs - as documented in function "plotZSAnalysis"
            outputTimes - as documented in function "plotZSAnalysis"
            eventId - as documented in function "plotZSAnalysis"
            tBefore - as documented in function "plotZSAnalysis"
            tAfter - as documented in function "plotZSAnalysis"
            threshold - as documented in function "plotZSAnalysis"
            window - as documented in function "plotZSAnalysis"
            nAverage - as documented in function "plotZSAnalysis"
            tMin - as documented in function "plotZSAnalysis"
            tMax - as documented in function "plotZSAnalysis"

        Variables
            nInput - number of input ADC values
            nOutput - number of output ADC values
            time - data product used as x axis for plotting the full ADC range
            tADC - time interval of the ADC clock tick [ns]
            tMaxFull - time at the end of the input ADC. Note time is modelled from 0 for the input ADCs
            tMinCutIndex - index of the minimum time
            tMaxCutIndex - index of the maximum time
            tMaxCut - maximum time plotted
            time - vector of ADC times for input data
            inputADCsA - array of input ADCs
            outputADCsA - array of output ADCs
            gradientsA - array of analysis gradients
            averagedGradientsA - array of analysis averaged gradients
            outputTimesA - array of result output itmes
            nInput - number of accepted input data entries
            nOutput - numnber of accepted output data entries
            c - canvas used to generate the plot
            gInputADCs - graph of input ADCs
            gGradients - graph of analysis gradients
            gAveragedGradient - graph of analysis averaged gradients
            splitIndexes - indexes of discontinuities in the output ADC times
            splitOutputTime - vector of times used to generate the output waveforms
            splitOutputADCs - vector of ADCs used to generate the output waveofrms
            nOutputPlots - number of output waveform plots to generate. Need to be separate as otherwise the TGraph connects the discontinuities displaying data that is not there
            outputPlots - plots of output waveforms
            x1, x2, y1, y1 - coordinates of legend as a percentage of the canvas area
            legend - legend for the plot
            cFileName - filename the plot is saved as
    */

    std::cout << "Max averaged gradient " << *std::max_element(averagedGradients.begin(), averagedGradients.end()) << std::endl; // deleteme
    // Count the number of entries used to generate this plot
    int nInput = inputADCs.size(), nOutput = outputADCs.size();

    // Construct the time data used for the ZS algoritm analysis stages
    const double tADC = 3.125, tMaxFull = tADC * nInput;

    // Sanity checks
    if (tMin > tMaxFull)
        Fatal("plot", "Selected minimum cut is greater than the maximum, exiting.");

    // Convert the times to the number of clock ticks for selecting the range of entries to plot
    const int tMinCutIndex = static_cast<int>(tMin / tADC);
    const int tMaxCutIndex = (tMaxFull < tMax || tMax < std::numeric_limits<double>::epsilon()) ? nInput : static_cast<int>(tMax / tADC); // if the time is zero or if the time is greater than the max, set the index to the number of entries in the dataset
    const double tMaxCut = tMaxCutIndex * tADC;
    // Sanity checks
    if (tMaxCutIndex - tMinCutIndex < 10)
        Fatal("plot", "Selected data range is such that there are very few entries to plot, exiting.");

    // Construct the time vector
    std::vector<double> time;
    for (int i = tMinCutIndex; i < tMaxCutIndex; i++)
        time.push_back(i * tADC);

    // Construct the variables used to store the data relevant for plotting
    std::vector<double> inputADCsA, outputADCsA, gradientsA, averagedGradientsA, outputTimesA;

    // Select and convert the relevant data for plotting
    for (int i = tMinCutIndex; i < tMaxCutIndex; i++) {
        inputADCsA.push_back(static_cast<double>(inputADCs[i]));
        gradientsA.push_back(static_cast<double>(gradients[i]));
        averagedGradientsA.push_back(averagedGradients[i]);
    };
    nInput = inputADCsA.size();
    double t = 0.0;
    for (int i = 0; i < nOutput; i++) {
        t = outputTimes[i];
        if ((tMin < std::numeric_limits<double>::epsilon() || tMin < t) && (tMax < std::numeric_limits<double>::epsilon() || t < tMax)) {
            outputADCsA.push_back(static_cast<double>(outputADCs[i]));
            outputTimesA.push_back(t);
        };
    };
    nOutput = outputADCsA.size();
    if (nOutput == 0)
        Fatal("plot", "Number of output points is zero, exiting");

    // Set up the canvas for the plot
    TCanvas* c = new TCanvas("c", "c", 1500, 1000);
    gPad->SetLeftMargin(0.14);
    gPad->SetRightMargin(0.08);

    // Plot the input ADCs
    TGraph *gInputADCs = new TGraph(nInput, time.data(), inputADCsA.data());
    gInputADCs->GetXaxis()->SetLimits(tMin, tMaxCut);
    gInputADCs->Draw("APL");
    gInputADCs->SetTitle("Zero suppression;Time [ns];ADC [arb. unit]");

    // Plot the gradients
    TGraph *gGradients = new TGraph(nInput, time.data(), gradientsA.data());
    gGradients->Draw("PL SAME");
    gGradients->SetLineColor(kRed);

    // Plot the analysis gradients
    TGraph *gAveragedGradients = new TGraph(nInput, time.data(), averagedGradientsA.data());
    gAveragedGradients->Draw("PL SAME");
    gAveragedGradients->SetLineColor(kBlue);

    // Plot the zero-suppressed waveforms
    // Determine the start indexes of separate waveforms
    std::vector<uint> splitIndexes;
    for (uint i = 1; i < nOutput; i++) {
        if ((outputTimesA[i] - outputTimesA[i - 1]) > tADC)
            splitIndexes.push_back(i);
    };

    // Construct variables to hold the plot data
    std::vector<double> splitOutputTimes, splitOutputADCs;
    const int nOutputPlots = splitIndexes.size();

    // Construct the plots of zero-suppressed waveforms
    std::vector<TGraph*> outputPlots;
    for (uint i = 1; i < nOutputPlots; i++) {
        splitOutputTimes.assign(outputTimesA.begin() + splitIndexes[i - 1], outputTimesA.begin() + splitIndexes[i]);
        splitOutputADCs.assign(outputADCsA.begin() + splitIndexes[i - 1], outputADCsA.begin() + splitIndexes[i]);
        outputPlots.emplace_back(new TGraph(splitIndexes[i] - splitIndexes[i - 1], splitOutputTimes.data(), splitOutputADCs.data()));
        outputPlots.back()->Draw("PL SAME");
        outputPlots.back()->SetLineColor(kGreen + 2);
        outputPlots.back()->SetLineWidth(2);
        splitOutputTimes.clear();
        splitOutputADCs.clear();
    };

    // Construct the legend, assign all the generated plots
    const double x1 = 0.7, x2 = 0.9, y1 = 0.4, y2 = 0.6;
    TLegend *legend = new TLegend(x1, y1, x2, y2);
    legend->SetBorderSize(1);
    legend->SetFillColor(0);
    legend->SetTextSize(0.04);
    legend->AddEntry(gInputADCs, "Input ADCs", "l");
    legend->AddEntry(gGradients, "Grad.", "l")->SetTextColor(kRed);
    legend->AddEntry(gAveragedGradients, "Avg. Grad.", "l")->SetTextColor(kBlue);
    legend->AddEntry(outputPlots.back(), "Output ADCs", "l")->SetTextColor(kGreen + 2);
    legend->Draw();
    c->Update();

    // Save the generated plot
    std::string cFileName = "ZS.Event" + std::to_string(eventId) + ".png";
    c->SaveAs(cFileName.c_str());

    // Cleanup
    gInputADCs->Delete();
    gGradients->Delete();
    gAveragedGradients->Delete();
    outputPlots.clear();
    legend->Delete();
    c->Close();

    return;
};

std::vector<unsigned int> makeUniqueEventIds(std::vector<unsigned int> &inputEventIds, std::vector<unsigned int> &outputEventIds) {
    /*
        Description
            Generates a unique intersection vector of the event IDs retrieved from the input files

        Argument
            inputEventIds - as documented in function "plotZSAnalysis"
            outputEventIds - as documented in function "plotZSAnalysis"

        Variables
            inputEventIdSet - set of event IDs from ZS input data
            uniqueInputEventIds - vector of unique event IDs from ZS input data
            outputEventIdSet - set of event IDs from ZS output data
            uniqueOutputEventIds - vector of unique event IDs from ZS output data
            overlap - vector of unique event IDs present in both input and output data
    */
    std::cout << "\nMaking the event IDs unique\n" << std::endl;

    // Get a unqiue set of input event IDs
    std::set<unsigned int> inputEventIdSet(inputEventIds.begin(), inputEventIds.end());

    // Convert the set to a vector
    std::vector<unsigned int> uniqueInputEventIds(inputEventIdSet.begin(), inputEventIdSet.end());

    // Get a unique set of output event IDs
    std::set<unsigned int> outputEventIdSet(outputEventIds.begin(), outputEventIds.end());

    // Convert the set to a vector
    std::vector<unsigned int> uniqueOutputEventIds(outputEventIdSet.begin(), outputEventIdSet.end());

    // Find the intersection of the two unique event ID vectors
    std::vector<unsigned int> overlap;
    std::set_intersection(uniqueInputEventIds.begin(), uniqueInputEventIds.end(), uniqueOutputEventIds.begin(), uniqueOutputEventIds.end(), std::back_inserter(overlap));

    // Return the unqiue intersection
    return overlap;
}

void plotZSAnalysis(const std::string analysisFileName, const std::string resultFileName, const std::string treeName, const unsigned int eventID = 0, const double tBefore = 2000, const double tAfter = 10000, const double threshold = -100, const int window = 100, const int nAverage = 5, const double tMin = 0.0, const double tMax = 0.0) {
    /*
        Description
            Generates a plot of the ZS waveform analysis with file name
                ZS.Event<EventID>.png
            <EventID> is allocated even if the parameter "eventID" is not used.

        Parameters
            fileName - name of the root file generated with STMZeroSuppression_module.cc
            treeName - name of the ttree containing all the data in fileName
            eventID - controls what eventID to plot for. If 0, generates plots for all events
            tBefore - amount of data kept before the averaged gradient crosses the threshold [ns]
            tAfter - amount of datta kept after the averaged graedient crosses the threshold [ns]
            threshold - averaged gradient threshold [ADC]
            window - window over which the average was applied
            nAverage - number of entries used for the average
            tMin - minimum ttime to plot [ns]. If zero, does not apply a time cut [ns]
            tMax - maximum time to plot [ns]. If zero, does not apply a time cut [ns]

        Variables
            inputADCs - vector of ADC values used as input to ZS algorithm
            outputADCs - vector of ADC values stored by ZS algorithm
            gradients - vector of gradients calculated with ZS algorithm
            averagedGradients - vector of averaged gradients calculated with ZS algorithm
            outputTimes - vector of times associated with outputADCs
            inputEventIds - vector of eventIds used as input to ZS algorithm
            outputEventIds - vector of eventIds saved by ZS algorithm
            nInput - number of ADC entries imported from the ZS algorithm input
            nOutput - number of ADC entried imported from the ZS algorithm output
            plotInputADCs - selected input ADCs to use when generating the input plot
            plotOutputADCs - selected output ADCs to use when generating the input plot
            plotGradients - selected gradients to use when generating the input plot
            plotAveragedGradients - selected averaged gradients to use when generating the input plot
            plotOutputTimes - selected output times to use when generating the input plot
            plotEventIds - selected event IDss to use when generating the input plot
            plotEventId - iterator variable for plotEventIds

    */

    // Update global parameters
    SetErrorHandler(customErrorHandler);
    gROOT->SetBatch(kTRUE);

    // Sanity checks
    if (tMin < 0)
        Fatal("plotZSAnalysis", "tMin must be positive!");
    if (tMax < 0)
        Fatal("plotZSAnalysis", "tMax must be positive!");
    if (tMax < tMin)
        Fatal("plotZSAnalysis", "tMax must be greater than tMin!");

    // Construct the variables used to hold the input data
    std::vector<int16_t> inputADCs, outputADCs, gradients;
    std::vector<double> averagedGradients, outputTimes;
    std::vector<unsigned int> inputEventIds, outputEventIds;

    // Read the data from the input files
    collectAnalysisData(analysisFileName, treeName, inputADCs, gradients, averagedGradients, inputEventIds);
    collectResultData(resultFileName, treeName, outputADCs, outputTimes, outputEventIds);
    int nInput = inputADCs.size(), nOutput = outputADCs.size();

    // Construct the variables used to hold the data used for plotting
    std::vector<int16_t> plotInputADCs, plotOutputADCs, plotGradients;
    std::vector<double> plotAveragedGradients, plotOutputTimes;
    std::vector<unsigned int> plotEventIds = makeUniqueEventIds(inputEventIds, outputEventIds);

    // Validate the event ID in the case that it is not present in the input files
    if (eventID != 0) {
        if (std::find(plotEventIds.begin(), plotEventIds.end(), eventID) == plotEventIds.end())
            Fatal("plotZSAnalysis", "Requested eventID is not one of those from the file, exiting.");
        plotEventIds.clear();
        plotEventIds.push_back(eventID);
    };

    // Generate the plots
    for (unsigned int plotEventId : plotEventIds) {
        std::cout << "\nGenerating plot for event with ID " << plotEventId << std::endl;

        // Select the relevant data for plotting
        for (int i = 0; i < nInput; i++) {
            if (inputEventIds[i] == plotEventId) {
                plotInputADCs.push_back(inputADCs[i]);
                plotGradients.push_back(gradients[i]);
                plotAveragedGradients.push_back(averagedGradients[i]);
            };
        };
        for (int i = 0; i < nOutput; i++) {
            if (outputEventIds[i] == plotEventId) {
                plotOutputADCs.push_back(outputADCs[i]);
                plotOutputTimes.push_back(outputTimes[i]);
            };
        }

        // Sanity checks
        if (plotInputADCs.empty())
            Fatal("plotZSAnalysis", "Empty input ADC vector, exiting!");
        if (plotGradients.empty())
            Fatal("plotZSAnalysis", "Empty gradient vector, exiting!");
        if (plotAveragedGradients.empty())
            Fatal("plotZSAnalysis", "Empty averaged gradient vector, exiting!");
        if (plotOutputADCs.empty())
            Fatal("plotZSAnalysis", "Empty output ADC vector, exiting!");
        if (plotOutputTimes.empty())
            Fatal("plotZSAnalysis", "Empty output time vector, exiting!");

        // Generate the plots
        plot(plotInputADCs, plotGradients, plotAveragedGradients, plotOutputADCs, plotOutputTimes, plotEventId, tBefore, tAfter, threshold, window, nAverage, tMin, tMax);

        // Prepare the buffer variables for the next event ID
        plotInputADCs.clear();
        plotGradients.clear();
        plotAveragedGradients.clear();
        plotOutputADCs.clear();
        plotOutputTimes.clear();
    };
    return;
};
