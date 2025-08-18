#ifndef DriftCalibPDF_C
#define DriftCalibPDF_C

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <ctime>
#include <limits>

#include "EventNtuple/utils/rooutil/inc/RooUtil.hh"
#include "EventNtuple/utils/rooutil/inc/common_cuts.hh"

#include "TH1.h"
#include "TFile.h"


void DriftCalibPDF(std::string filename, std::string calib_type)
{
#ifndef __CINT__

    TH1::AddDirectory(false);
    
    
    // Histograms
    double binsize = 0.025;
    double doff = binsize * 10.;
    unsigned int finefactor = 100;
    
    // double rdbinedges[2] = {0.0, 2.5};
    double rdbinedges[2] = {0.0, 3.0};
    int nrdbins = (rdbinedges[1] - rdbinedges[0]) / binsize;
    
    double rcbinedges[2] = {rdbinedges[0] - doff, rdbinedges[1] + doff};
    int nrcbins = (rcbinedges[1] - rcbinedges[0]) / binsize;
    
    TH1D* h_mc_dist = new TH1D("h_mc_dist", "", nrdbins, rdbinedges[0], rdbinedges[1]);  // Unsigned MC true DOCA between particle trajectory and wire axis
    TH1D* h_rdrift  = new TH1D("h_rdrift", "", nrdbins, rdbinedges[0], rdbinedges[1]);   // Cluster drift distance calibrated under old method
    TH1D* h_cdrift  = new TH1D("h_cdrift", "", nrcbins, rcbinedges[0], rcbinedges[1]);   // Not calibrated cluster drift distance
    
    TH1D* h_cdrift_fine = new TH1D("h_cdrift_fine", "", nrcbins*finefactor, rcbinedges[0], rcbinedges[1]);
    TH1D* h_correction  = new TH1D("h_correction", "", nrcbins, rcbinedges[0], rcbinedges[1]);   // Bin-by-bin corrections to not-calibrated cluster drift distance
    TH1D* h_cdrift_corr = new TH1D("h_cdrift_corr", "", nrcbins, rcbinedges[0], rcbinedges[1]);  // Cluster drift distance corrected
    
    
    // Set up RooUtil
    RooUtil util(filename);
    
    
    // Event and track counters
    int events_all = util.GetNEvents();
    int events_withtracks = 0;
    
    int tracks_all      = 0;
    int tracks_passcut1 = 0;
    int tracks_passcut2 = 0;
    int tracks_passcut3 = 0;
    int tracks_goodtrackhits = 0;
    
    int trackhits_all      = 0;
    int trackhits_passcut4 = 0;
    int trackhits_passcut5 = 0;
    int trackhits_passcut6 = 0;
    int trackhits_passcut7 = 0;
    
    
    // Loop over events
    // ----------------
    for ( int i_event = 0; i_event < util.GetNEvents(); ++i_event )
    {
        // Get event
        auto& event = util.GetEvent(i_event);
        
        
        // Get tracks
        auto tracks = event.GetTracks(is_e_minus);
        
        if ( event.CountTracks() <= 0 ) continue;
        events_withtracks++;
        
        
        // Loop over tracks
        // ----------------
        for ( auto& track : tracks ) {
            tracks_all++;
            
            // Get track particles
            auto particles = track.GetMCParticles(is_track_particle);
            
            
            // FIRST CUT: Adapted from 'gfit': ["dem.status>0"]
            if ( !(track.trk->status > 0) ) continue;
            tracks_passcut1++;
            
            // SECOND CUT: Also adapted from 'gfit' ["demmcsim[0].prirel._rel==0"]
            if ( !(particles.at(0).mcsim->prirel.relationship() == 0) ) continue;
            tracks_passcut2++;
            
            
            // Get track segments
            auto trksegments = track.GetSegments();
            
            // Loop over track segments
            bool bad_trksegment = false;
            for ( auto& trksegment : trksegments ) {
                if ( !has_mc_step(trksegment) ) continue;
                
                // SECOND CUT: Adapted from 'gmom': ["demmcsim[0].mom.R()-demmcvd.mom.R()<2.0"]
                if ( !(particles.at(0).mcsim->mom.R() - trksegment.trksegmc->mom.R() < 2.0) ) {
                    bad_trksegment = true;
                    break;
                }
            }
            if ( bad_trksegment ) continue;
            tracks_passcut3++;
            
            
            // Get track hits
            auto trkhits = track.GetHits();
            
            
            // Loop over track hits
            // --------------------
            bool passcut4 = false;
            bool passcut5 = false;
            bool passcut6 = false;
            bool passcut7 = false;
            
            for ( auto& trkhit : trkhits ) {
                trackhits_all++;
                
                // THIRD CUT: If track hit has no reconstructed object associated, skip it
                if ( !(trkhit.reco != nullptr) ) continue;
                bool passcut3 = true;
                trackhits_passcut4++;
                
                // FOURTH CUT: Adapted from 'ghit': ["demtsh.state>-2"]
                if ( !(trkhit.reco->state > -2) ) continue;
                bool passcut4 = true;
                trackhits_passcut5++;
                
                // FIFTH CUT: Adapted from 'thit': ["demtshmc.rel._rel==0"]
                if ( !(trkhit.mc->rel.relationship() == 0) ) continue;
                bool passcut5 = true;
                trackhits_passcut6++;
                
                // SIXTH CUT: Drift quality
                if ( trkhit.reco->driftqual <= 0.2 ) continue;
                bool passcut6 = true;
                trackhits_passcut7++;
                
                
                // Define variables
                double mc_dist = trkhit.mc->dist;
                double rdrift  = trkhit.reco->rdrift;            
                double cdrift  = trkhit.reco->cdrift;
                
                
                // Fill histograms
                h_mc_dist -> Fill(mc_dist);
                h_rdrift  -> Fill(rdrift);
                h_cdrift  -> Fill(cdrift);
                h_cdrift_fine -> Fill(cdrift);
                
            }  // End of loop over track hits
            
            tracks_goodtrackhits++;
            
        }  // End of loop over tracks
    }  // End of loop over events
    
    std::cout << std::endl;
    std::cout << " All events:         " << events_all << std::endl;
    std::cout << " Events with tracks: " << events_withtracks << std::endl;
    std::cout << std::endl;
    std::cout << " All tracks:            " << tracks_all << std::endl;
    std::cout << " Tracks passing cut #1: " << tracks_passcut1 << std::endl;
    std::cout << " Tracks passing cut #2: " << tracks_passcut2 << std::endl;
    std::cout << " Tracks passing cut #3: " << tracks_passcut3 << std::endl;
    std::cout << " Tracks with at least one good trackhit: " << tracks_goodtrackhits << std::endl;
    std::cout << std::endl;
    std::cout << " All track hits:            " << trackhits_all << std::endl;
    std::cout << " Track hits passing cut #4: " << trackhits_passcut4 << std::endl;
    std::cout << " Track hits passing cut #5: " << trackhits_passcut5 << std::endl;
    std::cout << " Track hits passing cut #6: " << trackhits_passcut6 << std::endl;
    std::cout << " Track hits passing cut #7: " << trackhits_passcut7 << std::endl;
    std::cout << std::endl;
    
    
    // Normalize histograms   
    h_mc_dist -> Scale(1.0 / h_mc_dist->Integral());
    h_rdrift  -> Scale(1.0 / h_rdrift->Integral());
    h_cdrift  -> Scale(1.0 / h_cdrift->Integral());
    h_cdrift_fine -> Scale(1.0 / h_cdrift_fine->Integral());
    
    
    // Find corrections to cluster drift distance 'cdrift'
    // ---------------------------------------------------
    
    std::vector<double> offsets(nrcbins, 0.0);  // Offset vector with length equal to number of bins of 'h_cdrift' all initialized to 0
    
    int jbin = 1;
    double ncalib = h_mc_dist->GetBinContent(jbin);  // Content of a particular bin of 'h_mc_dist'
    double cedge  = h_mc_dist->GetBinLowEdge(jbin);  // Lower edge of a particular bin of 'h_mc_dist'
    
    double calibsum = ncalib;  // Sum of bin contents of 'h_mc_dist'
    double nraw = 0.0;
    double rawsum = nraw;      // Sum of bin contents of 'h_cdrift'
    
    // Loop over bins of 'h_cdrift'
    for ( int ibin = 1; ibin <= nrcbins; ++ibin ) {
        double bincontent = h_cdrift->GetBinContent(ibin);  // Bin content
        double bincenter  = h_cdrift->GetBinCenter(ibin);   // Bin center
        
        rawsum += bincontent;  // Add bin content to 'rawsum'
        offsets[ibin-1] = bincenter - (cedge + binsize*((nraw + 0.5*bincontent)/ncalib));  // Fill offset vector entry corresponding to this bin
        
        // If 'nraw' plus 'h_cdrift' bin content is less than 'h_mc_dist' bin content...
        if ( nraw + bincontent < ncalib ) {
            nraw += bincontent;  // Add 'h_cdrift' bin content to 'nraw'
        }
        else {
            nraw += bincontent - ncalib;  // Otherwise, add to 'nraw' the difference between bin contents of 'h_cdrift' minus 'h_mc_dist'
            ++jbin;   // Go to the next bin of 'h_mc_dist'
            
            // If the bin number of 'h_mc_dist' is too big, or the corresponding bin content is less than zero, break the loop
            if ( jbin > nrdbins || h_mc_dist->GetBinContent(jbin) <= 0 ) break;
            
            // Otherwise, prepare 'ncalib' and 'cedge' with the info from the next bin of 'h_mc_dist'
            ncalib = h_mc_dist->GetBinContent(jbin);
            cedge  = h_mc_dist->GetBinLowEdge(jbin);
            
            calibsum += ncalib;  // Add to 'calibsum' the content of such bin of 'h_mc_dist'
        }
    }
    
    
    // Perform a safety check
    if ( jbin < nrcbins && (std::fabs(1.0 - calibsum) > 1.0e-6 || std::fabs(1.0 - rawsum) > 1.0e-6) ) {
        std::cout << " Early exit! " << std::endl;
        std::cout << " \tjbin                  = " << jbin << std::endl;
        std::cout << " \tRemaining MC fraction = " << 1.0 - calibsum << std::endl;
        std::cout << " \tcdrift fraction       = " << 1.0 - calibsum << std::endl;
    }
    
    
    // Fill bin-by-bin corrections and corrected drift histograms
    for ( unsigned int ibin = 0; ibin < offsets.size(); ++ibin ) {
        h_correction -> SetBinContent(ibin+1, offsets[ibin]);
    }
    
    for ( unsigned int ibin = 1; ibin <= finefactor*nrcbins; ++ibin ) {
        unsigned int jbin = floor(float(ibin)/float(finefactor));
        h_cdrift_corr -> Fill(h_cdrift_fine->GetBinCenter(ibin+1) - offsets[jbin], h_cdrift_fine->GetBinContent(ibin+1));
    }
    h_cdrift_corr -> Scale(1.0 / h_cdrift_corr->Integral());
    
    
    // // Save calibration file
    // std::string cfname = "DriftCalibPDF";
    // ofstream cfile(Form("%s_%s.txt", cfname.c_str(), calib_type.c_str()), ios::trunc);
    
    // time_t now = time(0);
    // char* dt = ctime(&now);
    
    // cfile << "# The following was produced by DriftCalibPDF.C with track hit selection 'trkhit.reco->driftqual > 0.2' on " << dt << std::endl;
    // cfile << std::setw(4) << std::setprecision(3);
    // cfile << "driftOffBins : [ ";
    // cfile << rcbinedges[0] << " , " << rcbinedges[1] << " ]" << endl;
    // cfile << "driftOffset : [ ";
    
    // bool first(true);
    // for( int ibin = 0; ibin < nrcbins; ++ibin ) {
    //     if ( !first ) cfile << " , ";
    //     first = false;
    //     cfile << offsets[ibin];
    // }
    // cfile << " ]" << endl;
    // cfile.close();
    
    
    // Save histograms
    TFile outputfile(Form("../root_files/DriftCalibPDF_%s.root", calib_type.c_str()), "RECREATE");
    
    h_mc_dist -> Write();
    h_rdrift  -> Write();
    h_cdrift  -> Write();
    
    h_cdrift_fine -> Write();
    h_correction  -> Write();
    h_cdrift_corr -> Write();
    
    outputfile.Close();
    
#endif  // __CINT__
}

#endif  // DriftCalibPDF_C