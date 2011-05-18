# $Id: MessageLogger_cfi.py,v 1.2 2011/05/18 02:27:15 wb Exp $
# $Author: wb $
# $Date: 2011/05/18 02:27:15 $
#
# Original author Rob Kutschke
#

import FWCore.ParameterSet.python.Config as mu2e

MessageLogger = mu2e.Service("MessageLogger",
    suppressInfo = mu2e.untracked.vstring(),
    debugs = mu2e.untracked.PSet(
        placeholder = mu2e.untracked.bool(True)
    ),
    suppressDebug = mu2e.untracked.vstring(),
    cout = mu2e.untracked.PSet(
        placeholder = mu2e.untracked.bool(True)
    ),
    warnings = mu2e.untracked.PSet(
        placeholder = mu2e.untracked.bool(True),
        threshold = mu2e.untracked.string('WARNING'),
        default = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(-1),
        ),
    ),
    default = mu2e.untracked.PSet(
        limit = mu2e.untracked.int32(-1)
    ),
    errors = mu2e.untracked.PSet(
        placeholder = mu2e.untracked.bool(True),
        threshold = mu2e.untracked.string('ERROR'),
        default = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(-1),
        ),
    ),
    cerr = mu2e.untracked.PSet(

        # Default imit for info level messages
        INFO = mu2e.untracked.PSet(
           limit = mu2e.untracked.int32(5)
        ),

        # Default imit for debug level messages
        DEBUG = mu2e.untracked.PSet(
           limit = mu2e.untracked.int32(5)
        ),

        noTimeStamps = mu2e.untracked.bool(False),

        # Per event heartbeat from the InputSource.
        FwkReport = mu2e.untracked.PSet(
            reportEvery = mu2e.untracked.int32(1),
            limit = mu2e.untracked.int32(5)
        ),

        # Defines default for all severities and categories, unless overridden.
        default = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(6)
        ),
        Root_NoDictionary = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(0)
        ),
        FwkJob = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(0)
        ),
        FwkSummary = mu2e.untracked.PSet(
            reportEvery = mu2e.untracked.int32(1),
            limit = mu2e.untracked.int32(10000000)
        ),

        # Open/close file messages.
        fileAction = mu2e.untracked.PSet(
            reportEvery = mu2e.untracked.int32(1),
            limit = mu2e.untracked.int32(-1)
        ),
        threshold = mu2e.untracked.string('INFO')
    ),
    FrameworkJobReport = mu2e.untracked.PSet(
        default = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(0)
        ),
        FwkJob = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(10000000)
        )
    ),
    suppressWarning = mu2e.untracked.vstring(),
    statistics = mu2e.untracked.vstring('cerr_stats'),
    cerr_stats = mu2e.untracked.PSet(
        threshold = mu2e.untracked.string('WARNING'),
        output = mu2e.untracked.string('cerr')
    ),
    infos = mu2e.untracked.PSet(
        Root_NoDictionary = mu2e.untracked.PSet(
            limit = mu2e.untracked.int32(0)
        ),
        placeholder = mu2e.untracked.bool(True)
    ),
    destinations = mu2e.untracked.vstring('warnings',
        'errors',
        'infos',
        'debugs',
        'cout',
        'cerr'),
    debugModules = mu2e.untracked.vstring(),
    categories = mu2e.untracked.vstring('FwkJob',
        'FwkReport',
        'FwkSummary',
        'fileAction',
        'Root_NoDictionary'),
    fwkJobReports = mu2e.untracked.vstring('FrameworkJobReport')
)


