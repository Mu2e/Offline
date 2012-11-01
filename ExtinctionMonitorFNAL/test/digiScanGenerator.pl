#!/usr/bin/perl -w
#
# Andrei Gaponenko, 2012

use strict;
use File::Basename;
use Getopt::Long;
use Fcntl qw(O_CREAT O_EXCL O_WRONLY);

#================================================================
sub usage() {
    my $self = basename($0);
    return <<EOF
Usage: 
	$self \\
	       --scanpar=name \\
	       --npoints=int \\
	       --start=float \\
	   {--step=float|--factor=float}


Where the {alternative1|alternative2} constructs denote mutually
exclusive alternatives.  The --step option requests a linear scan,
and the --factor option an exponential one.  The script will create
npoints fcl files in the current directory.

EOF
;
}

#================================================================
sub makeFileName($) {
    my $pars = shift;
    my $res = "scan";
    foreach my $key (sort(keys %$pars)) {
	#print $key, '=', $$pars{$key}, "\n";
	$res .= "_$key$$pars{$key}";
    }
    return $res;
}

#================================================================
sub header() {
    return <<EOF
#include "ExtinctionMonitorFNAL/test/digiTuning.fcl"
EOF
}

#================================================================
sub writeConfig($) {
    my $pars = shift;
    my $fn = makeFileName($pars) . ".fcl";
    print "*** Writing config $fn\n";

    my $threshold = $$pars{'threshold'};
    my $nclusters = $$pars{'nclusters'};

    my $current = $$pars{'current'};
    # Current is a combination of totCalib and qCalib.
    # We can fix one of them and vary another.
    my $totCalib = 1;
    my $qCalib = $threshold + $totCalib * $current;

    sysopen(my $out, $fn, O_CREAT|O_EXCL|O_WRONLY ) or die "Can't open $fn for writing: $!\n";
    print $out header();
    print $out "physics.producers.pixelDigitization.discriminatorThreshold : $threshold\n";
    print $out "physics.producers.pixelDigitization.qCalib : $qCalib\n";
    print $out "physics.producers.pixelDigitization.totCalib : $totCalib\n";
    print $out "physics.producers.pixelDigitization.numClustersPerHit : $nclusters\n";
}

#================================================================
# parameters we care about, with their default values
my %pars = (
    threshold => 4500.,
    current => 16000.,
    nclusters => 128,
    );

my %opt = (help => 0,);

# Process command line opts.
GetOptions(\%opt, 
	   "help",
	   "scanpar=s",
	   "npoints=i",
	   "start=f",
	   "step=f",
	   "factor=f",
	   ) 
    or die "\nError processing command line options.\n";

if($opt{'help'}) {
    print usage();
    exit 0;
}

# Check that all of the required args are present:
foreach my $k ('scanpar', 'npoints', 'start') {
    defined $opt{$k} or die "Error: --$k must be specified.  Try the --help option.\n";
}

my $scanpar = $opt{'scanpar'};
my $npoints = $opt{'npoints'};
my $start = $opt{'start'};

#----------------
my $logscan = 0;
if(defined $opt{'step'}) {
    my $step = $opt{'step'};

    die "Error: --step and --factor are mutually exclusive.\n"
	if defined $opt{'factor'};

    my $count = 0;
    my $val = $start;
    while($count++ < $npoints) {
	$pars{$scanpar} = $val;
	writeConfig(\%pars);
	$val += $step;
    }

}
else {
    die "Error: either --step or --factor must be specified.\n"
	if !defined $opt{'factor'};

    my $factor = $opt{'factor'};

    my $count = 0;
    my $val = $start;
    while($count++ < $npoints) {
	$pars{$scanpar} = $val;
	writeConfig(\%pars);
	$val *= $factor;
    }
}

exit 0;

