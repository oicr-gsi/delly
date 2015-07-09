#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use constant DEBUG=>0;

=head2 Preparing inputs for HMMcopy

 HMMcopy is a tool for analyzing copy-number variation
 This script prepares inputs for the main part (analysis)
 It is very simple, the main purpose of it - 
 make sure that we don't produce duplicate files

 ./convertHMMcopy.pl --input input.bam --output output.wig --read-counter /bin/HMMcopy/bin/readCounter

 will write into output.wig

=cut

my $USAGE = "convertHMMcopy.pl --input [input bam] --output [output wig] --read-counter [path to HMMcopy readCounter binary]";

# Required parameters
my($input,$output,$readcounter);
my $results = GetOptions ("read-counter=s"  =>  \$readcounter,
                          "input=s"    =>  \$input,
                          "output=s"   =>  \$output);

if (!$input || !$output || !$readcounter) { die $USAGE; }
if (-e $output && -s $output) {
    print STDERR "Output [$output] exists, won't overwrite\n";
    exit;
}

# Make output directory if it does not exists                            
my $outdir = dirname($output);
`mkdir -p $outdir` if (!-e $outdir || !-d $outdir);


#=====================================
# RUNNING ReadCounter here
#=====================================
print STDERR "Running $readcounter $input > $output\n" if DEBUG;
`$readcounter $input > $output`;

