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

 ./launchHMMcopy.pl --r-libdir Path/to/directory/with/R.modules --output output.wig --read-counter /bin/HMMcopy/bin/readCounter

 will write into output.wig

=cut

my $USAGE = "launchHMMcopy.pl --r-libdir [Path to directory with R modules] --wig-normal [input normal wig] --wig-tumor [input tumor wig] --cg-file [cg file] --map-file [map file] --output-base [basename for output]";

# Required parameters
my($hmm_script,$rlibdir,$wig_normal,$wig_tumor,$cg_file,$map_file,$output_base);
my $results = GetOptions( "r-libdir=s"      =>  \$rlibdir,
                          "normal-wig=s"    =>  \$wig_normal,
                          "tumor-wig=s"     =>  \$wig_tumor,
                          "cg-file=s"       =>  \$cg_file,
                          "map-file=s"      =>  \$map_file,
                          "hmm-script=s"    =>  \$hmm_script,
                          "output-base=s"   =>  \$output_base);

if (!$wig_normal || !$wig_tumor || !$rlibdir || !$cg_file || !$map_file || !$output_base) { die $USAGE; }

# Make output directory if it does not exists                            
my $outdir = dirname($output_base);
`mkdir -p $outdir` if (!-e $outdir || !-d $outdir);

$ENV{R_LIBS} = $rlibdir;

#=====================================
# RUNNING HMMcopy script here
#=====================================
print STDERR "Running Rscript $hmm_script $wig_normal $wig_tumor $cg_file $map_file $output_base\n" if DEBUG;
`Rscript $hmm_script $wig_normal $wig_tumor $cg_file $map_file $output_base`;

