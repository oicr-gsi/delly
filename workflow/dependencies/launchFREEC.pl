#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use constant DEBUG=>0;

=head2 Running FREEC

 FREEC is a tool for analyzing copy-number variation
 and may run with WholeGenome and Exome data (currently this script supports only WG)
 bare-bone run:

 ./launchFREEC.pl --input-normal normal.bam --input-tumor tumor.bam --freec /bin/freec --id someID --lenfile hg19.len --samtools /bin/samtools

 will write into ./data/FREEC_someID directory by default

=cut

my $USAGE = "launchFREEC.pl --input-tumor [tumor input] --input-normal [normal input] --data-type [optional, default is WG, EX supported] ".
            " --outdir [root data dir] --config-file [name of template config file] --samtools [path to samtools]\n";

# Required parameters
my($input_n,$input_t,$type,$id,$config,$lenfile,$samtools,$freec);
# Optional parameters
my($datadir,$ploidy,$makebedgraph,$matetype,$logfile,$targetFile,$cvar,$quiet,$window);
my $results = GetOptions ("input-normal=s" =>  \$input_n,
                          "input-tumor=s"  =>  \$input_t,
                          "samtools=s"     =>  \$samtools,
                          "lenfile=s"      =>  \$lenfile,
                          "freec=s"        =>  \$freec,
                          "id=s"           =>  \$id,
                          # Optional parameters
                          "quiet=s"        =>  \$quiet,
                          "data-type=s"    =>  \$type,
                          "var-coefficient=s"=> \$cvar,
                          "window=s"       => \$window,
                          "target-file=s"  =>\$targetFile,
                          "outdir=s"       =>\$datadir,
                          "config-file=s"  =>\$config,
                          "ploidy=s"       =>\$ploidy,
                          "make-bedgraph=s"=> \$makebedgraph,
                          "mate-type=s"    =>\$matetype);

if (!$input_t || !$input_n || !$samtools || !$freec || !$id) { die $USAGE; }
if ($type && $type eq "EX" && (!$targetFile || !-e $targetFile)) { die "Exome data passed, but no target .bed provided!"; }

# Set defaults
$datadir  ||= "data/";
$config   ||= "freec.".$id.".conf";
$ploidy   ||= 2;
$makebedgraph ||= "TRUE";
$matetype ||= &guessMate($input_n);
$type     ||= "WG";
$cvar     ||= "0.5";
$logfile  = "freec.".$id.".log";

$datadir.="/" if $datadir!~m!/$!;
$datadir.="FREEC_".$id."/";
`mkdir -p $datadir` if (!-e $datadir || !-d $datadir);

# Default config file
my $configDefault = <<CONF;
[general]

chrLenFile = LENGHTFILE_TAG
BedGraphOutput = MAKEBED_TAG
coefficientOfVariation = CVAR_TAG
outputDir = OUTDIR_TAG
ploidy = PLOIDY_TAG
samtools = SAMTOOLS_TAG
window   = WINDOW_TAG

[sample]

mateFile = TUMOR_TAG
inputFormat = BAM
mateOrientation = MATE_TAG

[control]

mateFile = NORMAL_TAG
inputFormat = BAM
mateOrientation = MATE_TAG

CONF

my $exomeConf = <<EXOME;

[target]

captureRegions = TARGETFILE_TAG

EXOME

my $configBody = $configDefault;

if ($type eq "EX") {
    $configBody .= $exomeConf;
    $configBody =~s/TARGETFILE_TAG/$targetFile/;
    $window ||= 500;
} else {
    $window ||= 50000;
}


#else {
#    print STDERR "THIS IS NOT EX, undeff-ing window param\n" if DEBUG;
#    undef($window);
#}

if ($window) {
    print STDERR "Replacing window param, setting it to [$window]\n" if DEBUG;
    $configBody =~s/WINDOW_TAG/$window/;
} else {
    print STDERR "window param not set, removing tag\n" if DEBUG;
    $configBody =~s/.*WINDOW_TAG//;
}

$configBody=~s/LENGHTFILE_TAG/$lenfile/;
$configBody=~s/MAKEBED_TAG/$makebedgraph/;
$configBody=~s/CVAR_TAG/$cvar/;
$configBody=~s/OUTDIR_TAG/$datadir/;
$configBody=~s/PLOIDY_TAG/$ploidy/;
$configBody=~s/SAMTOOLS_TAG/$samtools/;
$configBody=~s/TUMOR_TAG/$input_t/;
$configBody=~s/NORMAL_TAG/$input_n/;
$configBody=~s/MATE_TAG/$matetype/g;

my $configPath = $datadir.$config;
print STDERR "# Will write into Config file $configPath \n" if DEBUG;
open(CONF, ">$configPath") or die "Cannot write to Config file [$configPath]";
print CONF $configBody;
close CONF;

=head2 Logging

 By default, we log output from freec, however option quiet will disable this

=cut
#=====================================
# RUNNING FREEC here
#=====================================
if (!-e $freec || !-e $input_t || !-e $input_n || !-e $lenfile) { die "Invalid parameters, make sure all files exist and run again!";}

# Depending on our wish to log output:
if (!$quiet) {
  my $logPath = $datadir.$logfile;
  print STDERR "Would launch Command [$freec --conf $configPath] with logging\n" if DEBUG;
  `$freec --conf $configPath >> $logPath` unless DEBUG; 
} else {
  print STDERR "Would launch Command [$freec --conf $configPath] without logging\n" if DEBUG;
  `$freec --conf $configPath  1>/dev/null` unless DEBUG;
}

=head2 Sequencing type guessing

 The script cannot intelligently tell the difference between mate-pair and paired-end sequencing, therefore we 
 rely on parameters passed by user. However, we can attempt a guess by looking at samtools flagstat data
 if no mate information has been passed to us

=cut
#================================================================
# Subroutine for guessing mate type by looking at flagstat output
#================================================================

sub guessMate {

 my $bam = shift @_;
 my @flaglines = grep {/paired/} `samtools flagstat $bam`;
 if ($flaglines[0] =~/^0/) {return 0;}
 
 # Return FR (Illumina paired-end sequencing, SOLID is 'SS' and Illumina Mate Pair data 'RF'
 # this script will only guess between two options - single-end (0) or paired-end (FR)
 return "FR";

}
