#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($Bin);
use constant DEBUG=>0;

=head2 Running FREEC

 FREEC is a tool for analyzing copy-number variation
 and may run with WholeGenome and Exome data (currently this script supports only WG)
 bare-bone run:

 ./launchFREEC.pl --input-normal normal.bam --input-tumor tumor.bam --freec /bin/freec --id someID --lenfile hg19.len --samtools /bin/samtools --r-libdir My/R_modules

 will write into ./data/FREEC_someID directory by default

=cut

my $USAGE = "launchFREEC.pl --input-tumor [tumor input] --input-normal [normal input] --data-type [optional, default is WG, EX supported]".
            " --r-libdir [R lib dir] --outdir [root data dir] --config-file [name of template config file] --samtools [path to samtools] --prefix [prefix for result files]\n";

# Required parameters
my($input_n,$input_t,$type,$id,$config,$lenfile,$samtools,$freec,$prefix);
# Optional parameters
my($datadir,$ploidy,$makebedgraph,$matetype,$logfile,$targetFile,$cvar,$quiet,$window,$rlibdir);
my $results = GetOptions ("input-normal=s" =>  \$input_n,
                          "input-tumor=s"  =>  \$input_t,
                          "r-libdir=s"     =>  \$rlibdir,
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
                          "prefix=s"       => \$prefix,
                          "config-file=s"  =>\$config,
                          "ploidy=s"       =>\$ploidy,
                          "make-bedgraph=s"=> \$makebedgraph,
                          "mate-type=s"    =>\$matetype);

if (!$input_t || !$input_n || !$samtools || !$freec || !$id || !$rlibdir) { die $USAGE; }
if ($type && $type ne "WG" && (!$targetFile || !-e $targetFile)) { die "Exome data passed, but no target .bed provided!"; }

$ENV{R_LIBS} = $rlibdir;

# Set defaults
$datadir  ||= "data/";
$config   ||= "freec.".$id.".conf";
$ploidy   ||= 2;
$makebedgraph ||= "TRUE";
$matetype ||= &guessMate($input_n);
$type     ||= "WG";
$cvar     ||= "0.05";
$logfile  = "freec.".$id.".log";

$datadir.="/" if $datadir!~m!/$!;
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

if ($type ne "WG") {
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

=head2
 
 Adding Wilcoxon and Kolmogorov-Smirnov test p-values to identified CNVS
 File will be named [$input_t]_CNVs.p.value.txt
 
 First we need to find if there are inputs

=cut

my @files = `find $datadir -type f`;
my $freec_dir = undef;

if (@files && @files > 1) {
 $freec_dir = dirname($files[0]);
 print STDERR "Got [$freec_dir] for FREEC output dir\n" if DEBUG; 
}
 my $tumor_f   = fileparse($input_t);
 my $ratioFile = $tumor_f."_ratio.txt";
 my $cnvFile   = $tumor_f."_CNVs";

 if ($freec_dir) {
   print STDERR "Adding significance values to called variations in [$freec_dir/$ratioFile]\n";
   print STDERR "Command: Rscript $Bin/assess_significance.R $freec_dir/$ratioFile $freec_dir/$cnvFile\n" if DEBUG;
   `Rscript $Bin/assess_significance.R $freec_dir/$ratioFile $freec_dir/$cnvFile`;
 }


=head2
 
 We can also visualize the results using makeGraph.R 
 script removing NA first

=cut

 my $noNAfile = $ratioFile;
 $noNAfile=~s/_ratio.txt/_ratio_noNA.txt/;

 if (-e "$freec_dir/$ratioFile") {
  
  open(IN,"<$freec_dir/$ratioFile") or die "Couldn't read from ratio file";
  open(OUT,">$freec_dir/$noNAfile") or die "Couldn't write to output [$freec_dir/$noNAfile]";
  while (<IN>) {
   my @temp = split("\t");
   if ($temp[2] eq "-1"){next;}
   print OUT $_;
  }
 
  close IN;
  close OUT;

  print STDERR "Making CNV visualization for [$freec_dir/$ratioFile]\n" if DEBUG;
  `Rscript $Bin/makeGraph.R $ploidy $freec_dir/$noNAfile`;
 
 }

=head2

 It looks like FREEC cannot produce custom names, so
 if we want to prefix our files with something like freec_
 this is where we do it

=cut

 &prefixFiles($freec_dir, $prefix);

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
 my @flaglines = grep {/paired/} `$samtools flagstat $bam`;
 if ($flaglines[0] =~/^0/) {return 0;}
 
 # Return FR (Illumina paired-end sequencing, SOLID is 'SS' and Illumina Mate Pair data 'RF'
 # this script will only guess between two options - single-end (0) or paired-end (FR)
 return "FR";

}

#==================================================
# Subroutine for prefixing all files in a directory
#==================================================

sub prefixFiles {

 my($dir,$pr) = @_;
 return if !$pr || !-d $dir;

 opendir(DIR,$dir) or die "Couldn't read from directory [$dir]";
 my @files = readdir(DIR);
 closedir DIR;

 #actual renaming
 foreach my $file (@files) {
   next unless -e "$dir/$file" && -f "$dir/$file";
   my $prefixed_file = $pr.$file;
   `mv $dir/$file $dir/$prefixed_file`;
 } 
 
}
