#!/usr/bin/perl -w

=head2 Varscan2 wrapper - wrapping Varscan tool for CNV workflow

 Initial Data processing
 
 First, we need to produce pileup file using normal and tumor inputs
 It needs to be additionally parsed (using system awk) to remove zero-coverage
 entries. Samtools removes low-quality reads

 it is recommended to pass sample name as id b/c it will be used in plots as a title

 samtools mpileup -q 1 -f hg18.fa normal_sorted.bam tumor_sorted.bam | awk -F"\t" '$4 > 0 && $7 > 0' > normtumor_sorted.pileup

 The command will be run according to parameters supplied by user, will look like this:

 java -Xmx4G -jar VarScan.jar copynumber normtumor_sorted.pileup output_basename -mpileup 1 --p-value 0.01 --min-coverage 10

 Note that user may specify p-value (threshold) but the default is set to 0.05

=cut

use strict;
use File::Basename;
use File::Path;
use FindBin qw($Bin);
use Getopt::Long;
use constant DEBUG=>0;

my($normal, $tumor, $out_dir, $rlibdir, $ref_fasta, $samtools, $java, $varscan, $id, $min_coverage, $del_threshold, $min_region_size, $recenter_up, $recenter_down, $pvalue, $xmxmem);
my $USAGE = "launchVarscan2.pl --input-normal [normal .bam file] --input-tumor [tumor .bam file] --output-dir [Output directory] --ref-fasta [Reference fasta] --samtools [Path to samtools] --r-libdir [path to directory with R libs] --id [unique id]";
my $result = GetOptions('input-normal=s'    => \$normal,
                        'input-tumor=s'     => \$tumor,
                        'output-dir=s'      => \$out_dir,
                        'ref-fasta=s'       => \$ref_fasta,
                        'java=s'            => \$java,
                        'varscan=s'         => \$varscan,
                        'xmxmem=i'          => \$xmxmem,
                        'id=s'              => \$id,
                        'r-libdir=s'        => \$rlibdir,
                        'samtools=s'        => \$samtools,
                        'p-value=s'         => \$pvalue,
                        'min-coverage=s'    => \$min_coverage,
                        'del-coverage=s'    => \$del_threshold,
                        'min-region-size=s' => \$min_region_size,
                        'recenter-up=s'     => \$recenter_up,
                        'recenter-down=s'   => \$recenter_down);

if (!$varscan || !$normal || !$tumor || !$out_dir || !$ref_fasta || !$samtools || !$rlibdir ) {die $USAGE;}
if (!-e $out_dir || !-d $out_dir) {
 mkpath($out_dir);
}
$pvalue ||= 0.05; # Default threshold for p-value
$xmxmem ||= 6;

map{if(!/\.bam$/){die "Files are not in .bam format!";}} ($normal,$tumor);
my $pileup_command = "$samtools mpileup -q 1 -f $ref_fasta $normal $tumor | awk -F \"\t\" '\$4 > 0 && \$7 > 0' > $out_dir/normtumor_sorted.pileup";
print STDERR "Will run command $pileup_command\n" if DEBUG;

`$pileup_command`;
print STDERR "Pileup data produced, starting actual analysis...\n" if DEBUG;

$ENV{R_LIBS} = $rlibdir;

=head2 Running analysis

 Varscan would produce just one file, but some modulation of min-coverage parameter
 may be needed. The default is 20, but the script will attempt to lower this parameter
 if Varscan fails to produce results

=cut

my $resultsOK = undef;

if (-e "$out_dir/normtumor_sorted.pileup" && -s "$out_dir/normtumor_sorted.pileup") {
  my $cvg = undef;
  do {
      my $varscan_command = "$java -Xmx".$xmxmem."G -jar $varscan copynumber $out_dir/normtumor_sorted.pileup $out_dir/$id -mpileup 1 --p-value $pvalue";
      $varscan_command .= " --min-coverage $cvg" if $cvg;
      print STDERR "Command: $varscan_command\n" if DEBUG; 
      $result = `$varscan_command 2>&1`;
      print STDERR "Got results: [$result]\n" if DEBUG;
      $cvg ||=20;
      $cvg-=4; 		# This is our increment
 
     if ($result=~/(\d+) had sufficient coverage/) {
        if ($1 == 0) {
          if ($cvg < 2) {
            print STDERR "Unable to run Varscan even with min-coverage set to $cvg, aborting...\n" if DEBUG;
            last;
          }
          print STDERR "Coverage threshold too high, trying min-coverage $cvg...\n" if DEBUG;
        } else {
          $resultsOK = 1;
          if ($min_coverage > $cvg) {
              $min_coverage = $cvg;
          }
        }
      }
    
  } while (!$resultsOK && $cvg >= 2);
  
  # Log failed analysis
  my $timeStamp = time;
  `echo $timeStamp > $out_dir/VARSCAN_FAILED` if !$resultsOK;  

} else {
  die "Pileup file is not suitable for analysis, will exit now";
}


=head2 Additional Smoothing

 Circular Binary Segmentation aims to improve (potentially) noisy data
 Requires Bioconductors library DNAcopy. Smoothed data will be written
 into the same directory as input, file name will have .segmented 
 suffix added

=cut 

if (-e "$out_dir/$id.copynumber") {

 print STDERR "Will run Additional filtering on results\n" if DEBUG;
 
 my $filterCommand = "$java -Xmx".$xmxmem."G -jar $varscan copyCaller $out_dir/$id.copynumber --output-file $out_dir/$id.copynumber.filtered";

 # Include filtering options, if set
 if ($min_coverage) {
     $filterCommand.=" --min-coverage $min_coverage";
 }

 if ($del_threshold) {
     $filterCommand.=" --del-threshold $del_threshold";
 }

 if ($min_region_size) {
     $filterCommand.=" --min-region-size $min_region_size";
 }

 if ($recenter_up) {
     $filterCommand.=" --recenter-up $recenter_up";
 } 

 if ($recenter_down) {
     $filterCommand.=" --recenter-down $recenter_down";
 }

 `$filterCommand`;
 
 $ENV{R_LIBS} = $rlibdir;

 print STDERR "Will run CBS script to reduce noise in data\n" if DEBUG;
 my $message = `Rscript $Bin/smooth_varscan.r $out_dir/$id.copynumber $id`;
 print STDERR $message;
 $message    = `Rscript $Bin/smooth_varscan.r $out_dir/$id.copynumber.filtered $id.filtered`;
} else {
 print STDERR "File with copynumber calls does not exist, won't attempt CBS smoothing/visualization\n";
}
