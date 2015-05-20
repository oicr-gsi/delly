#!/usr/bin/perl -w

# =========================================================
# Varscan2 wrapper - wrapping Varscan tool for CNV workflow
# =========================================================

=head2 Initial Data processing
 
 First, we need to produce pileup file using normal and tumor inputs
 It needs to be additionally parsed (using system awk) to remove zero-coverage
 entries. Samtools removes low-quality reads

 samtools mpileup -q 1 -f hg18.fa normal_sorted.bam tumor_sorted.bam | awk -F"\t" '$4 > 0 && $7 > 0' > normtumor_sorted.pileup

=cut

use strict;
use File::Basename;
use File::Path;
use Getopt::Long;
use constant DEBUG=>0;

my($normal, $tumor, $out_dir, $rlibs_dir, $ref_fasta, $samtools, $java, $varscan, $id);
my $USAGE = "launchVarscan2.pl --input-normal [normal .bam file] --input-tumor [tumor .bam file] --output-dir [Output directory] --ref-fasta [Reference fasta] --samtools [Path to samtools] --rlibs-dir [path to directory with R libs] --id [unique id]";
my $result = GetOptions('input-normal=s' => \$normal,
                        'input-tumor=s'  => \$tumor,
                        'output-dir=s'   => \$out_dir,
                        'ref-fasta=s'    => \$ref_fasta,
                        'java=s'         => \$java,
                        'varscan=s'      => \$varscan,
                        'id=i'           => \$id,
                        'rlibs-dir=s'    => \$rlibs_dir,
                        'samtools=s'     => \$samtools);

if (!$normal || !$tumor || !$out_dir || !$ref_fasta || !$samtools || !$rlibs_dir) {die $USAGE;}
if (!-e $out_dir || !-d $out_dir) {
 mkpath($out_dir);
}

map{if(!/\.bam$/){die "Files are not in .bam format!";}} ($normal,$tumor);
my $pileup_command = "$samtools mpileup -q 1 -f $ref_fasta $normal $tumor | awk -F \"\t\" '\$4 > 0 && \$7 > 0' > $out_dir/normtumor_sorted.pileup";
print STDERR "Will run command $pileup_command" if DEBUG;

`$pileup_command`;
print STDERR "Pileup data produced, starting actual analysis...\n" if DEBUG;

=head2 Running analysis

 Varscan would produce just one file, but some modulation of min-coverage parameter
 may be needed. The default is 20, but the script will attempt to lower this parameter
 if Varscan fails to produce results

=cut

my $resultsOK = undef;

if (-e "$out_dir/normtumor_sorted.pileup" && -s "$out_dir/normtumor_sorted.pileup") {
  my $cvg = undef;
  do {
      my $varscan_command = "$java -jar $varscan copynumber $out_dir/normtumor_sorted.pileup $out_dir/varscan_out.$id -mpileup 1";
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

$ENV{R_LIBS} = $rlibs_dir;

print STDERR "Will run CBS script to reduce noise in data\n" if DEBUG;
`Rscript smooth_varscan.r $out_dir/varscan_out.$id.copynumber`;

