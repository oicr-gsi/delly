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
use constant DEBUG=>1;

my($normal, $tumor, $out_dir, $rlibs_dir, $ref_fasta, $samtools, $java, $varscan);
my $USAGE = "launchVarscan2.pl --input-normal [normal .bam file] --input-tumor [tumor .bam file] --output-dir [Output directory] --ref-fasta [Reference fasta] --samtools [Path to samtools] --rlibs-dir [path to directory with R libs]";
my $result = GetOptions('input-normal=s' => \$normal,
                        'input-tumor=s'  => \$tumor,
                        'output-dir=s'   => \$out_dir,
                        'ref-fasta=s'    => \$ref_fasta,
                        'java=s'         => \$java,
                        'varscan=s'      => \$varscan,
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

=head2 Running analysis

 Varscan would produce just one file, but some modulation of min-coverage parameter
 may be needed. The default is 20, but the script will attempt to lower this parameter
 if Varscan fails to produce results

=cut

my $resultsOK = undef;

if (-e "$out_dir/normtumor_sorted.pileup" && -s "$out_dir/normtumor_sorted.pileup") {
  my $cvg = undef;
  do {
      my $varscan_command = "$java -jar $varscan copynumber $out_dir/normtumor_sorted.pileup $out_dir/varscan_out -mpileup 1";
      $varscan_command .= " --min-coverage $cvg" if $cvg;
    
      $result = `$varscan_command`;
      if ($result=~/(\d+) had sufficient coverage/) {
        if ($1 == 0) {
          $cvg ||=20;
          $cvg-=4; # This is our increment
        } else {
          $resultsOK = 1;
        }
      }
    
  } while (!$resultsOK || $cvg >= 2);
  
  # Log failed analysis
  my $timeStamp = time;
  `echo $timeStamp > $out_dir/VARSCAN_FAILED` if !$resultsOK;  

} else {
  print STDERR "Pileup file is not suitable for analysis, will exit now\n";
}

=head2 Additional Smoothing

 Circular Binary Segmentation aims to improve (potentially) noisy data
 Requires Bioconductor and library DNAcopy

=cut 

$ENV{R_LIBS} = $rlibs_dir;

