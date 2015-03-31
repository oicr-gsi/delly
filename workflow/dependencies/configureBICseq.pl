#!/usr/bin/perl -w

use strict;
use Getopt:Long;
use constant DEBUG=>1;

my $USAGE = "configureBICseq.pl --input-t [tumor input] --input-n [normal input] --outdir [root data dir] --bicseq-interval [I] --bicseq-spread [S] --samtools [path to modified samtools]\n";
my($input_n,$input_t,$output,$samtools);
my $results = GetOption("--input-n=s"  => \$input_n,
                        "--input-t=s"  => \$input_t,
                        "--outdir=s"   => \$datadir,
                        "--bicseq-interval=i" => \$bicseq_interval,
                        "--bicseq-spread=i" => \$bicseq_spread,
                        "--samtools=s" => \$samtools);

my @inputs = ($input_n, $input_t);
my %info = ();
map{&process_bam} @inputs;








#
# Subroutine for parsing .bam into .seq
#

sub process_bam {
 my $bam = shift @_;
 my $seqfile = $bam;
 $seqfile=~s/.bam$/.seq/;

 if (!-e $bam && !-s $bam) {
  return;
 }
 
 $pg = grep {/^\@PG/} `$samtools view -H $bam`;
 if ($pg && $pg=~/ID:(\S+)/) {
   my $aligner = $1;
   if ($aligner=~/BWA/) {
     `$samtools view -U BWA,$datadir,`;
   } elsif ($aligner=~/Bowtie/) {
     ``;

   } else {


   }
 }
}
