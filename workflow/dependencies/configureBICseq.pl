#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use constant DEBUG=>1;

my $USAGE = "configureBICseq.pl --input-t [tumor input] --input-n [normal input] --outdir [root data dir] --config-file [name of config file] --samtools [path to modified samtools]\n";
my($input_n,$input_t,$datadir,$config,$samtools);
my $results = GetOptions ("input-n=s"  => \$input_n,
                          "input-t=s"  => \$input_t,
                          "outdir=s"   => \$datadir,
                          "config-file=s"=>\$config,
                          "samtools=s" => \$samtools);
if (!$input_t || !$input_n || !$datadir || !$config || !$samtools){die $USAGE;}
my @inputs = ($input_n, $input_t);
my %info   = ();
my @chroms = ();
map{&process_bam} @inputs;

#=====================================
# Configuring BICseq run here
#=====================================
my @chromlines = grep {/^\@SQ/} `$samtools view -H $input_n`;
# TODO handle ref-free experiments
print STDERR "Got ".scalar(@chromlines)." chromosomes from bam file\n" if DEBUG;

my $tumorSeq  = basename($input_t,(".bam")).".seq";
my $normalSeq = basename($input_n,(".bam")).".seq";

print STDERR "Got $tumorSeq and $normalSeq for analysis\n" if DEBUG;

open(CONF,">$datadir/$config") or die "Could not create config file for BICseq";
print CONF join("\t",("chrom","tumor","normal"))."\n";
foreach my $line (@chromlines) {
  if ($line!~/SN\:(\S+)/){next;}
    my $c = $1;
    print STDERR "Found chromosome $c\n" if DEBUG;
    print CONF join("\t",($c,$tumorSeq,$normalSeq))."\n";
}

close CONF;

#======================================
# Subroutine for parsing .bam into .seq
#======================================
sub process_bam {
 my $bam = shift @_;

 if (!-e $bam && !-s $bam) {
  return;
 }
 
 my $pg = grep {/^\@PG/} `$samtools view -H $bam`;
 if ($pg && $pg=~/ID:(\S+)/) {
   my $aligner = $1;
   if ($aligner=~/BWA/) {
     `$samtools view -U BWA,$datadir,N,N $bam`;
   } elsif ($aligner=~/Bowtie/) {
     `$samtools view -U Bowtie,$datadir,N,N $bam`;
   } else {
     print STDERR "BICseq supports Bowtie and BWA aligners only, will terminate\n";
     exit;
   }
 }
}

