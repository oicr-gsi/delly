#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use constant DEBUG=>0;

my $USAGE = "configureBICseq.pl --input-tumor [tumor input] --input-normal [normal input] --outdir [root data dir] --config-file [name of config file] --samtools [path to modified samtools]\n";
my($input_n,$input_t,$datadir,$config,$samtools);
my $results = GetOptions ("input-normal=s"  => \$input_n,
                          "input-tumor=s"  => \$input_t,
                          "outdir=s"   => \$datadir,
                          "config-file=s"=>\$config,
                          "samtools=s" => \$samtools);
if (!$input_t || !$input_n || !$datadir || !$config || !$samtools){die $USAGE;}
$datadir.="/" if $datadir!~m!/$!;
my @inputs = ($input_n, $input_t);
my %files   = ();
my @chroms = ();
map{&process_bam($_)} @inputs;

#=====================================
# Configuring BICseq run here
#=====================================
my @chromlines = grep {/^\@SQ/} `$samtools view -H $input_n`;
# TODO handle ref-free experiments
print STDERR "Got ".scalar(@chromlines)." chromosomes from bam file\n" if DEBUG;

my $tumorSeq  = basename($input_t,(".bam"));
my $normalSeq = basename($input_n,(".bam"));

print STDERR "Got $tumorSeq and $normalSeq for analysis\n" if DEBUG;

open(CONF,">$datadir/$config") or die "Could not create config file for BICseq";
print CONF join("\t",("chrom","tumor","normal"))."\n";
foreach my $line (@chromlines) {
  if ($line!~/SN\:(\S+)/){next;}
    my $c = $1;
    print STDERR "Found chromosome $c\n" if DEBUG;
    print CONF join("\t",($c,$datadir.$tumorSeq."_".$c.".seq",
                             $datadir.$normalSeq."_".$c.".seq"))."\n";
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
 print STDERR "Processing [$bam]...\n" if DEBUG;
 my @pg = grep {/^\@PG/} `$samtools view -H $bam`;
 if (@pg && $pg[0]=~/ID:(\S+)/) {
   print STDERR "Making seq file for aligner $1...\n" if DEBUG;
   my $aligner = $1;
   my $bamSeq = basename($bam,(".bam"))."_";

   opendir(DIR,"$datadir") or die "Couldn't read from directory [$datadir]";
   my @bfiles = grep {/$bamSeq/} readdir(DIR);
   if (@bfiles > 0) {
     chop($bamSeq);
     print STDERR "File(s) exist, will not create seq files for $bamSeq\n" if DEBUG;
   }

   if ($aligner=~/BWA/i && !-e $bamSeq) {
     `$samtools view -U BWA,$datadir/$bamSeq,N,N $bam`;
   } elsif ($aligner=~/Bowtie/i && !-e $bamSeq) {
     `$samtools view -U Bowtie,$datadir/$bamSeq,N,N $bam`;
   } else {
     print STDERR "BICseq supports Bowtie and BWA aligners only, will terminate\n";
     exit;
   }
 }
}

