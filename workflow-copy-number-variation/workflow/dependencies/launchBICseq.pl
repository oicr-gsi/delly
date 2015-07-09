#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use constant DEBUG=>0;

# PERL_pipeline/BICseq_1.1.2/BIC-seq/BIC-seq.pl --I 150,20 /u/pruzanov/Data/CNVtools/BICseq/test1.config /scratch2/users/pruzanov/Data/CNVTOOLS/BIC-seq.hn.test1 \"ResultID\"
my $USAGE = "launchBICseq.pl --bicseq-interval [bicseq interval] --bicseq-spread [bicseq spread] --outdir [output dir] --config-file [name of config file] --bicseq [path to BicSeq] --result-id [unique result id]\n";
my($bicseqi,$bicseqs,$bicseq,$outdir,$config,$id,$samtools);
my $results = GetOptions ("outdir=s"          => \$outdir,
                          "bicseq-interval=s" => \$bicseqi,
                          "bicseq-spread=s"   => \$bicseqs,
                          "config-file=s"     => \$config,
                          "result-id=s"       => \$id,
                          "bicseq=s"          => \$bicseq);

if ( !$id || !$bicseqi || !$bicseqs || !$outdir || !$config || !$bicseq){die $USAGE;}
$outdir.="/" if $outdir!~m!/$!;
if (-e $outdir) {die "Directory exists, won't run BICseq";}
#=====================================
# Launch BICseq
#=====================================
$id =~s/ /_/g; #Remove spaces
my $bicseqCommand = "$bicseq --I $bicseqi\,$bicseqs $config $outdir $id";
print STDERR "Command is: $bicseqCommand\n" if DEBUG;
`$bicseqCommand`;

