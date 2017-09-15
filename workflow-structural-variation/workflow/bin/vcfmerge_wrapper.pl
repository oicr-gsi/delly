#!/usr/bin/perl -w

# =====================================================
# wrapper for vcf-merge tool, meant to bypass 
# output as giant matrix (Version for workflow)
# =====================================================

use Getopt::Long;
use Data::Dumper;
use strict;


use constant DEBUG =>0;
# Below is the dafault for vcf_compare, should not be used when workflow runs

my($list,$vcf_merge,$outfile,$datadir,$path_to_tabix);
my $USAGE="vcfmerge_wrapper.pl --list=[req] --datadir=[optional] --vcf-merge=[req]";

my $result = GetOptions ('list=s'            => \$list, # list with filenames
                         'datadir=s'         => \$datadir, # directory with vcf files
                         'tabix=s'           => \$path_to_tabix, # path to tabix dir
                         'output=s'          => \$outfile, #output file
                         'vcf_merge=s'       => \$vcf_merge); # path to vcftools vcf-merge script

if (!$list || !$vcf_merge || !$datadir || !$path_to_tabix) {die $USAGE;}
$datadir ||=".";

#===========================================================
# If tabix is not in the PATH, add its location to $PATH
#===========================================================
my $PATH = `echo \$PATH`;
my $reset_path;
my $tabix_check = `which tabix`;

if (!$tabix_check) {
 print STDERR "tabix not found, adding [$path_to_tabix] to the PATH...\n" if DEBUG;
 $ENV{PATH} = "$path_to_tabix:$PATH";
 $reset_path = 1;
} else {
 print STDERR "Found tabix, proceeding...\n" if DEBUG;
}

chomp($list);
my @listed = grep {/\S+/} split(" ",$list);
my @existing = ();
map{push @existing, $_ if (-e $_ && -s $_)} (@listed);
print STDERR Dumper(@existing) if DEBUG;
my $vetted_list = join(" ",@existing);
#============================================================
# vcf-merge wrapping
#============================================================
if (@existing > 0) {
my $vetted_list = join(" ",@existing);
 my $c = $vcf_merge." ".$vetted_list." > ".$outfile;
 print STDERR "Command: $c\n" if DEBUG;
 `$c`;
}
