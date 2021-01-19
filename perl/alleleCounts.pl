#! /usr/bin/perl -w
# parse vcf file to get allele counts for jgil
# restrict to <= 3 alternative alleles (for a total of 4)
# ============================================================

use Getopt::Long;
use List::Util qw(first);
use warnings;

while (<>) {
  
  chomp $_;
  
  if ($_ =~ m/^\#\#/) { next; }
  if ($_ =~ m/^\#CHROM/) {
    
    @line = split /\t/, $_;
    @ids = @line[9..$#line];
    print "CHR POS REF ALT ", join(" ", @ids), "\n";
    next;
    
  }
  
  @line = split /\t/, $_;
  @alleles = split /,/, $line[4];
  
  if ($#alleles > 2) { next; }
  
  @format = split /:/, $line[8];
  $adidx = first { $format[$_] eq "AD" } 0..$#format;
  
  @ad = ();
  
  for (my $i = 0; $i <= $#line - 9; $i++) {
    
    @this = split /:/, $line[$i + 9];
    push(@ad, $this[$adidx]);
    
  }
  
  print join(" ", (@line[(0, 1, 3, 4)], @ad)), "\n";
  
}
