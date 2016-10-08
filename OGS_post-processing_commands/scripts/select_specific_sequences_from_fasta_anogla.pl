#!usr/bin/perl -w
use strict;

my $fasta = shift @ARGV or die;
#$list = list of sequence IDs that you want to retrieve from fasta file
my $list = shift @ARGV or die;
my $out = shift @ARGV or die;

my %list = ();
open ( my $LIST, $list ) or die;
while ( my $line = <$LIST> ){
    chomp $line;
    my $newline = $line."-RA";
    $list{$newline} = $line;
}

open ( my $FASTA, $fasta ) or die;
open ( my $OUT, ">$out" ) or die;
my $counter = 0;
while ( my $fline = <$FASTA> ){
    chomp $fline;
    if ( $fline =~ /^>(\S+)/ ){
	if ( defined $list{$1} ){
	    $counter++;
	    print $OUT "$fline\n";
	}
	else {
	    $counter = 0;
	}
    }
    else {
	if ( $counter > 0 ){
	    print $OUT "$fline\n";
	}
    }
}
