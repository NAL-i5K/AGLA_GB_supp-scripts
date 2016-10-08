#!/usr/bin/perl -w

use strict;

if( scalar ( @ARGV ) != 2 ){
	my @temp=split(/\//, $0);
	die
"

USAGE: $temp[-1] <list of sequence names> <fasta file>

This script reads a fasta file, and extracts the reads whose names appear in the small file

";
}

open ( NAMES, $ARGV[0] ) or die;

my %names_to_extract = ();

# this will hold the number of sequences to be extracted
my $seq_nr = 0;

while ( my $one_name = <NAMES> ) {
	
	chomp ($one_name);
	if ( !$names_to_extract{$one_name} ) {
		$names_to_extract{$one_name} = 1;
		$seq_nr++;
	}
}



open ( FASTA, $ARGV[1] ) or die;

# this will hold the number of extracted sequences
my $seq_ext = 0;

my $line = <FASTA>;

while ( defined($line) ) {
	if ( $line !~ /^>/ ) {
		$line = <FASTA>;
	}
	
	if ( $line =~ /^>/ ) {
		my $header = $line;
		
		my $name = $header;
		$name =~ s/^>(\S+).*$/$1/g;
		chomp ($name);

		my $sequence = "";
		$line = <FASTA>;
		until ( !defined ($line) || $line =~ /^>/ ) {
			$sequence .= $line;
			$line = <FASTA>;
		}
		
		foreach my $tmp (keys%names_to_extract) {
			if ( $header =~ /$tmp/ ) {
				print $header,$sequence;
				$seq_ext++;
			}
		}
	}
}

# see if you didn't find some sequences
if ( $seq_nr != $seq_ext ) {
	print STDERR "Didn't find ", ( $seq_nr - $seq_ext ) ," sequences\n";
}

