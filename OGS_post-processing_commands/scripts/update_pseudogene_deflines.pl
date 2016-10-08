#!usr/bin/perl -w
use strict;

#fasta file of cDNA sequences (exported from WA)
my $fasta = shift @ARGV or die;
#conversion file from Dan; format: 13175F52BB45029D6C64186F724BA8AD        OFAS027001
my $convert = shift @ARGV or die;
#gff file to get pseudogenic_transcript IDs for conversion
my $gff = shift @ARGV or die;
my $out = shift @ARGV or die;

open ( my $CONVERT, $convert ) or die;
my %convert = ();
while ( my $line = <$CONVERT> ){
    chomp $line;
    my @array = split "\t", $line;
    $convert{$array[0]} = $array[1];
}
close $CONVERT;

open ( my $GFF, $gff ) or die;
my %new_convert = ();
while ( my $line = <$GFF> ){
    chomp $line;
#ignore commented lines                                                                                                                                     
    if ( $line =~ /^#/ ){
        next;
    }
#if there's a fasta section at the end of the gff3, skip to the end of the file                                                                             
    if ( $line =~ /^>/ ){
        last;
    }
    else {
        my @array = split /\t/, $line;
        my @col9 = split ';', $array[8];
        my $id;                                                                                                                                          
	my $parent="NA";                                                                                                                                 
	foreach my $element (@col9){                                                                                                                         
	    if ( $element =~ /ID=(.*)/ ){
                $id = $1;
            }
	    elsif ( $element =~ /Parent=([A-Za-z0-9]{32})$/ ){
		$parent = $1;
	    }
	}
	if ( defined $convert{$parent} ){
#new conversion hash = WA transcript ID, DAn ID plus '-RA' suffix
	    $new_convert{$id} = $convert{$parent}."-RA";
	}
    }
}
close $GFF;

open ( my $FASTA, $fasta ) or die;
open ( my $OUT, ">$out" ) or die;
while ( my $fline = <$FASTA> ){
    chomp $fline;
    if ( $fline =~ /^>(\S+)/ ){
	if ( defined $new_convert{$1} ){
	    print $OUT ">$new_convert{$1}\n";
	}
	else {
	    warn "Defline $1 is not defined in conversion file\n";
	}
    }
    else {
	print $OUT "$fline\n";
    }
}
close $FASTA;
close $OUT;
exit;
