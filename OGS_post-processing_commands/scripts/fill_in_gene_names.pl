#!usr/bin/perl -w
use strict;

#purpose: most genes in the anogla OGS are without a name. This script finds these genes and sets the gene name equal to the gene ID. 

my $gff = shift @ARGV or die;
my $out = shift @ARGV or die;

open ( my $GFF, $gff ) or die;
open ( my $OUT, ">$out" ) or die;
while ( my $line = <$GFF> ){
    chomp $line;
    if ( $line =~ /^#/ ){
        print $OUT "$line\n";
    }
    else {
        my @array = split /\t/, $line;
        my @col9 = split ";", $array[8];
        my $type = $array[2];
	if ( $type =~ /gene/ ){
	    my $id;
	    my $name = "NA";
	    foreach my $element (@col9){
		if ( $element =~ /ID=(.*)/ ){
		    $id = $1;
		}
		elsif ( $element =~ /Name=(.*)/ ){
		    $name = $1;
		}
	    }
	    if ( $name =~ /^NA$/ ){
		if ( $id =~ /^NA$/ ){
		    warn "no ID defined in line:\n$line\n";
		}
		else {
#gene features without a name shouldn't have a name attribute, period, so you can add a new name attribute to the end of the line. Check for semi-colon at the end of line
		    if ( $line =~ /^.+;$/ ){
			$line = $line."Name=$id;";
		    }
		    else {
			$line = $line.";Name=$id";
                    }
		    print "added gene name to the following line:\n$line\n";
     		}
	    }
	    print $OUT "$line\n";
	}
	else {
	    print $OUT "$line\n";
	}
    }
}
close $GFF;
close $OUT;
