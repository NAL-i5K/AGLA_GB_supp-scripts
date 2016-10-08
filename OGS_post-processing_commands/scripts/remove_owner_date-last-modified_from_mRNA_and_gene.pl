#usr/bin/perl -w
use strict;

#script searches for gene ID; if there's a match to this ID in the list hash, then finds owner attribute and removes from line; prints line to stdout. If no matches, prints line to stdout. Prints 'removed ...' statement to stderror if owner attribute is removed

my %list = ('AGLA014663', '1', 'AGLA003801', '1', 'AGLA017751', '1', 'AGLA003809', '1','AGLA000919', '1','AGLA000103', '1');
my $gff = shift @ARGV or die;
my $out = shift @ARGV or die;

open ( my $GFF, $gff ) or die;
open ( my $OUT, ">$out" ) or die;
while ( my $line = <$GFF> ){
    chomp $line;
    my $id;
    if ( $line =~ /ID=(\w+\d+)-*\w*;*/ ){
	$id = $1;
    }
    if ( defined $list{$id} ){
	if ( $line =~ /owner=\w+;/ ){
	    $line =~ s/owner=\w+;//;
	    print "removed owner attribute from the following line:\n$line\n";
	}
	elsif ( $line =~ /owner=\w+/ ){
            $line =~ s/owner=\w+//;
            print "removed owner attribute from the following line:\n$line\n";
	}
	if ( $line =~ /date_last_modified=\d+-\d+-\d+;/ ){
            $line =~ s/date_last_modified=\d+-\d+-\d+;//;
            print "removed date_last_modified attribute from the following line:\n$line\n";
        }
	elsif ( $line =~ /date_last_modified=\d+-\d+-\d+/ ){
            $line =~ s/date_last_modified=\d+-\d+-\d+//;
            print "removed date_last_modified attribute from the following line:\n$line\n";
        }
#date_last_modified=2014-04-20
    }
    print $OUT "$line\n";
}

close $GFF;
close $OUT;
exit;
