#!/usr/bin/perl -w
use strict;

my $wa = shift @ARGV or die;
my $dan = shift @ARGV or die;
my $out = shift @ARGV or die;

my %wa = ();
my %exon_hash = ();
open ( my $WA, $wa ) or die;
while ( my $wline = <$WA> ){
    chomp $wline;
#ignore commented lines                                                                                                        
    if ( $wline =~ /^#/ ){
        next;
    }
#if there's a fasta section at the end of the gff3, skip to the end of the file                                                  
    if ( $wline =~ /^>/ ){
        last;
    }
    else {
        my @warray = split /\t/, $wline;
        my $wtype = $warray[2];
	my $wscaf = $warray[0];
	my $wstart = $warray[3];
	my $wstop = $warray[4];
	my $wdir = $warray[6];
	my $wid;
	my $wparent;
	my @wcol9 = split ";", $warray[8];
	foreach my $welement ( @wcol9 ){
	    if ( $welement =~ /ID=(.*)/ ){
		$wid = $1;
	    }
	    elsif ( $welement =~ /Parent=(.*)/ ){
                $wparent = $1;
            }
	}
	my $wprint = "$wscaf\t$wdir\t$wstart\t$wstop";
#	if ( $wtype =~ /gene|pseudogenic_transcript|RNA/ ){
       if ( $wtype =~ /pseudogenic_transcript|RNA/ ){                                                                                                                                
#	    if ( defined $wa{$wprint} ){
#		warn "multiple features in file $wa ($wid, $wa{$wprint}) with coordinates $wprint\n";
#	    }
#	    else {
		$wa{$wid} = $wprint;
#	    }
	}
	elsif ( $wtype =~ /exon|CDS/ ){
	    $exon_hash{$wparent} .= "\t$wstart-$wstop";
	}
    }
}

foreach my $ekey ( keys %exon_hash ){
    my @array = split '\t', $exon_hash{$ekey};
    my @sort = sort @array;
    my $newvalue = join '\t', @sort;
    $exon_hash{$ekey} = $newvalue;
    if ( defined $wa{$ekey} ){
	$wa{$ekey} .= $newvalue;
    }
    else {
	print "missing parent in wa hash for ID $ekey\n";
    }
}

my %dan = ();
my %dexon_hash = ();
open ( my $DAN, $dan ) or die;
while ( my $line = <$DAN> ){
    chomp $line;
    if ( $line =~ /^#/ ){
	next;
    }
    my @array = split /\t/, $line;
    my $type = $array[2];
    my $scaf = $array[0];
    my $start = $array[3];
    my $stop = $array[4];
    my $dir = $array[6];
    my $id;
    my $parent;
    my @col9 = split ";", $array[8];
    foreach my $element ( @col9 ){
	if ( $element =~ /ID=(.*)/ ){
	    $id = $1;
	}
	elsif ( $element =~ /Parent=(.*)/ ){
            $parent = $1;
        }
    }
    my $print = "$scaf\t$dir\t$start\t$stop";
#    if ( $type =~ /gene|pseudogenic_transcript|RNA/ ){
    if ( $type =~ /pseudogenic_transcript|RNA/ ){
#	if ( defined $dan{$print} ){
#	    warn "multiple features in file $dan ($id, $dan{$print}) with coordinates $print\n";
#	}
#	else {
	    $dan{$id} = $print;
#	}
    }
    elsif ( $type =~ /exon|CDS/ ){
	$dexon_hash{$parent} .= "\t$start-$stop";
    }
}

foreach my $dkey ( keys %dexon_hash ){
    my @array = split '\t', $dexon_hash{$dkey};
    my @sort = sort @array;
    my $newvalue = join '\t', @sort;
    $dexon_hash{$dkey} = $newvalue;
    if ( defined $dan{$dkey} ){
	$dan{$dkey} .= $newvalue;
    }
    else {
	print "missing parent in dan hash for ID $dkey\n";
    }
}

my %rev_wa = reverse %wa;
my %rev_dan = reverse %dan;

my $num_wa = scalar keys %wa;
my $num_revwa = scalar keys %rev_wa;
my $num_dan = scalar keys %dan;
my $num_revdan = scalar keys %rev_dan;

print "wa: $num_wa, reverse wa: $num_revwa, dan: $num_dan, reverse dan: $num_revdan\n";


#output format: old ID\t new ID
open ( my $OUT, ">$out" ) or die;
foreach my $key ( keys %rev_wa ){
    if ( defined $rev_dan{$key} ){
	print $OUT "$rev_wa{$key}\t$rev_dan{$key}\n";
    }
    else {
	print "WA entry $rev_wa{$key}\t$key not defined in Dan's file\n";
    }
}

my %test = ();
foreach my $key2 ( keys %wa ){
    if ( defined $test{$wa{$key2}} ){
	print "$key2 and $test{$wa{$key2}} are have identical values\n";
    }
    else {
	$test{$wa{$key2}} = $key2;
    }
}

my %test2 = ();
foreach my $key3 ( keys %dan ){
    if ( defined $test2{$dan{$key3}} ){
        print "$key3 and $test2{$dan{$key3}} are have identical values\n";
    }
    else {
	$test2{$dan{$key3}} = $key3;
    }
}
