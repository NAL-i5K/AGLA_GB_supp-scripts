#!usr/bin/perl -w                                                                                                       
use strict;

#find mismatches between isoform suffixes in name and ID
#put in hash
#change ID suffix
#update parent IDs of child features
#need to also ensure that all transcript IDs remain unique; current version assumes that isoforms within a gene model are all formatted similarly, but this is not necessarily the case. Added hash that checks for duplicate transcript ID entries. 

my $gff = shift @ARGV or die;
my $out = shift @ARGV or die;

#key: ID, value: new ID
my %mismatches = ();
my %typehash = ();
#added in transcript ID checker to ensure that all transcript IDs are unique in modified gff
my %transcript_id_hash = ();
open ( my $GFF, $gff ) or die;
open ( my $OUT, ">$out" ) or die;
while ( my $line = <$GFF> ){
    chomp $line;
    if ( $line =~ /^#/ ){
        print $OUT "$line\n";
    }
    else {
        my @array = split /\t/, $line;
	my $col9 = pop @array;
#rejoin truncated line for printing later                                                                                  
        my $line2 = join "\t", @array;
        my @col9 = split ";", $col9;
	my $col9_2 = $col9;
	my $type = $array[2];
	my $id = "NA";
	my $name = "NA";
	my $parent;
	my %col9 = ();
	foreach my $element (@col9){
	    my @temp = split "=", $element;
	    $col9{$temp[0]} = $temp[1];
	    if ( $element =~ /ID=(.*)/ ){
		$id = $1;
	    }
	    elsif ( $element =~ /Name=(.*)/ ){
		$name = $1;
	    }
            elsif ( $element =~ /Parent=(.*)/ ){
                $parent = $1;
            }
	}
        if ( $type =~ /RNA|transcript/ ){ 
#trim whitepace from end of name
	    $name =~ s/\s+$//; 
#save last letter from name and convert to upper case (if meets the next condition, should be isoform suffix)
	    my $last_letter = uc(substr($name, -1));
#if name ends with R[A-Z], isoform [A-Z], or N[a-z]
           if ( $name =~ /R[A-Z]$|[Ii]soform\s+[A-Z]$|\s+[0-9]+[a-z]$/ ){
		my $id_letter;
		my $id_base;
		if ( $id =~ /^(.+-R)([A-Z])$/ ){
		    $id_letter = $2;
		    $id_base = $1;
		}
		if ( $id_letter =~ /$last_letter/i ){
		    $transcript_id_hash{$id} += 1;
		}
		else {
                    my $newid = $id_base.$last_letter;
		    $transcript_id_hash{$newid} += 1;
		    $mismatches{$id} = $newid;
		    $col9{'ID'} = $newid;
		    $col9_2 = paste_col9(\%col9);
		    print "modified $type ID $id to $newid based on name $name\n";
		}
	    }
	    else {
		$transcript_id_hash{$id} += 1;
	    }
	}
#child features: assumes that these are *exon, CDS, or *_utr
	elsif ( $type =~ /exon|CDS|_utr/ ){
	    if ( defined $mismatches{$parent} ){
		$col9{'Parent'} = $mismatches{$parent};
	    }
	    #make new line to print out; turn %col9 into scalar, concatenate w $line2
	    $col9_2 = paste_col9(\%col9); 
	}
	else {
#would be nice to build in a type counter to print out all types and their quantities for debugging (in case you missed some types)	    
	    $typehash{$type} += 1;
	}
	my $newline = $line2."\t".$col9_2;                                                                                                      
	print $OUT "$newline\n";                                                                                                                
    }
}

print "number of SO types in addition to *exon*, *CDS*, *_utr*, *RNA*, *transcript*\n";
foreach my $key ( keys %typehash ){
    print "$key\t$typehash{$key}\n";
}

foreach my $tkey ( keys %transcript_id_hash ){
    if ( $transcript_id_hash{$tkey} > 1 ){
	print "duplicate transcript ID $tkey generated $transcript_id_hash{$tkey} times.\n";
    }
}

sub paste_col9 {
    my %in = %{$_[0]};
    my $newcol9;
    foreach my $key ( keys %in){
	$newcol9 .= "$key=$in{$key};";
    }
    return $newcol9;
    $newcol9 = ();
}

