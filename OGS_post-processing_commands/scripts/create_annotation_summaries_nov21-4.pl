#!usr/bin/perl -w
use strict;


my $gff = shift @ARGV or die;
#need fasta file to identify whether start + stop codons are present in CDS; and to adjust aa sequence length if stop codon is present
my $fasta = shift @ARGV or die;
#gene id conversion file between dan's and Web Apollo IDs                                                                                                                    
my $gene_convert = shift @ARGV or die;
#transcript id conversion file between dan's and Web Apollo IDs                                                                                                              
my $transcript_convert = shift @ARGV or die;

my $out = shift @ARGV or die;
#species code (e.g. anogla) - needed to generate URL
my $species_code = shift @ARGV or die;

my %convert = ();
open ( my $GCONVERT, $gene_convert ) or die;
while ( my $line = <$GCONVERT> ){
    chomp $line;
    my @array =split '\t', $line;
    $convert{$array[0]} = $array[1];
}
close $GCONVERT;

open ( my $TCONVERT, $transcript_convert ) or die;
while ( my $line = <$TCONVERT> ){
    chomp $line;
    my @array =split '\t', $line;
    $convert{$array[0]} = $array[1];
}
close $TCONVERT;

open ( my $GFF, $gff ) or die;
my %gene_ids = ();
my %transcript_ids = ();
my %aas = ();
my %num_cds_introns = ();
my %CDS_start = ();
my %CDS_stop = ();
my %num_exon_introns = ();
my %cds_true_stop_coordinate = ();
my %CDS_phase = ();
my %stop_codon_readthrough = ();
my %CDS = ();
my %transcript_extent = ();
my %seqmods = ();
my %non_canonical_splice_sites = ();
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
	my @col9 = split /;/, $array[8];
        my $scaffold = $array[0];
	my $strand = $array[6];
        my $start = $array[3];
	my $stop = $array[4];
	my $id;
	my $owner;
	my $name; 
	my $symbol="NA";
	my $creation_date;
	my $comments="NA";
	my $status="NA";
	my $parent="NA";
	foreach my $element (@col9){
	    if ( $element =~ /ID=(.*)/ ){
		$id = $1;
	    }
	    elsif ( $element =~ /owner=(.*)/ ){
		$owner = $1;
	    }
	    elsif ( $element =~ /Name=(.*)/ ){
		$name = $1;
	    }
            elsif ( $element =~ /symbol=(.*)/ ){
                $symbol = $1;
            }
            elsif ( $element =~ /date_creation=(.*)/ ){
                $creation_date = $1;
            }
            elsif ( $element =~ /Note=(.*)/ ){
                $comments = $1;
            }
	    elsif ( $element =~ /status=(.*)/ ){
		$status = $1;
	    }
	    elsif ( $element =~ /Parent=(.*)/ ){
		$parent = $1;
	    }
	}
#populate gene/pseudogene  hash
	if ( $array[2] =~ /gene|pseudogene/ ){
#check for new ID in conversion file
	    if ( defined $convert{$id} ){
#in gene hash, key: id; value: name, SO type, date created, notes.	    
		$gene_ids{$id} = "$name\t$convert{$id}\t$array[2]\t$creation_date\t$comments";
	    }
	    else {
		print "$id is not defined in the conversion file; not generating output for this ID. Here is the line:\n$line\n";
		die;
	    }
	}
#populate transcript hash
	elsif ( $array[2] =~ /mRNA|transcript|pseudogenic_transcript|rRNA/ ){
	    if ( defined $gene_ids{$parent} ){
		if ( defined $convert{$id} ){
		    my $link = "https://apollo.nal.usda.gov/".$species_code."/jbrowse/?loc=".$scaffold."%3A".$start."..".$stop."&tracks=DNA%2CAnnotations%2C".$species_code."_current_models&highlight=";
		    $transcript_ids{$id} = "$gene_ids{$parent}\t$owner\t$scaffold\t$start\t$stop\t$strand\t$array[2]\t$name\t$convert{$id}\t$comments\t$link";
		    $transcript_extent{$scaffold} -> {$strand} -> {$start} -> {$id}  = $stop;
		}
		else {
		    print "$id is not defined in the conversion file; not generating output for this ID. Here is the line:\n$line\n";
		}
	    }
	    else {
		warn "parents and children out of synch here:\n$parent\t$id\n";
	    }
	}
#populate transcript_ids hash w sequence mods (top-level, no parents, no children)	
	elsif ( $array[2] =~ /deletion|insertion|substitution/ ){
	    $seqmods{$scaffold} -> {$strand} -> {$start} = $stop;
	    my $link = "https://apollo.nal.usda.gov/".$species_code."/jbrowse/?loc=".$scaffold."%3A".$start."..".$stop."&tracks=DNA%2CAnnotations%2C".$species_code."_current_models&highlight=";
	    $transcript_ids{$id} = "$name\tNA\t$array[2]\t$creation_date\tNA\t$owner\t$scaffold\t$start\t$stop\t$strand\t$array[2]\t$id\tNA\tNA\t$link";
	}
	elsif ( $array[2] =~ /CDS/ ){
	    $CDS{$id} = $parent;
#align CDS feature with transcript parent ID
	    if ( defined $transcript_ids{$parent} ){
#add aa's
		$aas{$parent} += ( abs($start - $stop) +1)/3;
		#count CDS features to derive number of introns
		$num_cds_introns{$parent} += 1;
		#define CDS start and stop coordinate (Hugh Robertson's request); include phase
		if ( defined $CDS_start{$parent} ){
		    if ( $start < $CDS_start{$parent} ){
			$CDS_start{$parent} = $start;
			$CDS_phase{$parent}{$start} = $array[7];
		    }
		    if ( $stop > $CDS_stop{$parent} ){
			$CDS_stop{$parent} = $stop;
			$CDS_phase{$parent}{$stop} = $array[7];
		    }
		}
		else {
		    $CDS_start{$parent} = $start;
		    $CDS_stop{$parent} = $stop;
		    $CDS_phase{$parent}{$start} = $array[7];
		    $CDS_phase{$parent}{$stop} = $array[7];
		}
	    }
#already pre-populate true stop coordinate hash for proper aa length calculaton (accounting for stop codons; complete value in next section. id -> scaf  -> dir -> stop                                                                         
	    $cds_true_stop_coordinate{$parent}{$scaffold}{$array[6]} = 1;
	}
	elsif ( $array[2] =~ /exon/ ){
	    #define the total number of introns (within CDS and UTRs combined)
	    $num_exon_introns{$parent} += 1;
	}
#check for stop codon readthroughs to associate w CDS                                                                           
        elsif ( $array[2] =~ /stop_codon_read_through/ ){
	    if ( defined $CDS{$parent} ){
#hash: key, parent mRNA ID; value, CDS ID
		$stop_codon_readthrough{$CDS{$parent}} = "stop_codon_readthrough";
		print "stop codon readthrough: id: $id, parent: $parent, mRNA ID: $CDS{$parent}\n";
	    }
        }
	elsif ( $array[2] =~ /non_canonical_/ ){
	    $non_canonical_splice_sites{$parent} += 1;
	}
#grab bag of other features for accounting's sake
	else {
	    print "don't know what to do with feature $id of type $array[2]; therefore not reflected in output\n";
	}
    }
}

foreach my $key ( keys %transcript_ids ){
    if ( defined $non_canonical_splice_sites{$key} ){
	$transcript_ids{$key} .= "\t$non_canonical_splice_sites{$key}";
    }
    else {
	$transcript_ids{$key} .= "\tNA";
    }
}

#will have to slice out string from last 3 bases of CDS (note that this is specific to Web Apollo coding - not all gff3s code their stop codons as CDS) and see whether it matches stop codons (or their reverse complement, depending on strand)         
#TAG, TAA, TGA                                                                                                              
close $GFF;

# #Determine overlap between sequence mod and sequence feature
# foreach my $scafkey ( keys %transcript_extent ){
#     foreach my $strandkey ( keys %{$transcript_extent{$scafkey}} ){
# 	foreach my $startkey ( keys %{$transcript_extent{$scafkey}{$strandkey}} ){
# 	    foreach my $idkey ( keys %{$transcript_extent{$scafkey}{$strandkey}{$startkey}} ){
# 		if ( defined $seqmods{$scafkey}{$strandkey} ){
# 		    print "made it\n";
# 		    foreach my $start2key ( keys %{$seqmods{$scafkey}{$strandkey}} ){
# 			if ( $start2key > $startkey && $seqmods{$scafkey}{$strandkey}{$start2key} < $transcript_extent{$scafkey}{$strandkey}{$startkey}{$idkey} ){
# 			    print "$idkey has a sequence mod\n";
# 			    $transcript_ids{$idkey} .= "has-sequence-mod";
# 			}
# 			else {
# 			    print "$idkey didn't make it; seq mod start and stop: $start2key, $seqmods{$scafkey}{$strandkey}{$start2key}, mRNA start and stop: $startkey, $transcript_extent{$scafkey}{$strandkey}{$startkey}{$idkey}, scaffold $scafkey\n";
# 			}
# 		    }
# 		}
# 	    }
# 	}
#     }
# }
# die;


open ( my $FASTA, $fasta ) or die;
my %fasta = ();
my $defline;
while ( my $fline = <$FASTA> ){
    chomp $fline;
    if ( $fline =~ /^>(.+)/ ){
	$defline = $1;
    }
    else {
	$fasta{$defline} .= $fline;
    }
}

#id -> scaf -> dir -> 1 (should be populated with true stop coordinate  
#code to determine whether stop codon is present in CDS and calculate true aa sequence length
#WIP: STILL NO CODE TO INCORPORATE PHASE; will probably have to run gff through a translator
foreach my $cds ( keys %cds_true_stop_coordinate ){
    foreach my $scafkey ( keys %{$cds_true_stop_coordinate{$cds}} ){
	if ( defined $fasta{$scafkey} ){
	    foreach my $dirkey ( keys %{$cds_true_stop_coordinate{$cds}{$scafkey}} ){
		if ( $dirkey eq "+" ){		    
		    #for forward strand, stop coordinate in CDS_stop is true stop coordinate, and coordinate in CDS_start is true start coordinate
		    my $forward_stop = $CDS_stop{$cds};
		    my $stop_slice = substr( $fasta{$scafkey}, ($forward_stop -3), 3);
		    if ( $stop_slice =~ /^TAG$|^TAA$|^TGA$/ ){
			my $newlength1 =$aas{$cds} -1;
			$aas{$cds} = "$newlength1\ty";
		    }
		    else {
			$aas{$cds} = "$aas{$cds}\tn";
			print "transcript $cds (forward strand) does not have a stop codon: $stop_slice\n";
		    }
                    my $forward_start =$CDS_start{$cds};
		    my $start_slice = substr($fasta{$scafkey}, ($forward_start -1), 3);
		    if ( $start_slice =~ /^ATG$/ ){
#if ATG is present, add info that start codon is present in CDS to amino acid length hash
			$aas{$cds} .= "\ty";
		    }
		    else {
			$aas{$cds} .= "\tn";
			print "transcript $cds (forward strand) does not have a start codon: $start_slice\n";
		    }
		}
		elsif ( $dirkey eq "-" ){
		    #for reverse strand,  coordinate in CDS_start is true stop coordinate
		    my $reverse_stop =$CDS_start{$cds};
		    my $reverse_stop_slice = substr( $fasta{$scafkey}, ($reverse_stop -1 ), 3);
		    if ( $reverse_stop_slice =~ /^CTA$|^TTA$|^TCA$/ ){
			my $newlength = $aas{$cds} -1;
			$aas{$cds} = "$newlength\ty";
		    }
		    else {
			$aas{$cds} = "$aas{$cds}\tn";
			print "transcript $cds (reverse strand) does not have a stop codon: $reverse_stop_slice\n";
		    }
		    my $reverse_start = $CDS_stop{$cds};
		    my $reverse_start_slice = substr( $fasta{$scafkey}, ($reverse_start -3 ), 3);
		    if ( $reverse_start_slice =~ /^CAT$/ ){
			#if ATG is present, add info that start codon is present in CDS to amino acid length hash
                        $aas{$cds} .= "\ty";
		    }
		    else {
			$aas{$cds} .= "\tn";
			print "transcript $cds (reverse strand) does not have a start codon: $reverse_start_slice\n";
		    }
		}
	    }
#	    }
	}
    }
}

open ( my $OUT, ">$out" ) or die;
#print out header
print $OUT "gene_name\tgene_ID\tgene_type\tgene_creation_date\tgene_comments\ttranscript_owner\ttranscript_scaffold\ttranscript_start\ttranscript_end\ttranscript_strand\ttranscript_type\ttranscript_name\ttranscript_ID\ttranscript_comments\ttranscript_URL\t#non-canonical-splice-sites\t#amino_acids\tstop-codon-present?\tstart-codon-present?\tCDS_start_coordinate\tCDS_end_coordinate\t#introns_in_CDS\ttotal#introns\thas_stop_codon_read_through\n";
foreach my $key (keys %transcript_ids ){
#add information from CDS and exon hashes
##NEED CONDITION FOR AASAND CDS START STOP IF TRANSCRIPT TYPE IS NOT MRNA - think I covered this with a band-aid by  subbing in 'NA'. Should probably build in a better check for failing to meet conditions
#if it has exons
    if ( defined $num_exon_introns{$key} ){
#if it has CDS
	if ( defined $aas{$key} and defined $CDS_start{$key} and defined $CDS_stop{$key} ){
#if it has a readthrough stop codon	    
	    if ( defined $stop_codon_readthrough{$key} ){
		print $OUT "$transcript_ids{$key}\t$aas{$key}\t$CDS_start{$key}\t$CDS_stop{$key}\t$num_cds_introns{$key}\t$num_exon_introns{$key}\t$stop_codon_readthrough{$key}\n";
	    }
	    else {
                print $OUT "$transcript_ids{$key}\t$aas{$key}\t$CDS_start{$key}\t$CDS_stop{$key}\t$num_cds_introns{$key}\t$num_exon_introns{$key}\tNA\n";
	    }
	}
#if it doesn't have CDS
	else {
	    print $OUT "$transcript_ids{$key}\tNA\tNA\tNA\tNA\tNA\tNA\t$num_exon_introns{$key}\tNA\n";
	}
    }
#if it doesn't have exons    
    else {
	print $OUT "$transcript_ids{$key}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
    }
}
