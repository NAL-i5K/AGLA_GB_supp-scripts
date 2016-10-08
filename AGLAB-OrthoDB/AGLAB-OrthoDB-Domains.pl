# Author: Robert Waterhouse
#!/usr/local/bin/perl

use strict;
use DBI;

my $host = "localhost";
my $database = "orthodb8";  
my $user = "waterhouse";
my $pw = "password"; 
my $dbh = DBI->connect("DBI:mysql:database=$database;host=$host",$user,$pw);

my $node='Arthropoda';
my @species=qw(ZNEVA PHUMA APISU AMELL NVITR AGLAB DPOND TCAST APLAN OTAUR PXYLO DPLEX MDEST DMELA AGAMB);
my $allspecs=join("','",@species);

undef my %int2pub;
undef my %pub2int;
undef my %int2ipr;
undef my %tax2spe;
undef my %int2sp;

foreach my $spec (@species) {
    my @lines=split(/\n/, `zcat /OrthoDB8/i5K/IPRS/$spec.idmap.txt.gz`);
    foreach my $line (@lines) {
	chomp($line);
	my($int,$pub)=split(/\s+/,$line);
	$int2pub{$int}=$pub;
	$pub2int{"$spec\:$pub"}=$int;
	my($aa,$bb)=split(/:/,$int);
	$tax2spe{$aa}=$spec;
	$int2sp{$int}=$spec;
    }
    @lines=split(/\n/, `zcat /OrthoDB8/i5K/IPRS/$spec\_*.iprs.txt.gz`);
    foreach my $line (@lines) {
	chomp($line);
	my($int,$ipr)=split(/\t/,$line);
	$int2ipr{$int}=$ipr;
    }
}

open(IN,"entry.list") || die $!;
my @lines=<IN>;
close(IN);
undef my %iprs;
foreach my $line (@lines) {
    #if($line=~/Glycoside hydrolase/ || $line=~/glycoside hydrolase/) {
    if($line=~/Peptidase/ || $line=~/peptidase/ || $line=~/Proteinase/ || $line=~/proteinase/ || $line=~/Protease/ || $line=~/protease/) {
    #if($line=~/P450/ || $line=~/cytochrome P450/ || $line=~/Cytochrome P450/) {
    #if($line=~/Esterase/ || $line=~/esterase/) {
    #if($line=~/glucosyltransferase/ || $line=~/Glucosyltransferase/ || $line=~/glucuronosyl/) {

	# EXCLUDE
	if($line=~/^IPR000801/) { next; }     # Putative esterase
	if($line=~/^IPR001087/) { next; }     # GDSL lipase/esterase
	if($line=~/^IPR015955/) { next; }     # Lactate dehydrogenase/glycoside hydrolase, family 4, C-terminal
	if($line=~/^IPR022700/) { next; }     # Proteinase, regulatory CLIP domain
	if($line=~/^IPR009020/) { next; }     # Proteinase inhibitor, propeptide
	if($line=~/^IPR003146/) { next; }     # Proteinase inhibitor, carboxypeptidase propeptide

	chomp($line);
	if($line=~/^(IPR\d+)\s+(.+)/) { $iprs{$1}=$2; }
    }
}

# redundant iprs not for chart
undef my %redundants;

# glycoside hydrolases
$redundants{'IPR017853'}="Glycoside hydrolase superfamily";
$redundants{'IPR011330'}="Glycoside hydrolase/deacetylase, beta/alpha-barrel";
$redundants{'IPR015341'}="Glycoside hydrolase, family 38, central domain";
$redundants{'IPR028995'}="Glycoside hydrolase, families 57/38, central domain";
$redundants{'IPR025887'}="Glycoside hydrolase family 31, N-terminal domain";
$redundants{'IPR000111'}="Glycoside hydrolase family 27/36, conserved site";

# cyp450s
$redundants{'IPR001128'}="Cytochrome P450";

# esterases
$redundants{'IPR023088'}="3'5'-cyclic nucleotide phosphodiesterase";
$redundants{'IPR006683'}="Thioesterase superfamily";

# proteases
$redundants{'IPR001254'}="Serine proteases, trypsin domain";
$redundants{'IPR001314'}="Peptidase S1A, chymotrypsin-type";
$redundants{'IPR008737'}="Peptidase aspartic, putative";
$redundants{'IPR021109'}="Aspartic peptidase domain";
$redundants{'IPR018061'}="Peptidase A2A, retrovirus RVP subgroup";
$redundants{'IPR018497'}="Peptidase M13, C-terminal domain";
$redundants{'IPR006026'}="Peptidase, metallopeptidase";
$redundants{'IPR011650'}="Peptidase M20, dimerisation domain";
$redundants{'IPR015917'}="Peptidase C14A, caspase precursor p45, core";
$redundants{'IPR011249'}="Metalloenzyme, LuxS/M16 peptidase-like";
$redundants{'IPR000101'}="Gamma-glutamyltranspeptidase";
$redundants{'IPR002870'}="Peptidase M12B, propeptide";
$redundants{'IPR001461'}="Aspartic peptidase";
$redundants{'IPR002469'}="Dipeptidylpeptidase IV, N-terminal domain";
$redundants{'IPR011765'}="Peptidase M16, N-terminal";
$redundants{'IPR028980'}="Creatinase/Aminopeptidase P, N-terminal";
$redundants{'IPR012338'}="Beta-lactamase/transpeptidase-like";
$redundants{'IPR011356'}="Leucine aminopeptidase/peptidase B";
$redundants{'IPR022684'}="Peptidase C2, calpain family";
$redundants{'IPR022683'}="Peptidase C2, calpain, domain III";
$redundants{'IPR022682'}="Peptidase C2, calpain, large subunit, domain III";
$redundants{'IPR008283'}="Peptidase M17, leucyl aminopeptidase, N-terminal";
$redundants{'IPR000243'}="Peptidase T1A, proteasome beta-subunit";
$redundants{'IPR003111'}="ATP-dependent protease La (LON), substrate-binding domain";
$redundants{'IPR008969'}="Carboxypeptidase-like, regulatory domain";
$redundants{'IPR013273'}="Peptidase M12B, ADAM-TS";
$redundants{'IPR019759'}="Peptidase S24/S26A/S26B";
$redundants{'IPR000246'}="Peptidase T2, asparaginase 2";
$redundants{'IPR023302'}="Peptidase S9A, N-terminal domain";
$redundants{'IPR007865'}="Aminopeptidase P, N-terminal";
$redundants{'IPR015500'}="Peptidase S8, subtilisin-related";
$redundants{'IPR012599'}="Peptidase C1A, propeptide";


undef my %ipr2spe2cnt;
undef my %ipr2tot2cnt;
undef my %ipr2genes;
foreach my $int (keys %int2ipr) {
    foreach my $ipr (keys %iprs) {
	if($int2ipr{$int}=~/$ipr/) {
	    my($aa,$bb)=split(/:/,$int);
	    my $spe=$tax2spe{$aa};
	    if(defined($ipr2spe2cnt{$ipr}{$spe})) { $ipr2spe2cnt{$ipr}{$spe}++; }
	    else { $ipr2spe2cnt{$ipr}{$spe}=1; }
	    if(defined($ipr2tot2cnt{$ipr})) { $ipr2tot2cnt{$ipr}++; }
	    else { $ipr2tot2cnt{$ipr}=1; }
	    push(@{$ipr2genes{$ipr}},$int);
	}
    }
}

undef my %ipr2max;
foreach my $ipr (keys %ipr2tot2cnt) {
    foreach my $spe (@species) { 
	if(defined($ipr2max{$ipr})) {
	    if(defined($ipr2spe2cnt{$ipr}{$spe})) {
		if($ipr2spe2cnt{$ipr}{$spe}>$ipr2max{$ipr}) { $ipr2max{$ipr}=$ipr2spe2cnt{$ipr}{$spe}; }
	    }
	}
	else { $ipr2max{$ipr}=$ipr2spe2cnt{$ipr}{$spe}; }
    }
}


open(OUT1, ">iprs.txt") || die $!;
open(OUT2, ">iprs4chart.txt") || die $!;
print OUT1 "IPR\tIPR_name";
print OUT2 "IPR\tIPR_name";
foreach my $spe (@species) { print OUT1 "\t$spe"; print OUT2 "\t$spe";}
print OUT1 "\tTotal\tMax\n";
print OUT2 "\tTotal\tMax\n";
undef my %genes;
foreach my $ipr (sort { $ipr2max{$b} <=> $ipr2max{$a} } keys %ipr2max) {

    if($ipr2max{$ipr}<=5) { next; }
    print OUT1 "$ipr\t$iprs{$ipr}";
    if(!defined($redundants{$ipr})) { print OUT2 "$ipr\t$iprs{$ipr}"; }

    foreach my $spe (@species) { 
	if(defined($ipr2spe2cnt{$ipr}{$spe})) { 
	    print OUT1 "\t$ipr2spe2cnt{$ipr}{$spe}"; 
	    if(!defined($redundants{$ipr})) { print OUT2 "\t$ipr2spe2cnt{$ipr}{$spe}"; }
	}
	else { 
	    print OUT1 "\t0"; 
	    if(!defined($redundants{$ipr})) { print OUT2 "\t0"; }
	}
    }

    print OUT1 "\t$ipr2tot2cnt{$ipr}\t$ipr2max{$ipr}\n";
    if(!defined($redundants{$ipr})) { 
	print OUT2 "\t$ipr2tot2cnt{$ipr}\t$ipr2max{$ipr}\n";
	print OUT2 "$ipr\tDIFF\-$iprs{$ipr}";
	foreach my $spe (@species) {
	    print OUT2 "\t";
	    if(defined($ipr2spe2cnt{$ipr}{$spe})) { print OUT2 $ipr2max{$ipr}-$ipr2spe2cnt{$ipr}{$spe}; }
	    else { print OUT2 $ipr2max{$ipr}; }
	}
	print OUT2 "\t$ipr2tot2cnt{$ipr}\t$ipr2max{$ipr}\n";
	foreach my $gen (@{$ipr2genes{$ipr}}) { $genes{$gen}=1; }
    }
    
}
close(OUT1);
close(OUT2);





my $allgenes=join("','", sort keys %genes);
my $sth = $dbh->prepare("select intid, ogid, spec from i5k_groups where node='$node' and intid in ('$allgenes')");
$sth->execute();
undef my %ogids;
undef my %int2og;
while (my $ref = $sth->fetchrow_hashref()) {
    $ogids{$ref->{'ogid'}}=1;
    $int2og{$ref->{'intid'}}=$ref->{'ogid'};
}
$sth->finish();

my $allogs=join("','", sort keys %ogids);
$sth = $dbh->prepare("select * from i5k_groups where node='$node' and ogid in ('$allogs') and spec in ('$allspecs')");
$sth->execute();
undef my %og2sp2gn;
undef my %gn2og;
while (my $ref = $sth->fetchrow_hashref()) {
    $gn2og{$ref->{'pubid'}}=$ref->{'ogid'};
    if(defined($og2sp2gn{$ref->{'ogid'}}{$ref->{'spec'}})) { $og2sp2gn{$ref->{'ogid'}}{$ref->{'spec'}}++; }
    else { $og2sp2gn{$ref->{'ogid'}}{$ref->{'spec'}}=1; }
}
$sth->finish();

$dbh->disconnect();

undef my %spec2sc;
undef my %spec2mc;
undef my %spec2ot;
foreach my $gn (keys %genes) {
    if(defined($int2og{$gn})) {
	my $og=$int2og{$gn};
	my $sp=$int2sp{$gn};
	if($og2sp2gn{$og}{$sp}==1) { 
	    if(defined($spec2sc{$sp})) { $spec2sc{$sp}++ }
	    else { $spec2sc{$sp}=1; }
	}
	else {
	    if(defined($spec2mc{$sp})) { $spec2mc{$sp}++ }
	    else { $spec2mc{$sp}=1; }
	}
    }
    else {
	if(defined($spec2ot{$int2sp{$gn}})) { $spec2ot{$int2sp{$gn}}++; }
	else { $spec2ot{$int2sp{$gn}}=1; }
    }
}

open(OUT,">iprs-ortho.txt") || die $!;

print OUT "Spec\tSC\tMC\tOT\n";
foreach my $sp (@species) { 
    print OUT "$sp";
    print OUT "\t$spec2sc{$sp}";
    print OUT "\t$spec2mc{$sp}";
    if(defined($spec2ot{$sp})) { print OUT "\t$spec2ot{$sp}"; }
    else { print OUT "\t0"; }
    print OUT "\n";
}
close(OUT);

