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
my @set=qw(ZNEVA+PHUMA+APISU+AMELL+NVITR AGLAB+DPOND+TCAST+APLAN+OTAUR PXYLO+DPLEX+MDEST+DMELA+AGAMB);
my @bitMask=();
foreach my $spec (@set) { push(@bitMask,0); }

my %totals = (
    'PHUMA' => 10773,
    'ZNEVA' => 15860,
    'APISU' => 36195,
    'DPLEX' => 16254,
    'PXYLO' => 18073,
    'MDEST' => 20163,
    'AMELL' => 15314,
    'NVITR' => 24389,
    'DMELA' => 13954,
    'AGAMB' => 12810,
    'TCAST' => 16524,
    'DPOND' => 13088,
    'APLAN' => 15497,
    'OTAUR' => 17483,
    'AGLAB' => 22035
    );

my %totals2=%totals;

my @subsets=();
getSubset(\@bitMask, \@set) while ( genMask(\@bitMask, \@set) );

undef my %ss2inc;
undef my %ss2exc;
undef my %ss2sor;

foreach my $sset (@subsets) {
    $sset=~s/^:://;
    my @incs=split(/::/,$sset);
    @{$ss2inc{$sset}}=@incs;
    undef my %in;
    foreach my $inc (@incs) { $in{$inc}=1; }
    foreach my $spe (@set) {
	if(!defined($in{$spe})) { push(@{$ss2exc{$sset}},$spe); }
    }
    if(scalar(@incs)==scalar(@set)) {next;}
    push(@{$ss2sor{$sset}},scalar(@incs),$sset);
}


# ALL
open(OUT,">beetle.venn4.txt") || die $!;
print OUT "library(gplots)\n";
my $gcnt=0;
foreach my $inc (@set) {
    $gcnt++;
    my $pres=$inc;
    $pres=~s/\+/\>0\) \+ \(/g;
    my $sth = $dbh->prepare("select ogid from i5k_counts where node='$node' and ( ($pres>0) >1)");
    $sth->execute();
    undef my %ogids;
    while (my $ref = $sth->fetchrow_hashref()) {
	$ogids{$ref->{'ogid'}}=1;
    }
    print OUT "#select ogid from i5k_counts where node='$node' and ( ($pres>0) >1)\n";
    print OUT "grp$gcnt<-c('". join("','",sort keys %ogids) ."')\n";
}
print OUT "input <-list(grp1,grp2,grp3)\n";
print OUT "pdf(\"beetle.venn3.pdf\")\n";
print OUT "venn(input)\n";
print OUT "dev.off()\n";
close(OUT);

my @togets=();
foreach my $inc (@set) {
    if($inc=~/\+/) {
	my @toge=split(/\+/,$inc);
	my $bits="( (" . join('>0)+(',@toge) . ">0) >1 )";
	push(@togets,$bits);
    }
    else { push(@togets,"$inc>0"); }
}

my $allspec=join(', ',@set);
$allspec=~s/\+/, /g;
@set=split(', ',$allspec);
my $present=join('>0 and ', @set);

print "$node\tsetA\tsetB\tOGs=na\tGenes=";
foreach my $spe (@set) { print "\t$spe"; }
print "\n";

undef my %spe2cnt;
undef my %spe2cnt1;
undef my %spe2cnt2;
foreach my $spe (@set) { $spe2cnt{$spe}=0; }
foreach my $spe (@set) { $spe2cnt1{$spe}=0; }
foreach my $spe (@set) { $spe2cnt2{$spe}=0; }

my $sth = $dbh->prepare("select ogid, $allspec from i5k_counts where node='$node' and $present>0");
print "$node\t$present>0\tna\t";
$sth->execute();

undef my %og2ss;
while (my $ref = $sth->fetchrow_hashref()) {
    $og2ss{$ref->{'ogid'}}="@set";
    foreach my $spe (@set) { 
	$spe2cnt{$spe}+=$ref->{$spe}; 
	if($ref->{$spe}==1) { $spe2cnt1{$spe}+=$ref->{$spe}; }
	else  { $spe2cnt2{$spe}+=$ref->{$spe}; }
    }
}
$sth->finish();

print "OGs=" . scalar(keys %og2ss) . "\tGenes=";
foreach my $spe (@set) { print "\t$spe2cnt{$spe}"; }
print "\n";

print "$node\t$present>0\tna\t";
print "OGs=" . scalar(keys %og2ss) . "\tGenes1=";
foreach my $spe (@set) { print "\t$spe2cnt1{$spe}"; }
print "\n";

print "$node\t$present>0\tna\t";
print "OGs=" . scalar(keys %og2ss) . "\tGenes2=";
foreach my $spe (@set) { print "\t$spe2cnt2{$spe}"; }
print "\n";

# now in sets
$present=join(' and ', @togets);
foreach my $spe (@set) { $spe2cnt{$spe}=0; }
foreach my $spe (@set) { $spe2cnt1{$spe}=0; }
foreach my $spe (@set) { $spe2cnt2{$spe}=0; }

$sth = $dbh->prepare("select ogid, $allspec from i5k_counts where node='$node' and $present");
print "$node\t$present\tna\t";
$sth->execute();

undef %og2ss;
while (my $ref = $sth->fetchrow_hashref()) {
    $og2ss{$ref->{'ogid'}}="@set";
    foreach my $spe (@set) { 
	$spe2cnt{$spe}+=$ref->{$spe}; 
	if($ref->{$spe}==1) { $spe2cnt1{$spe}+=$ref->{$spe}; }
	else  { $spe2cnt2{$spe}+=$ref->{$spe}; }
    }
}
$sth->finish();

print "OGs=" . scalar(keys %og2ss) . "\tGenes=";
foreach my $spe (@set) { print "\t$spe2cnt{$spe}"; $totals{$spe}=$totals{$spe}-$spe2cnt{$spe}; }
print "\n";

print "$node\t$present\tna\t";
print "OGs=" . scalar(keys %og2ss) . "\tGenes1=";
foreach my $spe (@set) { print "\t$spe2cnt1{$spe}"; }
print "\n";

print "$node\t$present\tna\t";
print "OGs=" . scalar(keys %og2ss) . "\tGenes2=";
foreach my $spe (@set) { print "\t$spe2cnt2{$spe}"; }
print "\n";


# SUBSETS
foreach my $sset ( sort { $ss2sor{$b}[0]<=>$ss2sor{$a}[0] or $ss2sor{$a}[1] cmp $ss2sor{$b}[1] } keys %ss2sor) {
    my $allspec=join(', ',@set);
    my @togets=();
    foreach my $inc (@{$ss2inc{$sset}}) {
	if($inc=~/\+/) {
	    my @toge=split(/\+/,$inc);
	    my $bits="( (" . join('>0)+(',@toge) . ">0) >1 )";
	    push(@togets,$bits);
	}
	else { push(@togets,"$inc>0"); }
    }
    my $present=join(' and ', @togets);
    
    @togets=();
    foreach my $exc (@{$ss2exc{$sset}}) {
	if($exc=~/\+/) {
	    my @toge=split(/\+/,$exc);
	    my $bits="( (" . join('>0)+(',@toge) . ">0) =0 )";
	    push(@togets,$bits);
	}
	else { push(@togets,"$exc=0"); }
    }
    my $absents=join(' and ', @togets);
    
    foreach my $spe (@set) { $spe2cnt{$spe}=0; }
    foreach my $spe (@set) { $spe2cnt1{$spe}=0; }
    foreach my $spe (@set) { $spe2cnt2{$spe}=0; }

    my $sth = $dbh->prepare("select ogid, $allspec from i5k_counts where node='$node' and $present and $absents");
    print "$node\t$present\t$absents\t";
    $sth->execute();
    
    undef my %og2ss;
    while (my $ref = $sth->fetchrow_hashref()) {
	$og2ss{$ref->{'ogid'}}=$sset;
	foreach my $spe (@set) { 
	    $spe2cnt{$spe}+=$ref->{$spe}; 
	    if($ref->{$spe}==1) { $spe2cnt1{$spe}+=$ref->{$spe}; }
	    else  { $spe2cnt2{$spe}+=$ref->{$spe}; }
	}
    }
    $sth->finish();

    print "OGs=" . scalar(keys %og2ss) . "\tGenes=";
    foreach my $spe (@set) { print "\t$spe2cnt{$spe}"; $totals{$spe}=$totals{$spe}-$spe2cnt{$spe}; }
    print "\n";

    print "$node\t$present\t$absents\t";
    print "OGs=" . scalar(keys %og2ss) . "\tGenes1=";
    foreach my $spe (@set) { print "\t$spe2cnt1{$spe}"; }
    print "\n";
    
    print "$node\t$present\t$absents\t";
    print "OGs=" . scalar(keys %og2ss) . "\tGenes2=";
    foreach my $spe (@set) { print "\t$spe2cnt2{$spe}"; }
    print "\n";
}


# ANY ORTHOLOGY
$sth = $dbh->prepare("select ogid, $allspec from i5k_counts where node='$node'");
print "$node\tANY\tANY\t";
$sth->execute();

undef %og2ss;
foreach my $spe (@set) { $spe2cnt{$spe}=0; }
foreach my $spe (@set) { $spe2cnt1{$spe}=0; }
foreach my $spe (@set) { $spe2cnt2{$spe}=0; }

while (my $ref = $sth->fetchrow_hashref()) {
    $og2ss{$ref->{'ogid'}}='any';
    foreach my $spe (@set) { 
	$spe2cnt{$spe}+=$ref->{$spe}; 
	if($ref->{$spe}==1) { $spe2cnt1{$spe}+=$ref->{$spe}; }
	else  { $spe2cnt2{$spe}+=$ref->{$spe}; }
    }
}
$sth->finish();

# ANY ORTHOLOGY
print "OGs=" . scalar(keys %og2ss) . "\tGenes=";
foreach my $spe (@set) { print "\t$spe2cnt{$spe}"; }
print "\n";

print "$node\tANY\tANY\t";
print "OGs=" . scalar(keys %og2ss) . "\tGenes1=";
foreach my $spe (@set) { print "\t$spe2cnt1{$spe}"; }
print "\n";

print "$node\tANY\tANY\t";
print "OGs=" . scalar(keys %og2ss) . "\tGenes2=";
foreach my $spe (@set) { print "\t$spe2cnt2{$spe}"; }
print "\n";


# totals
print "$node\ttotal\ttotal\tOGs=na\tGenes=";
foreach my $spe (@set) { print "\t$totals2{$spe}"; }
print "\n";



# gene list orthologues
my $allspecs=join("','",@set);
$sth = $dbh->prepare("select intid, spec from i5k_groups where node='$node' and spec in ('$allspecs')");
$sth->execute();

undef my %orth2spec;
undef my %selected;
while (my $ref = $sth->fetchrow_hashref()) {
    $orth2spec{$ref->{'intid'}}=$ref->{'spec'};
    my @bits=split(/:/,$ref->{'intid'});
    $selected{$bits[0]}=1;
}
$sth->finish();


# get all specs
$sth = $dbh->prepare("select run_long_id from ortho8_runs where url_id='$node'");
$sth->execute();

undef my %allarth;
while (my $ref = $sth->fetchrow_hashref()) {
    my @bits=split(/\_/,$ref->{'run_long_id'});
    foreach my $bit (sort @bits) { $allarth{$bit}=1; }
}
$sth->finish();


# get taxid 2 spec
undef my %tax2spe;
$sth = $dbh->prepare("select tax_id, code from ortho8_org_names");
$sth->execute();
while (my $ref = $sth->fetchrow_hashref()) {
    $tax2spe{$ref->{'tax_id'}}=$ref->{'code'};
}
$sth->finish();

$dbh->disconnect();


$allarth{'APLAN'}=1;
$allarth{'OTAUR'}=1;

undef my %homologs;
undef my %selfhomologs;
foreach my $sp1 (keys %allarth) {
    foreach my $sp2 (keys %selected) {
	if(-f "/orthodb8/metazoa/PWC/$sp1/$sp1\_$sp2\.align") {
	    open(IN,"/orthodb8/metazoa/PWC/$sp1/$sp1\_$sp2\.align") || die $!;
	}
	elsif(-f "/orthodb8/metazoa/PWC/$sp2/$sp2\_$sp1\.align") {
	    open(IN,"/orthodb8/metazoa/PWC/$sp2/$sp2\_$sp1\.align") || die $!;
	}
	elsif(-f "/orthodb8/metazoa/mapping/PWC/$sp1/$sp1\_$sp2\.align") {
	    open(IN,"/orthodb8/metazoa/mapping/PWC/$sp1/$sp1\_$sp2\.align") || die $!;
	}
	elsif(-f "/orthodb8/metazoa/mapping/PWC/$sp2/$sp2\_$sp1\.align") {
	    open(IN,"/orthodb8/metazoa/mapping/PWC/$sp2/$sp2\_$sp1\.align") || die $!;
	}
	else { print "Cannot find: $sp1 <=> $sp2\n"; next; }
	my @lines=<IN>;
	close(IN);
	foreach my $line (@lines) {
	    if($line=~/^#/) { next; }
	    if($line=~/^(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+/) {
		my $g1=$1;
		my $g2=$2;
		my $ev=$3;
		if($sp1 ne $sp2) {
		    if($ev<1e-5) { $homologs{$g1}=1; $homologs{$g2}=1; }
		}
		else {
		    if($ev<1e-5 && $g1 ne $g2) { $selfhomologs{$g1}=1; $selfhomologs{$g2}=1; }
		}
	    }
	}
    }
}

undef my %SELFHOM;
foreach my $sh (keys %selfhomologs) {
    if(!defined($homologs{$sh}) && !defined($orth2spec{$sh})) { 
	my @bits=split(/:/,$sh);
	my $spec=$bits[0];
	if(defined($tax2spe{$bits[0]})) { $spec=$tax2spe{$bits[0]}; }
	push(@{$SELFHOM{$spec}},$sh); 
    }
}
undef my %OTHEHOM;
foreach my $sh (keys %homologs) {
    if(!defined($orth2spec{$sh})) { 
	my @bits=split(/:/,$sh);
	my $spec=$bits[0];
	if(defined($tax2spe{$bits[0]})) { $spec=$tax2spe{$bits[0]}; }
	push(@{$OTHEHOM{$spec}},$sh); 
    }
}


print "$node\thomology\tother-5\tOGs=na\tGenes=";
foreach my $spe (@set) { print "\t". scalar(@{$OTHEHOM{$spe}}); }
print "\n";

print "$node\thomology\tself-5\tOGs=na\tGenes=";
foreach my $spe (@set) { print "\t". scalar(@{$SELFHOM{$spe}}); }
print "\n";



sub getSubset {
  my ($bitMask, $set) = @_;
  my $subset='';

  for (0 .. @$bitMask-1) {
    $subset.="::$set->[$_]" if $bitMask->[$_] == 1;
  }
  push(@subsets,$subset);
}

sub genMask {
  my ($bitMask, $set) = @_;
  my $i;
  for ($i = 0; $i < @$set && $bitMask->[$i]; $i++) {
    $bitMask->[$i] = 0;
  }
  if ($i < @$set) {
    $bitMask->[$i] = 1;
    return 1;
  }
  return 0;
}
