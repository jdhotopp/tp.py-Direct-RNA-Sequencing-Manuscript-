#!/usr/local/bin/perl -w

## Function: take transcript prediction gff and merge with the region table and a genbank gff to create a new annotation table and gff file

sub Usage {
    print <<_EOB_;
merge_transcripts_de - a program for to take a merged file with transcripts and CDSs, 6 outputs based on CDSs, 6 outputs based on transcripts, and 2 outputs based on transcripts, and merge them
Usage: merge_transcripts_de.pl <switches/options>
Options: 
-g: Input combined gff file (e.g. from merge_transcripts_gff)
-a: Input DE file 1 for CDS from FADU
-b: Input DE file 2 for CDS from FADU
-c: Input DE file 3 for CDS from FADU
-d: Input DE file 4 for CDS from FADU
-e: Input DE file 5 for CDS from FADU
-f: Input DE file 6 for CDS from FADU
-i: Input DE file 1 for transcripts from Salmon
-j: Input DE file 2 for transcripts from Salmon
-k: Input DE file 3 for transcripts from Salmon
-l: Input DE file 4 for transcripts from Salmon
-m: Input DE file 5 for transcripts from Salmon
-n: Input DE file 6 for transcripts from Salmon
-p: Input DE file 7 for transcripts from Salmon (ONT)
-q: Input DE file 8 for transcripts from Salmon (ONT)
-o: output file
_EOB_
    exit;
}

use Getopt::Std;
getopts ('hg:o:a:b:c:d:e:f:i:j:k:l:m:n:p:q:'); #colon after letter means it is an option, not a switch
Usage if $opt_h;
$gff_file=$opt_g;
Usage if ($opt_g =~ /^$/);
$cds1=$opt_a;
Usage if ($opt_a =~ /^$/);
$cds2=$opt_b;
Usage if ($opt_b =~ /^$/);
$cds3=$opt_c;
Usage if ($opt_c =~ /^$/);
$cds4=$opt_d;
Usage if ($opt_d =~ /^$/);
$cds5=$opt_e;
Usage if ($opt_e =~ /^$/);
$cds6=$opt_f;
Usage if ($opt_f =~ /^$/);
$transcript_file1=$opt_i;
Usage if ($opt_i =~ /^$/);
$transcript_file2=$opt_j;
Usage if ($opt_j =~ /^$/);
$transcript_file3=$opt_k;
Usage if ($opt_k =~ /^$/);
$transcript_file4=$opt_l;
Usage if ($opt_l =~ /^$/);
$transcript_file5=$opt_m;
Usage if ($opt_m =~ /^$/);
$transcript_file6=$opt_n;
Usage if ($opt_n =~ /^$/);
$transcript_file7=$opt_p;
Usage if ($opt_p =~ /^$/);
$transcript_file8=$opt_q;
Usage if ($opt_q =~ /^$/);
$output=$opt_o;
Usage if ($opt_o =~ /^$/);

#The gene name will match the parent for the CDS, and the gene parent is the transcript name.
#Gene to CDS is 1-to-1, but gene to transcript is not, with multiple genes in a transcript and multiple transcripts for a gene
#FADU output will have a row for each CDS and Salmon will have a row for each transcripts.  

open (OUTPUT, ">$output") || die "Cannot open $output : $!\n";

open (CDS, "$cds1") || die "Cannot open $cds1 : $!\n";
while (<CDS>) {
    chomp;
    @line=split;
    $tpm_CDS1{$line[0]}=$line[4];
}
close CDS;
open (CDS, "$cds2") || die "Cannot open $cds2 : $!\n";
while (<CDS>) {
    chomp;
    @line=split;
    $tpm_CDS2{$line[0]}=$line[4];
}
close CDS;
open (CDS, "$cds3") || die "Cannot open $cds3 : $!\n";
while (<CDS>) {
    chomp;
    @line=split;
    $tpm_CDS3{$line[0]}=$line[4];
}
close CDS;
open (CDS, "$cds4") || die "Cannot open $cds4 : $!\n";
while (<CDS>) {
    chomp;
    @line=split;
    $tpm_CDS4{$line[0]}=$line[4];
}
close CDS;
open (CDS, "$cds5") || die "Cannot open $cds5 : $!\n";
while (<CDS>) {
    chomp;
    @line=split;
    $tpm_CDS5{$line[0]}=$line[4];
}
close CDS;
open (CDS, "$cds6") || die "Cannot open $cds6 : $!\n";
while (<CDS>) {
    chomp;
    @line=split;
    $tpm_CDS6{$line[0]}=$line[4];
}
close CDS;

open (TRANSCRIPT, "$transcript_file1") || die "Cannot open $transcript_file1 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran1{$line[0]}=$line[3];
}
close TRANSCRIPT;
open (TRANSCRIPT, "$transcript_file2") || die "Cannot open $transcript_file2 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran2{$line[0]}=$line[3];
}
close TRANSCRIPT;
open (TRANSCRIPT, "$transcript_file3") || die "Cannot open $transcript_file3 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran3{$line[0]}=$line[3];
}
close TRANSCRIPT;
open (TRANSCRIPT, "$transcript_file4") || die "Cannot open $transcript_file4 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran4{$line[0]}=$line[3];
}
close TRANSCRIPT;
open (TRANSCRIPT, "$transcript_file5") || die "Cannot open $transcript_file5 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran5{$line[0]}=$line[3];
}
close TRANSCRIPT;
open (TRANSCRIPT, "$transcript_file6") || die "Cannot open $transcript_file6 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran6{$line[0]}=$line[3];
}
close TRANSCRIPT;
open (TRANSCRIPT, "$transcript_file7") || die "Cannot open $transcript_file7 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran7{$line[0]}=$line[3];
}
close TRANSCRIPT;
open (TRANSCRIPT, "$transcript_file8") || die "Cannot open $transcript_file8 : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    @line=split;
    $line[0] =~ s/ID=//;
    $tpm_tran8{$line[0]}=$line[3];
}
close TRANSCRIPT;

open (GFF, "$gff_file") || die "Cannot open $gff_file : $!\n";
while (<GFF>) {
    if (!($_ =~ /^#/)) {
	chomp;
	@line=split;
	if ($line[2]=~/CDS/) {  #CDS will have one parent
	    if ($_=~/Parent=([\w\-]*);/) {
		$gene_CDS=$1;
	    }
	    if ($_=~/ID=([\w\-\.]*);/) {
		$cds=$1;
	    }
	    $cds_hash{$gene_CDS}=$cds;
	}
    }
}
close GFF;

open (GFF, "$gff_file") || die "Cannot open $gff_file : $!\n";
while (<GFF>) {
    if ((!($_ =~ /^#/))&&(!($_=~/gene_biotype=tRNA/))&&(!($_=~/gene_biotype=rRNA/))&&(!($_=~/gene_biotype=RNase_P_RNA/))&&(!($_=~/gene_biotype=ncRNA/))&&(!($_=~/gene_biotype=tmRNA/))&&(!($_=~/gene_biotype=antisense_RNA/))&&(!($_=~/gene_biotype=SRP_RNA/))) {
	chomp;
	@line=split;
	if ($line[2]=~/gene/) {  #genes will have multiple parents and some genes will not have a parent
	    if ($_=~/Parent=ID=([\w\-]*)/) {
		$transcript=$1;
	    }
	    if ($_=~/\tID=([\w\-\_]*);/) {
		$gene_tran=$1;
	    }
	    $cds=$cds_hash{$gene_tran};
	    print OUTPUT "$transcript\t$gene_tran\t$cds\t$tpm_tran1{$transcript}\t$tpm_tran2{$transcript}\t$tpm_tran3{$transcript}\t$tpm_tran4{$transcript}\t$tpm_tran5{$transcript}\t$tpm_tran6{$transcript}\t$tpm_tran7{$transcript}\t$tpm_tran8{$transcript}\t$tpm_CDS1{$cds}\t$tpm_CDS2{$cds}\t$tpm_CDS3{$cds}\t$tpm_CDS4{$cds}\t$tpm_CDS5{$cds}\t$tpm_CDS6{$cds}\n";
#	    print OUTPUT "$transcript\t$gene_tran\t$cds\t$tpm_tran{$transcript}\t$tpm_CDS{$cds_hash{$gene_tran}}\n";
	}
    }
}
close GFF;
close OUTPUT;

__END__

    
