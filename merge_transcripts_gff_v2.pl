#!/usr/local/bin/perl -w

## Function: take transcript prediction gff and merge with the region table and a genbank gff to create a new annotation table and gff file

sub Usage {
    print <<_EOB_;
merge_transcripts_gff - a program for to merge a gff file, transcript region predictions, and transcript predictions adding in UTRs
Usage: merge_transcripts_gff.pl <switches/options>
Options: 
-g: Input reference gff file
-r: Input regions file
-t: Transcript gff file
-o: output file
_EOB_
    exit;
}

use Getopt::Std;
getopts ('hg:o:r:t:'); #colon after letter means it is an option, not a switch
Usage if $opt_h;
$gff_file=$opt_g;
Usage if ($opt_g =~ /^$/);
$regions=$opt_r;
Usage if ($opt_r =~ /^$/);
$transcript_file=$opt_t;
Usage if ($opt_t =~ /^$/);
$output=$opt_o;
Usage if ($opt_o =~ /^$/);

open (GFF, "$gff_file") || die "Cannot open $gff_file : $!\n";
while (<GFF>) {
    chomp;
    if ($_=~/^\#/) {
	push (@combined_gff,$_);
    } elsif ($_ =~ /\w/) {
	($ref,$alg,$type,$left_coord,$right_coord,$score,$strand,$something,$annotation)=split(/\t/);
	$alg_hash{$type}{$left_coord}{$right_coord}{$ref}=$alg;
	$score_hash{$type}{$left_coord}{$right_coord}{$ref}=$score;
	$strand_hash{$type}{$left_coord}{$right_coord}{$ref}=$strand;
	$something_hash{$type}{$left_coord}{$right_coord}{$ref}=$something;
	$annotation{$type}{$left_coord}{$right_coord}{$ref}=$annotation;
	@array = ($left_coord, $right_coord, $type, 0, $strand, $ref);
	push @coords, [@array];
    }
}
close (GFF) || die "Cannot close $gff_file : $!\n";

$temp="## GFF file from Genbank merged with transcripts file from rdoperon.py with merge_transcripts_gff.pl";
push (@combined_gff,$temp);

open (TRANSCRIPT, "$transcript_file") || die "Cannot open $transcript_file : $!\n";
while (<TRANSCRIPT>) {
    chomp;
    if ($_=~/^\#/) {
	push (@combined_gff,$_);
    } else {
	($ref,$alg,$type,$left_coord,$right_coord,$score,$strand,$something,$annotation)=split(/\t/);
	$talg_hash{$type}{$left_coord}{$right_coord}{$ref}=$alg;
	$tscore_hash{$type}{$left_coord}{$right_coord}{$ref}=$score;
	$tstrand_hash{$type}{$left_coord}{$right_coord}{$ref}=$strand;
	$tsomething_hash{$type}{$left_coord}{$right_coord}{$ref}=$something;
	$tannotation{$type}{$left_coord}{$right_coord}{$ref}=$annotation;
	@array = ($left_coord, $right_coord, $type, 0, $strand, $ref);
	push @coords, [@array];
    }
}
close (TRANSCRIPT) || die "Cannot close $transcript_file : $!\n";

#@coords_sorted = sort {$a->[0] <=> $b->[0]} @coords;
@coords_sorted = sort {$a->[5] cmp $b->[5] || $a->[0] <=> $b->[0]} @coords;

open (OUTPUT, ">$output") || die "Cannot create $output : $!\n";
print OUTPUT join ("\n", @combined_gff);
print OUTPUT "\n";
open (REGIONS, "$regions") || die "Cannot open $regions : $!\n";
while (<REGIONS>) {
    chomp;
    ($ref,$reg_left_coord,$reg_right_coord,$strand,$size)=split(/\t/);
    $j++;
    for ($m = 0; $m <= $#coords_sorted; $m++) {
	$gff1_left=$coords_sorted[$m][0];
	$gff1_right=$coords_sorted[$m][1];
        $gff1_type=$coords_sorted[$m][2];
	$gff1_flag=$coords_sorted[$m][3];
	$gff1_strand=$coords_sorted[$m][4];
	$gff1_ref=$coords_sorted[$m][5];
	if (($gff1_flag==0)&&($gff1_strand eq $strand)&&($gff1_left<$reg_left_coord)&&($gff1_ref eq $ref)&&(!($gff1_type=~/transcript/))) {
#	    print STDERR "$gff1_ref\t$alg_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$gff1_type\t$gff1_left\n";
#	    print STDERR "$gff1_right\t$score_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
#	    print STDERR "$strand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
#	    print STDERR "$something_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
#	    print STDERR "$annotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
	    print OUTPUT "$gff1_ref\t$alg_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$gff1_type\t$gff1_left\t$gff1_right\t$score_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$strand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$something_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$annotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
	    $coords_sorted[$m][3]=1;
	}
    }
    print STDERR "Working on region $j\n";
    $parent="Parent=CT_region_".$j;
    $reg=$ref."\tmerge_transcripts_gff.pl\tCT_region\t".$reg_left_coord."\t".$reg_right_coord."\t.\t".$strand."\t.\tID=CT_region_".$j;
    print OUTPUT "$reg\n";
    for ($m = 0; $m <= $#coords_sorted; $m++) {
	$gff1_left=$coords_sorted[$m][0];
	$gff1_right=$coords_sorted[$m][1];
        $gff1_type=$coords_sorted[$m][2];
	$gff1_flag=$coords_sorted[$m][3];
	$gff1_strand=$coords_sorted[$m][4];
	$gff1_ref=$coords_sorted[$m][5];
#	print STDERR "$reg_left_coord\t$gff1_left\t$gff1_right\t$reg_right_coord\t$gff1_flag\n";
	if (($reg_left_coord<$gff1_left)&&($gff1_right<$reg_right_coord)&&($gff1_strand eq $strand)&&($gff1_ref eq $ref)) {
	    $k=0;
	    if($gff1_type=~/transcript/) {
		@annot = "";
		$feature="ncRNA";
		$coords_sorted[$m][3]=1;
		for ($n = 0; $n <= $#coords_sorted; $n++) {
		    $gff2_left=$coords_sorted[$n][0];
		    $gff2_right=$coords_sorted[$n][1];
		    $gff2_type=$coords_sorted[$n][2];
#		    $gff2_flag=$coords_sorted[$n][3];
		    $gff2_strand=$coords_sorted[$n][4];
		    $gff2_ref=$coords_sorted[$n][5];
		    if (($gff1_left<$gff2_left)&&($gff2_right<$gff1_right)&&($gff2_ref eq $ref)) {
			if (($gff2_type=~/gene/)&&($gff2_strand eq $gff1_strand)&&($gff2_ref eq $gff1_ref)) {
			    $gff2_right_UTR=$gff2_left-1;
			    if (($k==0) && ($strand_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}=~/\+/)&&($gff2_ref eq $ref)) {
				$firstUTR=$gff1_ref."\tmerge_transcripts_gff.pl\tfive_prime_UTR\t".$gff1_left."\t".$gff2_right_UTR."\t.\t".$tstrand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t.\t".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."_5UTR;Parent=".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref};
			    } elsif (($k==0) && ($strand_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}=~/\-/)&&($gff2_ref eq $ref)) {
				$firstUTR=$gff1_ref."\tmerge_transcripts_gff.pl\tthree_prime_UTR\t".$gff1_left."\t".$gff2_right_UTR."\t.\t".$tstrand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t.\t".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."_3UTR;Parent=".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref};
			    }
			    $coords_sorted[$n][3]=1;
			    $k=1;
			    $gene=$gff2_ref."\t".$alg_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$gff2_type."\t".$gff2_left."\t".$gff2_right."\t".$score_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$strand_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$something_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$annotation{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}.";Parent=".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref};
			    push (@annot, $gene);
			    $feature="mRNA";
			    $save_right=$gff2_right+1;
			} elsif (($gff2_type=~/CDS/)&&($gff2_strand eq $gff1_strand)&&($gff2_ref eq $gff1_ref)) {
			    $CDS=$gff2_ref."\t".$alg_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$gff2_type."\t".$gff2_left."\t".$gff2_right."\t".$score_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$strand_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$something_hash{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref}."\t".$annotation{$gff2_type}{$gff2_left}{$gff2_right}{$gff2_ref};
			    push (@annot, $CDS);
			    $coords_sorted[$n][3]=1;
			}
		    }
		}
		if (($k==1) && ($tstrand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}=~/\+/)&&($gff1_ref eq $ref)) {
		    $secondUTR=$gff1_ref."\tmerge_transcripts_gff.pl\tthree_prime_UTR\t".$save_right."\t".$gff1_right."\t.\t".$tstrand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t.\t".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."_3UTR;Parent=".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref};
		} elsif (($k==1) && ($tstrand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}=~/\-/)&&($gff1_ref eq $ref)) {
		    $secondUTR=$gff1_ref."\tmerge_transcripts_gff.pl\tfive_prime_UTR\t".$save_right."\t".$gff1_right."\t.\t".$tstrand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t.\t".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."_5UTR;Parent=".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref};
		}
		$transcript=$gff1_ref."\t".$talg_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t".$feature."\t".$gff1_left."\t".$gff1_right."\t".$tscore_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t".$tstrand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t".$tsomething_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}."\t".$tannotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}.";".$parent;
		print OUTPUT "$transcript\n";
		if ($feature eq "mRNA") {
		    print OUTPUT "$firstUTR";
		    print OUTPUT join ("\n", @annot);
		    print OUTPUT "\n$secondUTR\n";
		}
	    } elsif (($gff1_flag==0)&&($gff1_ref eq $ref)) {
		print OUTPUT "$gff1_ref\t$alg_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$gff1_type\t$gff1_left\t$gff1_right\t$score_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$strand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$something_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$annotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
		$coords_sorted[$m][3]=1;
	    }
	}
    }
}
for ($m = 0; $m <= $#coords_sorted; $m++) {
    $gff1_left=$coords_sorted[$m][0];
    $gff1_right=$coords_sorted[$m][1];
    $gff1_type=$coords_sorted[$m][2];
    $gff1_flag=$coords_sorted[$m][3];
    $gff1_strand=$coords_sorted[$m][4];
    $gff1_ref=$coords_sorted[$m][5];
    if (($gff1_flag==0)&&(!($gff1_type=~/transcript/))) {
#	    print STDERR "$gff1_right\t$score_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
#	    print STDERR "$strand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
#	    print STDERR "$something_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
#	    print STDERR "$annotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
	print OUTPUT "$gff1_ref\t$alg_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$gff1_type\t$gff1_left\t$gff1_right\t$score_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$strand_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$something_hash{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\t$annotation{$gff1_type}{$gff1_left}{$gff1_right}{$gff1_ref}\n";
	$coords_sorted[$m][3]=1;
    }
}

close (OUTPUT) || die "Cannot close $output : $!\n";

__END__
    
