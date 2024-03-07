# Predicting Bacterial Transcripts

Analysis run working on finalizing manuscript 1st submission

Run in tcsh; includes tcsh specific syntax that won't work in bash

# Table of Contents

[Download Data](https://github.com/jdhotopp/RDOPERON/blob/main/PUBLIC_COMMANDS.md#download-data)<br />

[Analysis Methods](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#analysis-methods)<br />
- [ ] [Minimap2](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#run-minimap2-on-all-fastq-files)<br />
- [ ] [Filter and sort](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#filter-and-sort-bam-file)<br />
- [ ] [Merge bam files](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#merge-all-bam-files)<br />
- [ ] [Generate readtable bed file](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#generate-readtable-bed-file)<br />
- [ ] [Filter on minimum number of reads](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#filtering-on-minimum-number-of-reads-in-this-case-20-reads)<br />
- [ ] [Stats on regions](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#stats-on-regions-for-table-1)<br />
- [ ] [Run tp.py](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#run-tppy)<br />
- [ ] [Counts for transcripts](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#counts-for-transcripts-predicted-for-table-1)<br />
- [ ] [Stats on transcirpts](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#stats-on-transcripts-for-table-1)<br />
- [ ] [Merge GFF files and identify UTRs](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#merge-gff-files-and-identify-utrs)<br />
- [ ] [Stats on merged GFF files](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#stats-on-merged-gff-files-for-table-1)<br />

[Organism Specific Results and Commands](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#organism-specific-results-and-commands)<br />

- [ ] [_E. coli_ K12](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#e-coli-k12-1)<br />
- [ ] [_E. coli_ E2348/69](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#e-coli-e234869-1)<br />
- [ ] [_Listeria monocytogenes_ Scott A](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#listeria-monoctyogenes-scott-a-cm0011591)<br />
- [ ] [_Listeria monocytogenes_ RO15](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#listeria-monocytogenes-ro15-cadehj0000000001-1)<br />
- [ ] [_Pseudomonas aeruginosa_ SG17M](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#pseudomonas-aeruginosa-sg17m-1)
- [ ] [_Pseudomonas aeruginosa_ NN2](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#pseudomonas-aeruginosa-nn2-1)<br />
- [ ] [_Haloferax volcanii_](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#haloferax-volcanii-1)<br />

[Plotting](https://github.com/jdhotopp/RDOPERON/blob/9e01af349e5f36570785c7d9c85d83116b3f8542/PUBLIC_COMMANDS.md#plotting)<br />

# Software and versions used

```
minimap2 version 2.24-r1122
samtools version 1.11 using htslib 1.11
bedtools version v2.27.1
perl version v5.30.2
python3 version 3.6.8, unless otherwise specified as version 3.11.4
seqkit version 0.7.2
bedops version 2.4.36 (convert2bed)
salmon version 1.10.2
Rscript version 3.6.3
bamtools version 2.5.1
stringtie version 1.3.4d
tama versión_date_2020_12_14
cupcake (downloaded and installed Nov 9 2022)
```

# Download Data

## _E. coli_ K12

REF_FILE=GCF_000005845.2_ASM584v2_genomic.fna (NC_000913.3)

REF_GFF_FILE=GCF_000005845.2_ASM584v2_genomic.gff

SRA files: SRR18061005, SRR18061002, SRR27982845, SRR23886068, SRR27982844, SRR27982843, SRR27982842, SRR27982841, SRR27982840

## _E. coli_ E2348/69

REF_FILE=GCF_014117345.2_ASM1411734v2_genomic.fna

REF_GFF_FILE=GCF_014117345.2_ASM1411734v2_genomic.gff

SRA files: SRR18061003, SRR18061004

## _Listeria monocytogenes_ Scott A (CM001159)

https://www.researchsquare.com/article/rs-1530110/v1

REF_FILE=CM001159.fsa

REF_GFF_FILE=CM001159.gff

SRA files: RR9606780, ERR9606777

## _Listeria monocytogenes_ RO15 (CADEHJ000000000.1)

https://www.researchsquare.com/article/rs-1530110/v1

SRA: ERR9606554, ERR9606779

REF_FILE=CADEHJ01.1.fsa_nt

REF_GFF_FILE=CADEHJ01.1_fixed_ref_name.gff3

## _Pseudomonas aeruginosa_ SG17M

https://journals.asm.org/doi/full/10.1128/jb.00418-21

SRA: ERR7451774,ERR7451775,ERR7451777,ERR7451778

REF_FILE=GCF_020978345.1_ASM2097834v1_genomic.fna

## _Pseudomonas aeruginosa_ NN2

https://journals.asm.org/doi/full/10.1128/jb.00418-21

REF_FILE: GCF_900185255.1_NN2_genomic.fna

REF_GFF_FILE: GCF_900185255.1_NN2_genomic.gff

SRA: ERR7451779, ERR7451780, ERR7451781, ERR7451782, ERR7451783, ERR7451784, ERR7451785 

## _Haloferax volcanii_

https://journals.asm.org/doi/full/10.1128/jb.00418-21

`curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000025685.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_000025685.1.zip" -H "Accept: application/zip"`

REF_FILE: GCF_000025685.1_ASM2568v1_genomic.fna

REF_GFF_FILE: GCF_000025685.1_ASM2568v1_genomic.gff

SRA: SRR11991303, SRR11991304, ERR7451781, SRR11991308, SRR11991309 

# Analysis methods

## Run minimap2 on all fastq files 

```
FASTQ_FILE=FILE.fastq
SAM_FILE=FILE.sam
BAM_FILE=FILE.bam
minimap2 -ax map-ont -t 2 "$REF_FILE" "$FASTQ_FILE" > "$SAM_FILE"
```

## Filter and sort bam file

`samtools view -bhF 2308 $SAM_FILE | samtools sort -o $BAM_FILE`

## Merge all bam files

`samtools merge merged_all.bam $BAM_FILE.1 $BAM_FILE.2 ... $BAM_FILE.n`

## Generate readtable bed file

`bamToBed -i merged_all.bam > merged_all.readtable.bed`

## Filtering on minimum number of reads (in this case 20 reads)

`bedtools merge -i merged_all.readtable.bed -s -c 6,4 -o distinct,count | awk '$5>20 {print $0}' > merged_all.readtable_20reads.merged.bed` 

## Stats on regions for Table 1

    foreach line (/`ls merged_all.readtable_20reads.merged.bed/`)
    echo "Number of regions="
    wc $line
    echo "Number of regions on (+) strand="
    perl -ne 'chomp;@line=split;$diff=$line[2]-$line[1];if ($line[3]=~/\+/) {$sum=$sum+$diff;print "$_\t$diff\t$sum\n"}' $line | wc
    echo "Number of regions on (-) strand="
    perl -ne 'chomp;@line=split;$diff=$line[2]-$line[1];if ($line[3]=~/\-/) {$sum=$sum+$diff;print "$_\t$diff\t$sum\n"}' $line | wc
    echo "Span of regions in bp on (+) strand="
    perl -ne 'chomp;@line=split;$diff=$line[2]-$line[1];if ($line[3]=~/\+/) {$sum=$sum+$diff;print "$_\t$diff\t$sum\n"}' $line | tail -1
    echo "Span of regions in bp on (-) strand="
    perl -ne 'chomp;@line=split;$diff=$line[2]-$line[1];if ($line[3]=~/\-/) {$sum=$sum+$diff;print "$_\t$diff\t$sum\n"}' $line | tail -1
    echo "Average span of regions in bp on (+) strand="
    perl -ne 'chomp;@line=split;$diff=$line[2]-$line[1];if ($line[3]=~/\+/) {$a++;$average=$sum/$a;$sum=$sum+$diff;print "$_\t$diff\t$sum\t$average\n"}' $line | tail -1
    echo "Average span of regions in bp on (-) strand="
    perl -ne 'chomp;@line=split;$diff=$line[2]-$line[1];if ($line[3]=~/\-/) {$a++;$average=$sum/$a;$sum=$sum+$diff;print "$_\t$diff\t$sum\t$average\n"}' $line | tail -1
    end

## Run tp.py

`python3 tp.py -r merged_all.readtable.bed -d merged_all.readtable_20reads.merged.bed --max_depth=25000 -o transcripts`

## Counts for transcripts predicted for Table 1

`head -2 transcripts.gff3`

## Stats on transcripts for Table 1

    foreach line (`ls transcripts.gff3 `)
    echo "Number of transcripts="
    grep -v "^#" $line | wc
    echo "Number of transcripts on (+)-strand="
    perl -ne '@line=split;if($line[6] eq "+") {print $_}' $line | wc
    echo "Number of transcripts on (-)-strand="
    perl -ne '@line=split;if($line[6] eq "-") {print $_}' $line | wc
    echo "Number of regions with just 1 transcript="
    sed 's/_T\w*//g' $line | awk '{print $9}' | sort | uniq -c | sort -n | awk '{print $1}' | sort -n | uniq -c | head -1
    echo "Maximum number of transcripts predicted per region="
    sed 's/_T\w*//g' $line | awk '{print $9}' | sort | uniq -c | sort -n | awk '{print $1}' | sort -n | uniq -c | tail -1
    end

## Merge GFF files and identify UTRs

`merge_transcripts_gff_v2.pl -g $REF_GFF_FILE -r merged_all.readtable_20reads.merged.bed -t transcripts.gff3 -o merge.gff`

## Stats on merged GFF files for Table 1

    foreach line (`ls merge.gff`)
    echo "Number of 3'-UTRs ="
    grep three $line | grep UTR | wc
    echo "Mean 3'-UTR in bp = "
    grep three $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$mean\n"' | tail -1
    echo "Minimum 3'-UTR in bp ="
    grep three $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$size\n"' | sort -n | head -1
    echo "Maximum 3'-UTR in bp ="
    grep three $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$size\n"' | sort -nr | head -1
    echo "Instances of the mode of the 3'-UTR and the mode in bp = "
    grep three $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$size\n"' | sort -nr | uniq -c | sort -n | tail -1
    echo "Number of 5'-UTRs ="
    grep five $line | grep UTR | wc
    echo "Mean 5'-UTR="
    grep five $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$mean\n"' | tail -1
    echo "Minimum 5'-UTR in bp = "
    grep five $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$size\n"' | sort -n | head -1
    echo "Maximum 5'-UTR in bp = "
    grep five $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$size\n"' | sort -nr | head -1
    echo "Number of instances of the mode and the Mode of the 5'-UTR in bp ="
    grep five $line | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; $sum=$sum+$size; $a++; $mean=$sum/$a; print "$size\n"' | sort -nr | uniq -c | sort -n | tail -1
    echo "Number of genes total = "
    perl -ne '@line=split;if($line[2] eq gene) {print $_}' $line | perl -ne 'if($_=~/ID=(.*);Name/) {print "$1\n"}' | sort | uniq | wc
    echo "Number of genes in the annotation file that are in an annotated transcript ="
    perl -ne '@line=split;if($line[2] eq gene) {print $_}' $line | grep Parent= | perl -ne 'if($_=~/ID=(.*);Name/) {print "$1\n"}' | sort | uniq |wc
    echo "Number of genes that aren't in an annotated transcript ="
    perl -ne '@line=split;if($line[2] eq gene) {print $_}' $line | grep -v Parent= | sort | uniq | wc
    echo "Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)"
    perl -ne '@line=split;if($line[2] eq gene) {print $_}' $line | grep Parent= | perl -ne 'if($_=~/ID=(.*);Name/) {print "$1\n"}' | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c | sort -n -k2
    echo "Distribution of the number of genes (col 1) in a number of transcripts (col 2)"
    perl -ne '@line=split;if($line[2] eq gene) {print $_}' $line | grep Parent= | perl -ne 'if($_=~/Parent=(.*)$/) {print "$1\n"}' | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c
    echo "Number of predicted mRNAs = "
    perl -ne '@line=split;if($line[2] eq mRNA) {print $_}' $line | wc
    echo "Average predicted mRNA in bp = "
    perl -ne '@line=split;if($line[2] eq mRNA) {$size=$line[4]-$line[3];$sum=$sum+$size;$a++;$mean=$sum/$a;print "$mean\n"}' $line | tail -1
    echo "Largest mRNA in bp ="
    perl -ne '@line=split;if($line[2] eq mRNA) {$size=$line[4]-$line[3];print "$size\n"}' $line | sort -n | tail -1
    echo "Smallest mRNA in bp = "
    perl -ne '@line=split;if($line[2] eq mRNA) {$size=$line[4]-$line[3];print "$size\n"}' $line | sort -nr | tail -1
    echo "Number of predicted ncRNAs = "
    perl -ne '@line=split;if($line[2] eq ncRNA) {print $_}' $line | wc
    echo "Average predicted ncRNA in bp ="
    perl -ne '@line=split;if($line[2] eq ncRNA) {$size=$line[4]-$line[3];$sum=$sum+$size;$a++;$mean=$sum/$a;print "$mean\n"}' $line | tail -1
    echo "Largest ncRNA in bp = "
    perl -ne '@line=split;if($line[2] eq ncRNA) {$size=$line[4]-$line[3];print "$size\n"}' $line | sort -n | tail -1
    echo "Smallest ncRNA in bp = "
    perl -ne '@line=split;if($line[2] eq ncRNA) {$size=$line[4]-$line[3];print "$size\n"}' $line | sort -nr | tail -1
    end

# Organism Specific Results and Commands 

## _E. coli_ K12

### Table 1A 

##### Filter for primary alignments, remove soft clip regions, and calculate stats for N50 and max read length
```
foreach line (`cat K12_list.txt`)
samtools view -bhF 2308 -o $line.bam $line.sam &
java -Xmx10g -jar jvarkit/dist/biostar84452.jar -o $line.nosoftclip.bam $line.bam &
samtools fastq $line.nosoftclip.bam > $line.nosoftclip.bam.fastq &
seqkit stats -T -a -o $line.nosoftclip.bam.stats < $line.nosoftclip.bam.fastq &
end
```

##### Calculate the number of reads that map to the rRNA

`awk '$3 ~ /rRNA/ { print }' ../020624_Ecoli_K12_merge.gff > gff_RNA`

`bedops-2.4.36/convert2bed --input=gff --output=bed -d < gff_RNA > rRNA.bed`

```
foreach line ( `ls *.nosoftclip.bam` )
echo $line
samtools view -cL rRNA.bed $line
end
```

### Stats on regions for Table 1

```
Number of regions=
   1055  5275 35092 merged_all.readtable_20reads.merged.bed
Number of regions on (+) strand=
   521    3647   23853
Number of regions on (-) strand=
   534    3738   24397
Span of regions in bp on (+) strand=
NC_000913.3     4640403 4641652 +       25      1249    2068709
Span of regions in bp on (-) strand=
NC_000913.3     4639466 4640788 -       581     1322    2135707
Average span of regions in bp on (+) strand=
NC_000913.3     4640403 4641652 +       25      1249    2068709 3968.25335892514
Average span of regions in bp on (-) strand=
NC_000913.3     4639466 4640788 -       581     1322    2135707 3996.97565543071
```

### Counts for transcripts predicted for Table 1

```
##gff-version 3
NC_000913.3     tp.py   transcript      159     306     1.0     +       .       ID=R_0_T3
```

### Stats on transcripts for Table 1

```
Number of transcripts=
   3618   32562  240558
Number of transcripts on (+)-strand=
   1465   13185   96859
Number of transcripts on (-)-strand=
   2153   19377  143699
Number of regions with just 1 transcript=
    389 1
Maximum number of transcripts predicted per region=
      1 254
```

### Stats on merged GFF files for Table 1 and Figure 2

```
Number of 3'-UTRs =
   2484   22356  280870
Mean 3'-UTR in bp =
149.545088566828
Minimum 3'-UTR in bp =
1
Maximum 3'-UTR in bp =
2715
Instances of the mode of the 3'-UTR and the mode in bp =
     56 36
Number of 5'-UTRs =
   2484   22356  278378
Mean 5'-UTR=
133.164251207729
Minimum 5'-UTR in bp =
1
Maximum 5'-UTR in bp =
2121
Number of instances of the mode and the Mode of the 5'-UTR in bp =
     67 14
Number of genes total =
   4494    4494  282007
Number of genes in the annotation file that are in an annotated transcript =
   2357    2357  148597
Number of genes that aren't in an annotated transcript =
   2137   19233  448465
Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)
   1341 1
    670 2
    182 3
     72 4
     38 5
     21 6
      3 7
     10 8
      7 9
      6 10
      4 11
      1 12
      1 13
      1 15
Distribution of the number of genes (col 1) in a number of transcripts (col 2)
   1564 1
    535 2
    216 3
     87 4
     38 5
     17 6
     10 7
      4 8
      4 9
      1 10
      2 12
      1 13
      1 15
      1 17
Number of predicted mRNAs =
   2484   22356  201803
Average predicted mRNA in bp =
1617.98711755233
Largest mRNA in bp =
13305
Smallest mRNA in bp =
131
Number of predicted ncRNAs =
   1233   11374  113071
Average predicted ncRNA in bp =
516.680454176805
Largest ncRNA in bp =
2947
Smallest ncRNA in bp =
52
```

#### Calculating the median 5'-UTR for Table 1

`grep five Ecoli_K12_merge.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -1243 | tail -2`

```
 53
 53
```

#### Calculating the median 3'-UTR for Table 1

`grep three Ecoli_K12_merge.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -1243 | tail -2`

```
 72
 72
```

#### Composition of the largest transcript for Table 1

`grep mRNA Ecoli_K12_merge.gff | perl -ne 'chomp;@line=split;$diff=$line[4]-$line[3];if ($diff >13000) {print "$_\n"}'`

```
NC_000913.3     tp.py   mRNA    2099804 2113109 0.23    -       .       ID=R_779/1_T43;Parent=CT_region_505
```

`grep gene Ecoli_K12_merge.gff | grep ID=R_779/1_T43 | perl -ne 'if ($_=~/gene=(\w*);gene/) { print "$1\n"}' | sort`

```
glf
gnd
insH7
rfbA
rfbB
rfbC
rfbD
rfbX
wbbH
wbbI
wbbJ
wbbK
wbbL
```

### Generated LB/DMEM strand-specific bam files for viewing in IGV and for plotting below

`samtools merge merged_all_LB.bam FILE.1.bam FILE.2.bam ... FILE.n.bam`

`samtools merge merged_all_DMEM.bam FILE.1.bam FILE.2.bam ... FILE.n.bam`

`samtools view -b -F 16 merged_all_LB.bam > merged_all_LB.f.bam`

`samtools view -b -f 16 merged_all_LB.bam > merged_all_LB.r.bam`

`samtools view -b -F 16 merged_all_DMEM.bam > merged_all_DMEM.f.bam`

`samtools view -b -f 16 merged_all_DMEM.bam > merged_all_DMEM.r.bam`

`bamToBed -i merged_all_LB.f.bam > merged_all_LB.f.bed`

`bamToBed -i merged_all_LB.r.bam > merged_all_LB.r.bed`

`bamToBed -i merged_all_DMEM.f.bam > merged_all_DMEM.f.bed`

`bamToBed -i merged_all_DMEM.r.bam > merged_all_DMEM.r.bed`

## _E. coli_ E2348/69

### Table 1A 

##### Filter for primary alignments, remove soft clip regions, and calculate stats for N50 and max read length
```
foreach line (`ls E*.bam | sed 's/.bam//' `)
samtools view -bhF 2308 -o $line.bam $line.sam &
java -Xmx10g -jar jvarkit/dist/biostar84452.jar -o $line.nosoftclip.bam $line.bam &
samtools fastq $line.nosoftclip.bam > $line.nosoftclip.bam.fastq &
seqkit stats -T -a -o $line.nosoftclip.bam.stats < $line.nosoftclip.bam.fastq &
end
```

##### Calculate the number of reads that map to the rRNA

`awk '$3 ~ /rRNA/ { print }' ../020624_E2348_69.v3.gff > gff_RNA_E2348_69`

`bedops-2.4.36/convert2bed --input=gff --output=bed -d < gff_RNA_E2348_69 > rRNA_E2348_69.bed`

foreach line ( `ls E*.bam` )
echo $line
samtools view -cL rRNA_E2348_69.bed $line
end

### Stats on regions for Table 1

```
Number of regions=
 1071  5355 37490 merged_all_E2348_69.readtable_20reads.merged.bed
Number of regions on (+) strand=
    528    3696   25066
Number of regions on (-) strand=
    543    3801   25675
Span of regions in bp on (+) strand=
NZ_CP059840.2   4930520 4932974 +       109     2454    1951551
Span of regions in bp on (-) strand=
NZ_CP059840.2   4939925 4944414 -       1295    4489    1827581
Average span of regions in bp on (+) strand=
NZ_CP059840.2   4930520 4932974 +       109     2454    1951551 3777.31976744186
Average span of regions in bp on (-) strand=
NZ_CP059840.2   4939925 4944414 -       1295    4489    1827581 3446.29867674858
```

### Counts for transcripts predicted for Table 1

```
##gff-version 3
NZ_CP059840.2   tp.py   transcript      1       2662    1.0     +       .       ID=R_0_T1
```

### Stats on transcripts for Table 1

```
Number of transcripts=
   2248   20232  152182
Number of transcripts on (+)-strand=
   1101    9909   74207
Number of transcripts on (-)-strand=
   1147   10323   77975
Number of regions with just 1 transcript=
    429 1
Maximum number of transcripts predicted per region=
      1 141
```

### Stats on merged GFF files for Table 1

```
Number of 3'-UTRs =
   1843   16587  210241
Mean 3'-UTR in bp =
127.38741182854
Minimum 3'-UTR in bp =
1
Maximum 3'-UTR in bp =
1452
Instances of the mode of the 3'-UTR and the mode in bp =
     47 43
Number of 5'-UTRs =
   1843   16587  208404
Mean 5'-UTR=
120.507867607162
Minimum 5'-UTR in bp =
1
Maximum 5'-UTR in bp =
2817
Number of instances of the mode and the Mode of the 5'-UTR in bp =
     50 13
Number of genes total =
   4809    4809   91371
Number of genes in the annotation file that are in an annotated transcript =
   2034    2034   38646
Number of genes that aren't in an annotated transcript =
   2775   24975  421928
Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)
   1299 1
    494 2
    157 3
     39 4
     22 5
      5 6
      4 7
      4 8
      4 9
      3 10
      2 11
      1 12
Distribution of the number of genes (col 1) in a number of transcripts (col 2)
   1101 1
    412 2
    171 3
     70 4
     39 5
     22 6
      9 7
      2 8
      1 9
      3 10
      1 11
      2 13
      1 14
Number of predicted mRNAs =
   1843   16587  152298
Average predicted mRNA in bp =
1733.01899077591
Largest mRNA in bp =
15256
Smallest mRNA in bp =
129
Number of predicted ncRNAs =
    407    3672   34463
Average predicted ncRNA in bp =
645.159705159705
Largest ncRNA in bp =
2916
Smallest ncRNA in bp =
80
```

#### Calculating the median 5'-UTR for Table 1

`grep five E2348_69.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -922 | tail -1`

```
49
```

#### Calculating the median 3'-UTR for Table 1

`grep three E2348_69.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -922 | tail -1`

```
62
```

### Composition of the largest transcript for Table 1

`grep mRNA E2348_69.gff | perl -ne 'chomp;@line=split;$diff=$line[4]-$line[3];if ($diff >15000) {print "$_\n"}'`

```
NZ_CP059840.2   tp.py   mRNA    1654359 1669615 0.25    +       .       ID=R_178_T14;Parent=CT_region_338

```

`grep gene E2348_69.gff | grep ID=R_178_T14 | perl -ne 'if ($_=~/gene=(\w*);gene/) { print "$1\n"}' | sort`

```
nuoA
nuoB
nuoC
nuoE
nuoF
nuoG
nuoH
nuoI
nuoJ
nuoK
nuoL
nuoM
nuoN
```

### Generated LB/DMEM strand-specific bam files for viewing in IGV and for plotting below

`samtools view -b -F 16 E2348_69_LB.bam > E2348_69_LB.f.bam`

`samtools view -b -f 16 E2348_69_LB.bam > E2348_69_LB.r.bam`

`samtools view -b -F 16 E2348_69_DMEM.bam > E2348_69_DMEM.f.bam`

`samtools view -b -f 16 E2348_69_DMEM.bam > E2348_69_DMEM.r.bam`

`bamToBed -i E2348_69_LB.f.bam > E2348_69_LB.f.bed`

`bamToBed -i E2348_69_LB.r.bam > E2348_69_LB.r.bed`

`bamToBed -i E2348_69_DMEM.f.bam > E2348_69_DMEM.f.bed`

`bamToBed -i E2348_69_DMEM.r.bam > E2348_69_DMEM.r.bed`

### _glmY_ and _glmZ_

#### glmY: (plus strand, 1260087-1260262 bp) 

`grep 1260 E2348_69_transcripts.bed.bed`

```
NZ_CP059840.2   1260086 1260262 R_138_T18       .       +
```

`cutFasta -f NZ_CP059840.2 -x 1260087 -y 1260262 $REF_FILE`

Tried a BLASTN search against NT but useless as 100% identical to a ridiculous number of genomes.  But manually checked against the sequences on EcoCyc, confirming that they are attributed correctly. From EcoCyc (on 02/08/24):
>gnl|ECOLI|TKE1-RNA gn=glmY small regulatory RNA GlmY (complement(2691157..2691340)) Escherichia coli K-12 substr. MG1655
AGUGGCUCAU UCACCGACUU AUGUCAGCCC CUUCGGGACG UGCUACAUAA AAUACGAAUG
ACGCACAACA AGGUGCCUGC CGUCCAACUU CUGAUAUCAG CGUAGCUAUA UCAACCAUCG
GGCGAAACGU CGAGUUAGGC ACCGCCUUAU UCCAUAACAA AGCCGGGUAA UUCCCGGCUU
UGUU

#### glmZ: (minus strand, 4843198-4843393 bp)

`grep 484319 E2348_69_transcripts.bed.bed`

```
NZ_CP059840.2   4843197 4843393 R_1034/1_T0     .       -
```

`cutFasta -f NZ_CP059840.2 -x 4843198 -y 4843393 $REF_FILE`

Tried a BLASTN search against NT but useless as 100% identical to a ridiculous number of genomes.  But manually checked against the sequences on EcoCyc, confirming that they are attributed correctly. From EcoCyc (on 02/08/24): (Have to Reverse Complement though)
>gnl|ECOLI|SRAJ-RNA gn=glmZ small regulatory RNA GlmZ 3986432..3986603 Escherichia coli K-12 substr. MG1655
GUAGAUGCUC AUUCCAUCUC UUAUGUUCGC CUUAGUGCCU CAUAAACUCC GGAAUGACGC
AGAGCCGUUU ACGGUGCUUA UCGUCCACUG ACAGAUGUCG CUUAUGCCUC AUCAGACACC
AUGGACACAA CGUUGAGUGA AGCACCCACU UGUUGUCAUA CAGACCUGUU UU

### Differential expression analysis for Figure 5

#### Download from SRA

SRA: ERR3932853, ERR3932854, ERR3932855, ERR3932847, ERR3932848, ERR3932849

#### Run Salmon

`module load salmon/1.10.2`

`perl -ne 'chomp;@line=split;$a=$line[1]+1;print "$line[3]\t$a\t$line[2]\t$line[0]\n"' ../020624_E2348_69_transcripts.bed.bed | sed 's/\//_/g' > instructions_transcript_fasta`

`cutFasta -i instructions_transcript_fasta ../ncbi_dataset/data/GCF_014117345.2/GCF_014117345.2_ASM1411734v2_genomic.fna > E2348_predicted_transcripts.fsa`

`salmon index -i E2348_salmon_index -t E2348_predicted_transcripts.fsa`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932853_1.fastq -2 ERR3932853_2.fastq --validateMappings -o quants/ERR3932853`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932854_1.fastq -2 ERR3932854_2.fastq --validateMappings -o quants/ERR3932854`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932855_1.fastq -2 ERR3932855_2.fastq --validateMappings -o quants/ERR3932855`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932847_1.fastq -2 ERR3932847_2.fastq --validateMappings -o quants/ERR3932847`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932848_1.fastq -2 ERR3932848_2.fastq --validateMappings -o quants/ERR3932848`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932849_1.fastq -2 ERR3932849_2.fastq --validateMappings -o quants/ERR3932849`

Salmon detected they were all library type IU, which is unstranded.

`salmon quant -i E2348_salmon_index -l A -1 ERR3932853_1.fastq -2 ERR3932853_2.fastq --validateMappings -o quants/ERR3932853 &`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932854_1.fastq -2 ERR3932854_2.fastq --validateMappings -o quants/ERR3932854 &`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932855_1.fastq -2 ERR3932855_2.fastq --validateMappings -o quants/ERR3932855 &`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932847_1.fastq -2 ERR3932847_2.fastq --validateMappings -o quants/ERR3932847 &`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932848_1.fastq -2 ERR3932848_2.fastq --validateMappings -o quants/ERR3932848 &`

`salmon quant -i E2348_salmon_index -l A -1 ERR3932849_1.fastq -2 ERR3932849_2.fastq --validateMappings -o quants/ERR3932849 &`

# Need commands from Chris for running FADU

`minimap2 -ax map-ont -t 2 E2348_predicted_transcripts.fsa SRR18061003.fastq > SRR18061003.transcripts.bam`

`minimap2 -ax map-ont -t 2 E2348_predicted_transcripts.fsa SRR18061004.fastq > SRR18061004.transcripts.bam`

`salmon quant --ont --libType SF -a SRR18061003.transcripts.bam -t E2348_predicted_transcripts.fsa -o quants/SRR18061003 &`

`salmon quant --ont --libType SF -a SRR18061004.transcripts.bam -t E2348_predicted_transcripts.fsa -o quants/SRR18061004 &`

`sed 's/\//_/g' ../020624_E2348_69.v3.gff > 020624_E2348_69_just_HERE.gff`

`sed 's/#/_/g' quants/ERR3932853/quant.sf > quants/ERR3932853/quant.sf.fixed`

`sed 's/#/_/g' quants/ERR3932854/quant.sf > quants/ERR3932854/quant.sf.fixed`

`sed 's/#/_/g' quants/ERR3932855/quant.sf > quants/ERR3932855/quant.sf.fixed`

`sed 's/#/_/g' quants/ERR3932847/quant.sf > quants/ERR3932847/quant.sf.fixed`

`sed 's/#/_/g' quants/ERR3932848/quant.sf > quants/ERR3932848/quant.sf.fixed`

`sed 's/#/_/g' quants/ERR3932849/quant.sf > quants/ERR3932849/quant.sf.fixed`

`merge_transcripts_de_v2.pl -g 020624_E2348_69_just_HERE.gff -a FADU/ERR3932853.sorted.counts.txt -b FADU/ERR3932854.sorted.counts.txt -c FADU/ERR3932855.sorted.counts.txt -d FADU/ERR3932847.sorted.counts.txt -e FADU/ERR3932848.sorted.counts.txt -f FADU/ERR3932849.sorted.counts.txt -i quants/ERR3932853/quant.sf.fixed -j quants/ERR3932854/quant.sf.fixed -k quants/ERR3932855/quant.sf.fixed -l quants/ERR3932847/quant.sf.fixed -m quants/ERR3932848/quant.sf.fixed -n quants/ERR3932849/quant.sf.fixed -p quants/SRR18061003/quant.sf -q quants/SRR18061004/quant.sf -o 021524_E2348_69_6salmonIlluminashort_2salmonaONTlong_6FADUIlluminashort.txt &`

#### This table can be parsed to make Figure 5G

#### Figure 5 ABCDEF

`Rscript plot_alignments_prediction_seqid.r merged_all_E2348_69.readtable.bed merged_all_E2348_69.readtable.bed Figure5 Figure5 4728000 4735000 + - E2348_69_transcripts.bed.bed GCF_014117345.2_ASM1411734v2_genomic.gff NZ_CP059840.2 jpeg`

#### Figure 5H

`RScript Figure5H.R`

# Need Chris's script for Figure 5I

## _Listeria monoctyogenes_ Scott A (CM001159.1)

### Stats on regions for Table 1

```
Number of regions=
  525  2625 16732 merged_all_CM.readtable_20reads.merged.bed
Number of regions on (+) strand=
    238    1666   10303
Number of regions on (-) strand=
    287    2009   12563
Span of regions in bp on (+) strand=
CM001159.1      2931939 2934333 +       39      2394    703660
Span of regions in bp on (-) strand=
CM001159.1      3018219 3021560 -       169     3341    821637
Average span of regions in bp on (+) strand=
CM001159.1      2931939 2934333 +       39      2394    703660  2946.49579831933
Average span of regions in bp on (-) strand=
CM001159.1      3018219 3021560 -       169     3341    821637  2851.20557491289
```

### Counts for transcripts predicted for Table 1

```
##gff-version 3
CM001159.1      tp.py   transcript      4612    5065    1.0     +       .       ID=R_1_T0
```

### Stats for transcripts predicted for Table 1

```
Number of transcripts=
    881    7929   56480
Number of transcripts on (+)-strand=
    402    3618   25477
Number of transcripts on (-)-strand=
    479    4311   31003
Number of regions with just 1 transcript=
    218 1
Maximum number of transcripts predicted per region=
      1 32
```

### Stats on merged GFF files for Table 1

```
Number of 3'-UTRs =
    536    4824   59102
Mean 3'-UTR in bp =
121.744402985075
Minimum 3'-UTR in bp =
3
Maximum 3'-UTR in bp =
1306
Instances of the mode of the 3'-UTR and the mode in bp =
     24 38
Number of 5'-UTRs =
    536    4824   58564
Mean 5'-UTR=
143.389925373134
Minimum 5'-UTR in bp =
1
Maximum 5'-UTR in bp =
2303
Number of instances of the mode and the Mode of the 5'-UTR in bp =
     17 12
Number of genes total =
   3038    3038   50434
Number of genes in the annotation file that are in an annotated transcript =
    762     762   12623
Number of genes that aren't in an annotated transcript =
   2276   20497  322554
Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)
    634 1
    103 2
     16 3
      6 4
      2 5
      1 6
Distribution of the number of genes (col 1) in a number of transcripts (col 2)
    351 1
    103 2
     45 3
     19 4
      7 5
      3 6
      3 7
      1 8
      2 9
      1 22
      1 38
Number of predicted mRNAs =
    536    4824   42267
Average predicted mRNA in bp =
1660.26492537313
Largest mRNA in bp =
29034
Smallest mRNA in bp =
224
Number of predicted ncRNAs =
    345    3105   27572
Average predicted ncRNA in bp =
496.515942028986
Largest ncRNA in bp =
2585
Smallest ncRNA in bp =
95
```

#### Calculating the median 5'-UTR for Table 1

`grep five *.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n |head -269 | tail -2`

```
42
42
```

#### Calculating the median 3'-UTR for Table 1

`grep three *.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n |head -269 | tail -2`

```
49
49
```

### Calculations for Figure 2, made in Excel

`grep five GFF_FILE | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"'`

`grep three GFF_FILE | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"'`

### Composition of the largest transcript for Table 1

```
CM001159.1      tp.py   mRNA    101227  130261  0.69    +       .       ID=R_8/2_T3;Parent=CT_region_13`
```

`grep gene merged_CM001159.gff | grep ID=R_8/2_T3 | perl -ne 'if ($_=~/gene-(\w*);Name/) { print "$1\n"}' | sort`

```
LMOSA_9400
LMOSA_9410
LMOSA_9420
LMOSA_9430
LMOSA_9440
LMOSA_9450
LMOSA_9460
LMOSA_9470
LMOSA_9480
LMOSA_9490
LMOSA_9500
LMOSA_9510
LMOSA_9520
LMOSA_9530
LMOSA_9540
LMOSA_9550
LMOSA_9560
LMOSA_9570
LMOSA_9580
LMOSA_9590
LMOSA_9600
LMOSA_9610
LMOSA_9620
LMOSA_9630
LMOSA_9640
LMOSA_9650
LMOSA_9660
LMOSA_9670
LMOSA_9680
LMOSA_9690
LMOSA_9700
LMOSA_9710
LMOSA_9720
LMOSA_9730
LMOSA_9740
LMOSA_9750
LMOSA_9760
LMOSA_9770
```

## _Listeria monocytogenes_ RO15 (CADEHJ000000000.1)

### Stats on regions for Table 1

```
Number of regions=
  464  2320 18038 merged_all_CADEH.readtable_20reads.merged.bed
Number of regions on (+) strand=
    206    1442   10347
Number of regions on (-) strand=
    258    1806   13093
Span of regions in bp on (+) strand=
CADEHJ010000001.1       2956519 2958908 +       33      2389    589005
Span of regions in bp on (-) strand=
CADEHJ010000001.1       3038687 3042028 -       167     3341    759698
Average span of regions in bp on (+) strand=
CADEHJ010000001.1       2956519 2958908 +       33      2389    589005  2847.65048543689
Average span of regions in bp on (-) strand=
CADEHJ010000001.1       3038687 3042028 -       167     3341    759698  2931.61627906977
```

### Counts for transcripts for Table 1

```
##gff-version 3
CADEHJ010000001.1       tp.py   transcript      4295    4748    1.0     +       .       ID=R_1_T1
```

### Stats on transcripts for Table 1

```
Number of transcripts=
    793    7137   56368
Number of transcripts on (+)-strand=
    361    3249   25395
Number of transcripts on (-)-strand=
    432    3888   30973
Number of regions with just 1 transcript=
    199 1
Maximum number of transcripts predicted per region=
      1 31
```

### Stats for merged GFF for Table 1

```
Number of 3'-UTRs =
    491    4419   57481
Mean 3'-UTR in bp =
115.590631364562
Minimum 3'-UTR in bp =
1
Maximum 3'-UTR in bp =
1685
Instances of the mode of the 3'-UTR and the mode in bp =
     20 45
Number of 5'-UTRs =
    491    4419   56990
Mean 5'-UTR=
116.767820773931
Minimum 5'-UTR in bp =
1
Maximum 5'-UTR in bp =
2303
Number of instances of the mode and the Mode of the 5'-UTR in bp =
     18 12
Number of genes total =
   3149    3149   59831
Number of genes in the annotation file that are in an annotated transcript =
    677     677   12863
Number of genes that aren't in an annotated transcript =
   2472   22248  454280
Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)
    552 1
    104 2
     15 3
      4 4
      1 5
      1 7
Distribution of the number of genes (col 1) in a number of transcripts (col 2)
    317 1
     96 2
     41 3
     19 4
      6 5
      3 6
      3 7
      2 8
      2 9
      1 22
Number of predicted mRNAs =
    491    4419   42075
Average predicted mRNA in bp =
1606.95519348269
Largest mRNA in bp =
10791
Smallest mRNA in bp =
209
Number of predicted ncRNAs =
    303    2729   26518
Average predicted ncRNA in bp =
522.435643564356
Largest ncRNA in bp =
2588
Smallest ncRNA in bp =
136
```

#### Calculating the median 5'-UTR for Table 1

`grep five CADEH_merged.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -246 | tail -1`

```
33
```

#### Calculating the median 3'-UTR for Table 1

`grep three CADEH_merged.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -246 | tail -1`

```
47
```

## _Pseudomonas aeruginosa_ SG17M

### Stats on regions for Table 1

```
Number of regions=
  391  1955 13634 merged_all_SG17M.readtable_20reads.merged.bed
Number of regions on (+) strand=
    181    1267    8410
Number of regions on (-) strand=
    210    1470    9792
Span of regions in bp on (+) strand=
NZ_CP080369.1   6863494 6866263 +       22      2769    530329
Span of regions in bp on (-) strand=
NZ_CP080369.1   6868764 6873119 -       26      4355    589348
Average span of regions in bp on (+) strand=
NZ_CP080369.1   6863494 6866263 +       22      2769    530329  2914.69613259668
Average span of regions in bp on (-) strand=
NZ_CP080369.1   6868764 6873119 -       26      4355    589348  2785.68095238095
```

### Counts for transcripts for Table 1

```
##gff-version 3
NZ_CP080369.1   tp.py   transcript      124417  124698  1.0     +       .       ID=R_3/1_T1
```

### Stats for transcripts for Table 1

```
Number of transcripts=
    274    2466   18698
Number of transcripts on (+)-strand=
     79     711    5327
Number of transcripts on (-)-strand=
    195    1755   13371
Number of regions with just 1 transcript=
     85 1
Maximum number of transcripts predicted per region=
      1 68
```
   
### Stats for merged GFF for Table 1

```
Number of 3'-UTRs =
    132    1188   15100
Mean 3'-UTR in bp =
163.44696969697
Minimum 3'-UTR in bp =
2
Maximum 3'-UTR in bp =
2234
Instances of the mode of the 3'-UTR and the mode in bp =
      7 43
Number of 5'-UTRs =
    132    1188   14968
Mean 5'-UTR=
185.015151515152
Minimum 5'-UTR in bp =
1
Maximum 5'-UTR in bp =
1835
Number of instances of the mode and the Mode of the 5'-UTR in bp =
      4 95
Number of genes total =
   6349    6349  120631
Number of genes in the annotation file that are in an annotated transcript =
    208     208    3952
Number of genes that aren't in an annotated transcript =
   6141   55269 1092033
Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)
    167 1
     25 2
     14 3
      2 4
Distribution of the number of genes (col 1) in a number of transcripts (col 2)
     78 1
     25 2
     12 3
      8 4
      1 5
      2 6
      1 7
      4 8
      1 15
Number of predicted mRNAs =
    132    1188   10907
Average predicted mRNA in bp =
1601
Largest mRNA in bp =
14168
Smallest mRNA in bp =
186
Number of predicted ncRNAs =
    143    1289   12216
Average predicted ncRNA in bp =
572.363636363636
Largest ncRNA in bp =
6361
Smallest ncRNA in bp =
97
```

#### Calculating the median 5'-UTR for Table 1

`grep five 111023_SG17M.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -67 | tail -2`

```
93
93
```

#### Calculating the median 3'-UTR for Table 1

`grep three 111023_SG17M.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -67 | tail -2`

```
59
59
```

## _Pseudomonas aeruginosa_ NN2

### Stats on regions for Table 1

```
Number of regions=
 1209  6045 42254 merged_all_NN2.readtable_20reads.merged.bed
Number of regions on (+) strand=
    612    4284   28975
Number of regions on (-) strand=
    597    4179   28240
Span of regions in bp on (+) strand=
NZ_LT883143.1   6859572 6861564 +       54      1992    1944294
Span of regions in bp on (-) strand=
NZ_LT883143.1   6900332 6902966 -       115     2634    1886100
Average span of regions in bp on (+) strand=
NZ_LT883143.1   6859572 6861564 +       54      1992    1944294 3173.69607843137
Average span of regions in bp on (-) strand=
NZ_LT883143.1   6900332 6902966 -       115     2634    1886100 3154.88442211055
```

### Counts for transcripts for Table 1

```
##gff-version 3
NZ_LT883143.1   tp.py   transcript      40343   40864   0.64    +       .       ID=R_7_T4
```

### Stats for transcripts for Table 1

```
Number of transcripts=
   1103    9927   75660
Number of transcripts on (+)-strand=
    495    4455   33889
Number of transcripts on (-)-strand=
    608    5472   41771
Number of regions with just 1 transcript=
    258 1
Maximum number of transcripts predicted per region=
      1 63
```

### Stats on merged GFF files for Table 1

```
Number of 3'-UTRs =
    599    5391   69341
Mean 3'-UTR in bp =
236.838063439065
Minimum 3'-UTR in bp =
1
Maximum 3'-UTR in bp =
2809
Instances of the mode of the 3'-UTR and the mode in bp =
     16 48
Number of 5'-UTRs =
    599    5391   68742
Mean 5'-UTR=
204.791318864775
Minimum 5'-UTR in bp =
1
Maximum 5'-UTR in bp =
1942
Number of instances of the mode and the Mode of the 5'-UTR in bp =
     13 10
Number of genes total =
   6380    6380  114840
Number of genes in the annotation file that are in an annotated transcript =
    764     764   13752
Number of genes that aren't in an annotated transcript =
   5616   50544  976073
Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)
    571 1
    143 2
     33 3
      8 4
      4 5
      5 6
Distribution of the number of genes (col 1) in a number of transcripts (col 2)
    397 1
    117 2
     30 3
     19 4
     15 5
      5 6
      5 7
      4 8
      2 9
      1 10
      1 12
      1 14
      1 15
Number of predicted mRNAs =
    599    5391   50183
Average predicted mRNA in bp =
1739.78130217028
Largest mRNA in bp =
12709
Smallest mRNA in bp =
146
Number of predicted ncRNAs =
    505    4547   42867
Average predicted ncRNA in bp =
536.237623762376
Largest ncRNA in bp =
2851
Smallest ncRNA in bp =
77
```

#### Calculating the median 5'-UTR for Table 1

`grep five 111023_NN2.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -300 | tail -1`

```
85
```

#### Calculating the median 3'-UTR for Table 1

`grep three 111023_NN2.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -300 | tail -1`
```
78
```

## _Haloferax volcanii_

### Stats on regions for Table 1

```
Number of regions=
  640  3200 20648 merged_all.readtable_20reads.merged.bed
Number of regions on (+) strand=
    318    2226   13981
Number of regions on (-) strand=
    322    2254   14203
Span of regions in bp on (+) strand=
NC_013966.1     630775  631431  +       38      656     893429
Span of regions in bp on (-) strand=
NC_013966.1     633981  634926  -       95      945     974115
Average span of regions in bp on (+) strand=
NC_013966.1     630775  631431  +       38      656     893429  2807.46226415094
Average span of regions in bp on (-) strand=
NC_013966.1     633981  634926  -       95      945     974115  3022.26708074534
```

### Counts for transcripts for Table 1

```
##gff-version 3
NC_013964.1     tp.py   transcript      130869  133531  1.0     +       .       ID=R_1_T0
```

### Stats for transcripts for Table 1

```
Number of transcripts=
    613    5517   39857
Number of transcripts on (+)-strand=
    241    2169   15674
Number of transcripts on (-)-strand=
    372    3348   24183
Number of regions with just 1 transcript=
    226 1
Maximum number of transcripts predicted per region=
      1 27
```

### Stats on merged GFF files for Table 1

```
Number of 3'-UTRs =
    257    2313   28635
Mean 3'-UTR in bp =
180.824902723735
Minimum 3'-UTR in bp =
1
Maximum 3'-UTR in bp =
2039
Instances of the mode of the 3'-UTR and the mode in bp =
      9 31
Number of 5'-UTRs =
    257    2313   28378
Mean 5'-UTR=
389.494163424125
Minimum 5'-UTR in bp =
1
Maximum 5'-UTR in bp =
2955
Number of instances of the mode and the Mode of the 5'-UTR in bp =
      7 1
Number of genes total =
   3956    3956  154525
Number of genes in the annotation file that are in an annotated transcript =
    380     380   14830
Number of genes that aren't in an annotated transcript =
   3576   32184  672549
Distribution of the number of transcripts (col 1) associated with a number of genes (col 2)
    298 1
     64 2
     14 3
      3 4
      1 10
Distribution of the number of genes (col 1) in a number of transcripts (col 2)
    162 1
     49 2
     18 3
     11 4
      4 5
      5 6
      2 7
      1 8
      1 9
      1 10
      1 11
      2 15
Number of predicted mRNAs =
    257    2313   20537
Average predicted mRNA in bp =
1977.92996108949
Largest mRNA in bp =
10463
Smallest mRNA in bp =
136
Number of predicted ncRNAs =
    356    3204   28786
Average predicted ncRNA in bp =
722.890449438202
Largest ncRNA in bp =
3045
Smallest ncRNA in bp =
81
```

#### Calculating the median 5'-UTR for Table 1

`grep five 021324_transcripts_merge.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -129 | tail -1`

```
209
```

#### Calculating the median 3'-UTR for Table 1

`grep three 021324_transcripts_merge.gff | grep UTR | perl -ne 'chomp;@line=split;$size=$line[4]-$line[3]+1; print "$size\n"' | sort -n | head -129 | tail -1`

```
80
```

# Plotting

## _E. coli_ K12

### Troubleshooting - first 30 kbp plus and minus strand

`Rscript plot_alignments_prediction.r merged_all_LB.r.bed merged_all_DMEM.r.bed LB DMEM_0_30_minus 1 30000 - - transcripts.bed $REF_GFF_FILE jpeg`

`Rscript plot_alignments_prediction.r merged_all_LB.f.bed merged_all_DMEM.f.bed LB DMEM_0_30 1 30000 + + transcripts.bed $REF_GFF_FILE jpeg`

### Testing detection of a known small RNA

`Rscript plot_alignments_prediction.r merged_all_LB.f.bed merged_all_DMEM.f.bed LB DMEM_Region1551 3932999 3941000 + + transcripts.bed $REF_GFF_FILE jpeg`

### Classical operons but no transcription -- trp and lac

`Rscript plot_alignments_prediction.r merged_all_LB.r.bed merged_all_DMEM.r.bed LB DMEM_trp 1314000 1324000 - - transcripts.bed $REF_GFF_FILE jpeg`

### Ribosomal protein operon, frequently seen in papers, many predicted transcripts, potential overcalling

`Rscript plot_alignments_prediction.r merged_all_LB.r.bed merged_all_DMEM.r.bed LB DMEM_ribosome_v1 3430000 3460000 - - transcripts.bed $REF_GFF_FILE jpeg`

### Figure 1 - Threonine operon (attenutation)

`Rscript plot_alignments_prediction.r merged_all_LB.f.bed merged_all_DMEM.f.bed LB DMEM_0_6 1 6000 + + transcripts.bed $REF_GFF_FILE jpeg`

### Figure 3 - fdhD/fdoGHI/fdhE - visualizing problem with CDS focused analyses

`Rscript plot_alignments_prediction.r merged_all_LB.r.bed merged_all_DMEM.r.bed LB DMEM_Region_3404 4079000 4088000 - - transcripts.bed $REF_GFF_FILE jpeg`

### Figure 6AB

`samtools view merged_all.bam -b -h -o merged_all_wo_rRNA.bam -U rRNA.bed`

`bamToBed -i merged_all_wo_rRNA.bam > merged_all_wo_rRNA.readtable.bed`

`RScript Figure6AB.R`

## _E. coli_ 2348/69

### Figure 4 - Uses a different version of plotting script for genomes with plasmids or >1 contig, LEE4 (+)- or (-)-strand all reads

`Rscript plot_alignments_prediction_seqid.r merged_all_E2348_69.readtable.bed merged_all_E2348_69.readtable.bed plus_region_3 minus_region_3 72000 83000 + - E2348_69_transcripts.bed.bed $REF_GFF_FILE NZ_CP059840.2 jpeg`

### LEE4 (+)-strand, LB/DMEM

`Rscript plot_alignments_prediction_seqid.r E2348_69_LB.f.bed E2348_69_DMEM.f.bed LEE_LB_plus LEE_DMEM_plus 72000 83000 + + E2348_69_transcripts.bed.bed $REF_GFF_FILE NZ_CP059840.2 jpeg`

### LEE4 (-)-strand, LB/DMEM

`Rscript plot_alignments_prediction_seqid.r E2348_69_LB.r.bed E2348_69_DMEM.r.bed LEE_LB_minus LEE_DMEM_minus 72000 83000 - - E2348_69_transcripts.bed.bed $REF_GFF_FILE NZ_CP059840.2 jpeg`

## _L. monocytogenes_ RO15

### First 30 kbp of genome with + and - strands

`Rscript plot_alignments_prediction.r merged_all_CADEH.readtable.bed merged_all_CADEH.readtable.bed plus_CADEH minus_CADEH 1 30000 + - transcripts_CADEH.bed $REF_GFF_FILE jpeg`

## _P. aeruginosa_ SG17M

### Region of SG17M for verification

`Rscript plot_alignments_prediction.r merged_all_SG17M.readtable.bed merged_all_SG17M.readtable.bed SG17M_plus minus_v2 990000 1020000 + - transcripts_SG17M.bed $REF_GFF_FILE jpeg`

## _P. aeruginosa_ NN2

### Region of NN2 for verification

`Rscript plot_alignments_prediction.r merged_all_NN2.readtable.bed merged_all_NN2.readtable.bed NN2_plus minus_v2 765000 800000 + - transcripts_NN2.bed.bed $REF_GFF_FILE jpeg`

## _Haloferax volcanii_

### Region for verification

`Rscript ../../plot_alignments_prediction_seqid.r merged_all.readtable.bed merged_all.readtable.bed plus minus 1430000 1450000 + - 021324_transcripts.bed ncbi_dataset/data/GCF_000025685.1/genomic.gff NC_013967.1 jpeg`

## Figure 6DE

### Pull fastq from SRA: SRR23886069

`fasterq-dump SRR23886069`

### Combine SINV and ENO2 references for mapping -- TE12 reference from IU, ENO2 from guppy

`cat TE12_Barcode_4.10_Theoretical.fa YHR174W.fasta > SINV_ENO2_combined_ref.fasta`

### Use fasta formatter tool to make sure lines are the same width and format won't interfere with mapping

`fasta_formatter -i SINV_ENO2_combined_ref.fasta -w 80 > SINV_ENO2_combined_ref_formatted.fasta`

### Map to combined, formatted reference with minimap2

`minimap2 -ax map-ont -t 2 SINV_ENO2_combined_ref_formatted.fasta SRR23886069.fastq > SRR23886069_SINV_ENO2.sam`

### Filter the sam file to only include primary mappings

`samtools view -bhF 2308 SRR23886069_SINV_ENO2.sam | samtools sort -o SRR23886069_SINV_ENO2_primary.bam`

### Split the bam file based on the reference

`bamtools split -in SRR23886069_SINV_ENO2_primary.bam -reference`

### Number of reads in each split bam file

`samtools view SRR23886069_SINV_ENO2_primary.REF_TE12_Barcode_4.10_Theoretical.bam | wc -l`

#### Output: 1394

`samtools view SRR23886069_SINV_ENO2_primary.REF_YHR174W.bam | wc -l`

#### Output: 1,589,752

### Pull the read ends from each bam file with custom script

```
foreach line (`ls *.bam`)
Rscript read_end.r $line
mv read_ends.txt $line.read_ends.txt
end
```

### Filter the *read_ends.txt files for only reads that end at 11710 or higher for SINV, or 1310 or higher for ENO2

`awk '{if($1 >= 1310) print}' SRR23886069_SINV_ENO2_primary.REF_YHR174W.bam.read_ends.txt > ENO2_read_end_1310.txt`

`awk '{if($1 >= 11710) print}' SRR23886069_SINV_ENO2_primary.REF_TE12_Barcode_4.10_Theoretical.bam.read_ends.txt > SINV_read_end_11710.txt`

### Number of reads after filtering

`wc -l SINV_read_end_11710.txt`

#### Output: 231

`wc -l ENO2_read_end_1310.txt`

#### output: 1,444,997

### Place only read IDs from filtered read_end.txt files into separate file

`awk '{print $2}' ENO2_read_end_1310.txt > ENO2_read_end_1310_ids.txt`

`awk '{print $2}' SINV_read_end_11710.txt > SINV_read_end_11710_ids.txt`

### Filter split bam files for reads included in the filtered read_end.txt files

`samtools view SRR23886069_SINV_ENO2_primary.REF_TE12_Barcode_4.10_Theoretical.bam | grep -wF -f SINV_read_end_11710_ids.txt > SINV_read_end_11710.sam`

`samtools view SRR23886069_SINV_ENO2_primary.REF_YHR174W.bam | grep -wF -f ENO2_read_end_1310_ids.txt > ENO2_read_end_1310.sam`

### Pull the header from the unfiltered split bam files and add them to the new, filtered bam files

`samtools view -H SRR23886069_SINV_ENO2_primary.REF_TE12_Barcode_4.10_Theoretical.bam > SINV_header.txt`

`samtools view -H SRR23886069_SINV_ENO2_primary.REF_YHR174W.bam > ENO2_header.txt`

`cat SINV_header.txt SINV_read_end_11710.sam | samtools view -bho SINV_read_end_11710.bam`

`cat ENO2_header.txt ENO2_read_end_1310.sam | samtools view -bho ENO2_read_end_1310.bam`

### Convert bam to bed files

`bedtools bamtobed -i ENO2_read_end_1310.bam > ENO2_read_end_1310.bed`

`bedtools bamtobed -i SINV_read_end_11710.bam > SINV_read_end_11710.bed`

### Plot the reads, left-sorted

`Rscript plot_alignments_single_panel.r SINV_read_end_11710.bed SINV_IVT 0 11800 + jpeg`

`Rscript plot_alignments_single_panel.r ENO2_read_end_1310.bed ENO2 0 1320 + jpeg`

## Figure 6C

`samtools depth -d 0 -r NC_000913.3:225728-228692 merged_all.bam > merged_all.23S.sorted.depth`

`samtools depth -d 0 -r NC_000913.3:223740-225342 merged_all.bam > merged_all.16S.sorted.depth`

`samtools depth -d 0 SRR23886069_SINV_ENO2_primary.REF_YHR174W.bam > SINV_IVT_20220512_allpos.txt`

`Rscript Figure6C.R`

## Figure 6F

### Pull FASTQ files for ONT reads from SRA

`fasterq-dump SRR23886071`

### Map reads to reference, which is only transcripts that end with .1 and do not have a,b,c,d, etc. variants

`minimap2 -ax map-ont -t 2 Bm.v4.single.transcripts.fa SRR23886071.fastq > SRR23886071.transcripts.sam`

### Sort, index, and calculate idxstats

`samtools sort SRR23886071.sam -o SRR23886071.bam`

`samtools index SRR23886071.bam`

`samtools idxstats SRR23886071.bam > SRR23886071.transcripts.idxstats`

### Pull FASTQ files for Illumina reads from SRA

`fasterq-dump SRR3111494`

`module load python/3.11.4`

### Map reads to reference, which is only transcripts that end with .1 and do not have a,b,c,d, etc. variants

`python3 hisat2-build -c Bm.v4.single.transcripts.fa Bm.v4.single.transcripts`
`hisat2 -x Bm.v4.single.transcripts.fa -1 SRR3111494_1.fastq -2 SRR3111494_2.fastq`

hisat2 doesn't work though. Switch to bowtie2.

`bowtie2-build Bm.v4.single.transcripts.fa Bm.v4.single.transcripts`
`bowtie2 -x Bm.v4.single.transcripts -1 SRR3111494_1.fastq -2 SRR3111494_2.fastq -S SRR3111494.transcripts.bowtie.sam`

### Sort, index, and calculate idxstats

`samtools sort SRR3111494.transcripts.bowtie.sam -o SRR3111494.bam`

`samtools index SRR3111494.bam`

`samtools idxstats SRR3111494.bam > SRR3111494.transcripts.idxstats`

### Combine the idxstats for the Illumina and ONT reads

`paste SRR3111494.transcripts.idxstats SRR23886071.transcripts.idxstats | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | perl -ne 'chomp;@line=split/\t/;$length=$line[1]/1000;if (($line[1]!=0)&&($line[2]!=0)) {$rpkI=$line[2]/$length;$rpkO=$line[6]/$length;$sumrpkI=$sumrpkI+$rpkI;$sumrpkO=$sumrpkO+$rpkO;$rpkfactorI=$sumrpkI/1000000;$rpkfactorO=$sumrpkO/1000000;print "$line[0]\t$length\t$line[2]\t$line[6]\t$rpkI\t$rpkO\t$sumrpkI\t$sumrpkO\t$rpkfactorI\t$rpkfactorO\n"}' | tail -1`

`paste SRR3111494.transcripts.idxstats SRR23886071.transcripts.idxstats | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | perl -ne 'chomp;@line=split/\t/;$length=$line[1]/1000;if (($line[1]!=0)&&($line[2]!=0)) {$rpkI=$line[2]/$length;$rpkO=$line[6]/$length;$tpmI=$rpkI/31.7342963734812;$tpmO=$rpkO/1.40533802127267;$ratio=$tpmO/$tpmI;if($ratio!=0) {$log2ratio=log($ratio)/log(2);print "$line[0]\t$length\t$line[2]\t$line[6]\t$rpkI\t$rpkO\t$tpmI\t$tpmO\t$log2ratio\n"}}' > Table_for_Figure6F.txt`

`cat > header
Transcript      Length  Illumina        ONT     rpkI    rpkO    tpmIllumina     tpmONT    log2Ratio`

`cat header Table_for_Figure6F.txt > Table_for_Figure6F_header.txt`

`Rscript Figure6F.R`

## Supplementary File 1 - Testing transcript assembly methods (Kaylee's notes)

### FASTQ Files used:

SRR27982843 LB
SRR18070402 DMEM

### Stringtie LB with default parameters

#### K12 LB

`stringtie SRR27982843.bam -o SRR27982843.gtf -p 8`

#### K12 DMEM

`stringtie SRR18070402.bam -o SRR18070402.gtf -p 8`

### Tama run 1 with default parameters:

```
module unload python
module load tama
```

#### K12 LB

`tama_collapse.py -b BAM -s SRR27982843.bam -f GCF_000005845.2_ASM584v2_genomic.fna -p SRR27982843_tama_default -x no_cap -rm low_mem`

#### K12 DMEM

`tama_collapse.py -b BAM -s SRR18070402.bam -f GCF_000005845.2_ASM584v2_genomic.fna -p SRR18070402_tama_default -x no_cap -rm low_mem`

### Tama run 2 with “longest ends” parameter, 100 wobble threshold for both 5’ and 3’ ends (default is 10), coverage requirement 80% (default is 99), and identity calculation method “ident_map” to exclude hard and soft clipping from % identity calculation

#### K12 LB

`tama_collapse.py -b BAM -s SRR27982843.bam -f GCF_000005845.2_ASM584v2_genomic.fna -p SRR27982843_tama_100wobble_c80 -x no_cap -e longest_ends -c 80 -icm ident_map -rm low_mem -a 100 -z 100`

#### K12 DMEM

`tama_collapse.py -b BAM -s SRR18070402.bam -f GCF_000005845.2_ASM584v2_genomic.fna -p SRR18070402_tama_100wobble_c80 -x no_cap -e longest_ends -c 80 -icm ident_map -rm low_mem -a 100 -z 100`

### Tama LB 3, “common ends” parameter, the rest are the same as ‘Tama LB 2’

`tama_collapse.py -b BAM -s SRR27982843.bam -f GCF_000005845.2_ASM584v2_genomic.fna -p SRR27982843_tama_commonends -x no_cap -e common_ends -c 80 -icm ident_map -rm low_mem -a 100 -z 100`

`tama_collapse.py -b BAM -s SRR18070402.bam -f GCF_000005845.2_ASM584v2_genomic.fna -p SRR18070402_tama_commonends -x no_cap -e common_ends -c 80 -icm ident_map -rm low_mem -a 100 -z 100`

### Cupcake LB 1, default (this one was empty in IGV because the default filters don’t work with nanopore data):

`collapse_isoforms_by_sam.py --input SRR18070402.fastq --fq -b SRR18070402.bam -o SRR18070402_cupcake_default --cpus 4`

### Cupake LB 2, 95% coverage and 85% identity requirements:

`collapse_isoforms_by_sam.py -c 0.95 -i 0.85 --input SRR18070402.fastq --fq -b SRR18070402.bam -o SRR18070402_cupcake_c95i85 --cpus 4`
