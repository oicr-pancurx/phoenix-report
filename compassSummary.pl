#!/usr/bin/perl

use strict;
use warnings;
use JSON;
use Tabix;
use Time::Local;

my $pcsi = $ARGV[0];
my $compass = $ARGV[1];

my $analysisPath = "/oicr/data/archive/projects/PCSI/COMPASS/";
my $htmlDir = "$analysisPath/summaryHtml/";
my $biopsyFile = "/oicr/data/archive/projects/PCSI/COMPASS/COMPASS_Biopsy_Times.txt";
my $reportFile = "/oicr/data/archive/projects/PCSI/COMPASS/COMPASS_DNA_Reported.txt";
my $historyFile = "/oicr/data/archive/projects/PCSI/COMPASS/pdac_20160805_matrix.csv";

my $version = "v0.1";
my $currentTime = localtime;

my $alexandrovClass = "";
my $waddellClass = "";


my %geneNotes;

my $l;
my @f;

open (FILE, "$analysisPath/geneNotes.txt") or die "Couldn't open geneNotes.txt\n";
while ($l = <FILE>)
{
	chomp $l;
	@f = split(/\t/, $l);

	if (exists $geneNotes{$f[0]})
	{
		$geneNotes{$f[0]} .= ";" . $f[1];
	}
	else
	{
		$geneNotes{$f[0]} = $f[1];
	}
}
close FILE;

## Metadata
my ($id, $date);
my $biopsyDate = "null";
my $reportDate = "null";

my %limsTissue = (
	"Ad" => "Adipose",
	"Ap" => "Appendix",
	"Ag" => "Adrenal Gland",
	"Bm" => "Bone Marrow",
	"Bn" => "Brain",
	"Br" => "Breast",
	"Bu" => "Buccal cells",
	"Cb" => "Cored Blood",
	"Cn" => "Central Nervous System",
	"Du" => "Duodenum",
	"Ep" => "Esophagus",
	"Es" => "Esophagus",
	"Fs" => "Foreskin",
	"Gb" => "Gallbladder",
	"Hr" => "Heart",
	"Ki" => "Kidney",
	"Le" => "Leukocyte",
	"Li" => "Large Intestine",
	"Ln" => "Lymph Node",
	"Lu" => "Lung",
	"Lv" => "Liver",
	"Lx" => "Larynx",
	"Ly" => "Lymphocyte",
	"Md" => "Mediastinum",
	"Me" => "Mesenchyme",
	"Nk" => "Neck",
	"Oc" => "Oral Cavity",
	"Om" => "Omentum",
	"Ov" => "Ovary",
	"Pa" => "Pancreas",
	"Pb" => "peripheral blood",
	"Pr" => "Prostate",
	"Sa" => "Saliva",
	"Sg" => "Salivary Gland",
	"Si" => "Small Intestine",
	"Sk" => "Skin",
	"Sm" => "Skeletal Muscle",
	"Sp" => "Spleen",
	"St" => "Stomach",
	"Ta" => "Tail (Referring to Models)",
	"Tr" => "Trachea",
	"Mu" => "Muscle",
	"Wm" => "Worm (Referring to Models)",
	"nn" => "Unknown"
);

open (FILE, $biopsyFile) or die "Couldn't open $biopsyFile\n";
while ($l = <FILE>)
{
	chomp $l;
	($id, $date) = split(/\t/, $l);
	if ($id eq $compass)
	{
		$biopsyDate = $date;
	}
}
close FILE;

open (FILE, $reportFile) or die "Couldn't open $reportFile\n";
while ($l = <FILE>)
{
	chomp $l;
	($id, $date) = split(/\t/, $l);
	if ($id eq $compass)
	{
		$reportDate = $date;
	}
}
close FILE;

my ($day,$mon,$year);
$day = 1;
$mon = 1;
$year = 2000;
unless ($biopsyDate eq "null")
{
	($day,$mon,$year) = split(/\//, $biopsyDate);
	$mon--;
}
my $biopsyTickcount = timelocal(0,0,0,$day,$mon,$year);
my $daysElapsed;
unless ($reportDate eq "null")
{
	($day,$mon,$year) = split(/\//, $reportDate);
	$mon--;
	my $reportTickcount = timelocal(0,0,0,$day,$mon,$year);
	
	$daysElapsed = ($reportTickcount - $biopsyTickcount) / (60*60*24);
}
else
{
	$daysElapsed = -1;
}

my $tumTissue = "null";
if ($pcsi =~ /^PCSI_...._(..)_./)
{
	$tumTissue = $limsTissue{$1};
}

print "COMPASS DNA Summary\n";
print "$compass ($pcsi)\n";
print "Biopsied: $biopsyDate\n";
print "Reported: $reportDate\n";
print "Days Elapsed: $daysElapsed\n";
print "Tissue: $tumTissue\n";


### QC metrics
my ($tumour, $normal);
my ($tumCoverage, $norCoverage);
my ($tumError, $norError);
my ($tumSoft, $norSoft);
my $ls;

$tumour = $pcsi;

$ls = `ls $analysisPath/*/$pcsi/wgs/bwa/0.6.2/final_strelka-mutect/$pcsi.final.vcf`;
chomp $ls;
open (FILE, $ls) or die "Couldn't open $ls\n";
while ($l = <FILE>)
{
	chomp $l;
	if ($l =~ /matched_sample_id="(.*?)"/)
	{
		$normal = $1;
		last;
	}
}
close FILE;

my %jsonHash;
$ls = `ls $analysisPath/*/$tumour/wgs/bwa/0.6.2/json/$tumour.json`;
chomp $ls;
open (FILE, $ls) or die "Couldn't open $ls\n";
$l = <FILE>;
$jsonHash{j} = decode_json($l);
$tumCoverage = ($jsonHash{j}{"aligned bases"} * ($jsonHash{j}{"reads on target"} / $jsonHash{j}{"mapped reads"}) ) / $jsonHash{j}{"target size"};
$tumError = (($jsonHash{j}{"mismatch bases"} + $jsonHash{j}{"inserted bases"} + $jsonHash{j}{"deleted bases"}) / $jsonHash{j}{"aligned bases"}) * 100;
$tumSoft = $jsonHash{j}{"soft clip bases"} / ($jsonHash{j}{"aligned bases"} + $jsonHash{j}{"soft clip bases"} + $jsonHash{j}{"hard clip bases"}) * 100;
close FILE;

$ls = `ls $analysisPath/*/$normal/wgs/bwa/0.6.2/json/$normal.json`;
chomp $ls;
open (FILE, $ls) or die "Couldn't open $ls\n";
$l = <FILE>;
$jsonHash{j} = decode_json($l);
$norCoverage = ($jsonHash{j}{"aligned bases"} * ($jsonHash{j}{"reads on target"} / $jsonHash{j}{"mapped reads"}) ) / $jsonHash{j}{"target size"};
$norError = (($jsonHash{j}{"mismatch bases"} + $jsonHash{j}{"inserted bases"} + $jsonHash{j}{"deleted bases"}) / $jsonHash{j}{"aligned bases"}) * 100;
$norSoft = $jsonHash{j}{"soft clip bases"} / ($jsonHash{j}{"aligned bases"} + $jsonHash{j}{"soft clip bases"} + $jsonHash{j}{"hard clip bases"}) * 100;
close FILE;

print "\n";
print "Tumour Coverage: " . sprintf("%.1f",$tumCoverage) . "x\n";
print "Normal Coverage: " . sprintf("%.1f",$norCoverage) . "x\n";

print "Tumour Error Rate: " . sprintf("%.2f",$tumError) . "%\n";
print "Normal Error Rate: " . sprintf("%.2f",$norError) . "%\n";

print "Tumour Soft Clip Rate: " . sprintf("%.2f",$tumSoft) . "%\n";
print "Normal Soft Clip Rate: " . sprintf("%.2f",$norSoft) . "%\n";


# read cellularity and ploidy from celluloid
$ls = `ls $analysisPath/*/$pcsi/wgs/bwa/0.6.2/celluloid/v11/solution/parameters_$pcsi.txt`;
chomp $ls;
open (FILE, $ls) or die "Couldn't open $ls\n";
$l = <FILE>;	# header (value S N T1 Ploidy)
$l = <FILE>;
@f = split(/ /, $l);

my $cellularity = sprintf("%.2f",$f[3] * 100);
my $ploidy = sprintf("%.2f",$f[4]);

print "\n";
print "Cellularity: " . $cellularity . "%\n";
print "Ploidy: " . $ploidy . "\n";

# call patient sex from qc files
my ($xReads, $yReads, $chr);
my @XYs;
my $ratioString = "";
$ls = `ls $analysisPath/*/$normal/wgs/bwa/0.6.2/sampleIdentity/gendertype/*.gendertype`;
chomp $ls;
for my $file (split/\n/, $ls)
{
	open (FILE, $file) or die "Couldn't open $file\n";
	$l = <FILE>;
	chomp $l;
	$l =~ s/^ *//;
	($xReads,$chr) = split(/ /, $l);
	$l = <FILE>;
	warn "$l\n";
	$l =~ s/^ *//;
	($yReads,$chr) = split(/ /, $l);

	close FILE;
	warn "$xReads / $yReads\n";
	
	if (($xReads > 0) and ($yReads > 0))
	{
		push (@XYs, $xReads / $yReads);
		$ratioString .= int($xReads / $yReads) . ",";
	}
}
$ratioString =~ s/,$//;

my $inferredSex = "null";
my $diffRatio = 30;
for my $xy (@XYs)
{
	if ($xy > $diffRatio)
	{
		if (($inferredSex eq "null") or ($inferredSex eq "Female"))
		{
			$inferredSex = "Female";
		}
		else
		{
			$inferredSex = "Unclear";
		}
	}
	else
	{
		if (($inferredSex eq "null") or ($inferredSex eq "Male"))
		{
			$inferredSex = "Male";
		}
		else
		{
			$inferredSex = "Unclear";
		}
	}
}

## read file for population numbers
my %history;
my @historySamps;

my @header;
my $sample;
open (FILE, $historyFile) or die "Couldn't open $historyFile\n";

while ($l = <FILE>)
{
    chomp $l;
    if ($l =~ /^donor/)
    {
        @header = split(/,/, $l);
    }
    else
    {
        @f = split(/,/, $l);
        $sample = $f[1];

        for (my $i = 0; $i < scalar(@f); $i++)
        {
            $history{$sample}{$header[$i]} = $f[$i];
        }
		push (@historySamps, $sample);

    }
}
close FILE;


### Germline with filtering

my $germVars = `grep -v "#" $analysisPath/*/$tumour/wgs/bwa/0.6.2/final_gatk-germline/*.final.vcf | wc -l`;
chomp $germVars;

my $germNonsil = `grep -v "#" $analysisPath/*/$tumour/wgs/bwa/0.6.2/final_gatk-germline/*.final.vcf | grep -wf /u/aconnor/lists/non_silent_variants | wc -l`;
chomp $germNonsil;

my $germNonsilGenes = `grep -v "#" $analysisPath/*/$tumour/wgs/bwa/0.6.2/final_gatk-germline/*.final.vcf | grep -wf /u/aconnor/lists/non_silent_variants | grep -wf /u/aconnor/lists/hpcs_genes | wc -l`;
chomp $germNonsilGenes;

my $germNonsilGenesRare = `cat $analysisPath/*/$tumour/wgs/bwa/0.6.2/final_gatk-germline/annovar/*.genes_nonsil_rare.vcf | wc -l`;
chomp $germNonsilGenesRare;

my $germPathogenic = `cat $analysisPath/*/$tumour/wgs/bwa/0.6.2/final_gatk-germline/*path* | wc -l`;
chomp $germPathogenic;

print "\n";
print "Germline Variants: $germVars\n";
print "Germline Non-silent Variants: $germNonsil\n";
print "Germline Non-silent Variants in HPCS Genes: $germNonsilGenes\n";
print "Germline Rare Non-silent Variants in HPCS Genes: $germNonsilGenesRare\n";
print "Germline Pathogenic Variants: $germPathogenic\n";


### mutational load

my $nonsynCount = 0;
my $stopgainCount = 0;
my $stoplossCount = 0;
my $splicingCount = 0;
my $deleteriousSNVCount = 0;

my $framedelCount = 0;
my $frameinsCount = 0;
my $nondelCount = 0;
my $noninsCount = 0;
my $deleteriousIndelCount = 0;

my $delBPCount = 0;
my $dupBPCount = 0;
my $invBPCount = 0;
my $traBPCount = 0;
my $deleteriousSVCount = 0;

my $snvCount = `grep -v "^#" $analysisPath/*/$pcsi/wgs/bwa/0.6.2/final_strelka-mutect/$pcsi.final.vcf | grep -v TIR | wc -l`;
chomp $snvCount;

my $indelCount = `grep -v "^#" $analysisPath/*/$pcsi/wgs/bwa/0.6.2/final_strelka-mutect/$pcsi.final.vcf | grep TIR | wc -l`;
chomp $indelCount;

my $svCount = `grep -v chrom1 $analysisPath/*/$pcsi/wgs/bwa/0.6.2/final_crest-delly/$pcsi.annotatedSV.tsv | wc -l`;
chomp $svCount;

print "\n";
print "Somatic SNV Count: $snvCount\n";
print "Somatic Indel Count: $indelCount\n";
print "Somatic SV Count: $svCount\n";


### neoantigen load

my $hlaTypes = `cat $analysisPath/*/$pcsi/wgs/bwa/0.6.2/polysolver/1.0/*.txt`;
chomp $hlaTypes;
$hlaTypes =~ s/\n/, /g;

my $neoantigens = `cat $analysisPath/*/$pcsi/wgs/bwa/0.6.2/netMHC/pan-2.8/polysolver/1.0/*.txt | grep Protein`;
chomp $neoantigens;
my $neoCount = 0;
for $l (split(/\n/, $neoantigens))
{
	print $l . "\n";
	if ($l =~ /Number of high binders (.*?)\. Number of weak binders (.*?)\./)
	{
		$neoCount += $1;
		$neoCount += $2;
	}
}

print "\n";
print "Neoantigen Load: $neoCount\n";
print "HLA Types: $hlaTypes\n";


### somatically altered gene table

my ($gene,$consequence,$annovar,$cosmic);
my ($pos,$ref,$alt,$alts,$qual,$filter,$info,$format,$nGT,$tGT);
my %somaticGenes;

$ls = `ls $analysisPath/*/$pcsi/wgs/bwa/0.6.2/celluloid/v11/solution/segments_$pcsi.txt.sorted.bed.gz`;
chomp $ls;
my $tabix = Tabix->new(-data => $ls, -index => "$ls.tbi");
my $iter;
my @tab;

# get variants from sv

my $svFile = `ls $analysisPath/*/$pcsi/wgs/bwa/0.6.2/final_crest-delly/$pcsi.annotatedSV.tsv`;
chomp $svFile;

my %svLong = (
	"DEL" => "deletion-BP",
	"DUP" => "duplication-BP",
	"INV" => "inversion-BP",
	"TRA" => "translocation-BP",
);

open (FILE, $svFile) or die "Couldn't open $svFile\n";
while ($l = <FILE>)
{
	chomp $l;
	unless ($l =~ /^chrom/)
	{
		@f = split(/\t/, $l);

		$consequence = "$svLong{$f[4]} ($f[0]:$f[1]-$f[2]:$f[3])";

		unless ($f[5] eq ".")
		{
			$deleteriousSVCount++;
			for $gene (split/,/, $f[5])
			{
				if ($f[4] eq "DEL")
				{
					$delBPCount++;
				}
				elsif ($f[4] eq "DUP")
				{
					$dupBPCount++;
				}
				elsif ($f[4] eq "INV")
				{
					$invBPCount++;
				}
				else		# TRA
				{
					$traBPCount++;
				}
				$somaticGenes{$gene}{variant}{$consequence}{notes} = "";
				$somaticGenes{$gene}{variant}{$consequence}{freq} = "";

				$somaticGenes{$gene}{cn} = "";
				$somaticGenes{$gene}{af} = "";

				$iter = $tabix->query($f[0],$f[1] - 1,$f[1]);
	
		        if (defined $iter->{"_"})       # happens if the contig isn't in the bed file
		        {
		            while ($l = $tabix->read($iter))		# not handling a gene split over multiple segments!
		            {
		                @tab = split(/\t/, $l);
						$somaticGenes{$gene}{cn} = sprintf("%.2f", $tab[11]);
						if ($f[10] eq "NA")
						{
							$somaticGenes{$gene}{af} = sprintf("%.2f", $tab[9]);
						}
						else
						{
							$somaticGenes{$gene}{af} = sprintf("%.2f", $tab[10]);
						}
		            }
		        }
			}
		}
		unless ($f[6] eq ".")
		{
			for $gene (split/,/, $f[6])
			{
				$somaticGenes{$gene}{variant}{$consequence}{notes} = "";
				$somaticGenes{$gene}{variant}{$consequence}{freq} = "";

				$somaticGenes{$gene}{cn} = "";
				$somaticGenes{$gene}{af} = "";

				$iter = $tabix->query($f[2],$f[3] - 1,$f[3]);
	
		        if (defined $iter->{"_"})       # happens if the contig isn't in the bed file
		        {
		            while ($l = $tabix->read($iter))		# not handling a gene split over multiple segments!
		            {
		                @tab = split(/\t/, $l);
						$somaticGenes{$gene}{cn} = sprintf("%.2f", $tab[11]);
						if ($f[10] eq "NA")
						{
							$somaticGenes{$gene}{af} = sprintf("%.2f", $tab[9]);
						}
						else
						{
							$somaticGenes{$gene}{af} = sprintf("%.2f", $tab[10]);
						}
					}
				}

			}
		}
	}
}



# get variants from ssm

my ($totBases,$altBases,$depth,$freq);



$ls = `ls $analysisPath/*/$pcsi/wgs/bwa/0.6.2/final_strelka-mutect/$pcsi.final.vcf`;
chomp $ls;

open (FILE, $ls) or die "Couldn't open $ls\n";

while ($l = <FILE>)
{
	chomp $l;

	if ($l =~ /^#/)
	{
	}
	else
	{
		($chr,$pos,$id,$ref,$alts,$qual,$filter,$info,$format,$nGT,$tGT) = split(/\t/, $l);
		$gene = "";
		$consequence = "";
		$cosmic = "";
		
		$info .= ";";
		if ($info =~ /ANNOVAR=exonic,(.*?);/)
		{
			$gene = $1;
		}
		if ($info =~ /ANNOVAR=splicing,(.*?)\((.*?)\);/)
		{
			$gene = $1;
			$consequence = "splicing ($2)";
			$consequence =~ s/\(.*:/(/;
		}
		if ($info =~ /ANNOVAR_EXONIC=(.*?);/)
		{
			$consequence = $1;
			$consequence =~ s/,$gene:.*?:.*?:/ (/;
			$consequence =~ s/,.*/)/;
		}
		if ($info =~ /COSMIC=(.*?);/)
		{
			$cosmic = $1;
		}

		$alt = $alts;		# need to split lines with multiple alleles to annotate properly
		$alt =~ s/,.*//;

            if ($format eq "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")     # strelka format
            {
                if ($tGT =~ /^(.*?):.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
                {
                    $totBases = $1;

                    if ($alt eq "A")
                    {
                        $altBases = $2;
                    }
                    elsif ($alt eq "C")
                    {
                        $altBases = $3;
                    }
                    elsif ($alt eq "G")
                    {
                        $altBases = $4;
                    }
                    elsif ($alt eq "T")
                    {
                        $altBases = $5;
                    }

                    $depth = $totBases;
                    $freq = $altBases / $depth;
                }
                else
                {
                    die "Assumed Strelka format, couldn't parse frequencies from $tGT\n";
                }
            }
            elsif ($format eq "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50")
            {
                if ($tGT =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
                {
                    $totBases = $1;
                    $altBases = $2;

                    $depth = $totBases;
                    $freq = $altBases / $depth;
                }
                else
                {
                    die "Assumed Strelka format, couldn't parse frequencies from $tGT\n";
                }
            }
			else
			{
				die "Couldn't determine format of $format\n";
			}

			if ($consequence =~ /splicing/)
			{
				$splicingCount++;
				$deleteriousSNVCount++;
			}
			elsif ($consequence =~ /stopgain/)
			{
				$stopgainCount++;
				$deleteriousSNVCount++;
			}
			elsif ($consequence =~ /stoploss/)
			{
				$stoplossCount++;
				$deleteriousSNVCount++;
			}
			elsif ($consequence =~ /nonsynonymous/)
			{
				$nonsynCount++;
				$deleteriousSNVCount++;
			}
			elsif ($consequence =~ /^frameshift-deletion/)
			{
				$framedelCount++;
				$deleteriousIndelCount++;
			}
			elsif ($consequence =~ /^frameshift-insertion/)
			{
				$frameinsCount++;
				$deleteriousIndelCount++;
			}
			elsif ($consequence =~ /^nonframeshift-deletion/)
			{
				$nondelCount++;
				$deleteriousIndelCount++;
			}
			elsif ($consequence =~ /^nonframeshift-insertion/)
			{
				$noninsCount++;
				$deleteriousIndelCount++;
			}



		if (($consequence =~ /splicing/) or ($consequence =~ /frameshift/) or ($consequence =~ /stop/) or ($consequence =~ /nonsyn/))
		{
			$somaticGenes{$gene}{variant}{$consequence}{notes} = "";
			if ($id ne ".")
			{
				$somaticGenes{$gene}{variant}{$consequence}{notes} .= "$id ";
			}
			if ($cosmic ne "")
			{
				$somaticGenes{$gene}{variant}{$consequence}{notes} .= "$cosmic ";
			}
			if (exists $geneNotes{$gene})
			{
				$somaticGenes{$gene}{variant}{$consequence}{notes} .= $geneNotes{$gene} . " ";
			}

			$somaticGenes{$gene}{variant}{$consequence}{freq} = sprintf("%.0f", $freq * 100) . "% ($altBases/$totBases)";


			$iter = $tabix->query($chr,$pos - 1,$pos);

	        if (defined $iter->{"_"})       # happens if the contig isn't in the bed file
	        {
	            while ($l = $tabix->read($iter))		# not handling a gene split over multiple segments!
	            {
	                @f = split(/\t/, $l);
					$somaticGenes{$gene}{cn} = sprintf("%.2f", $f[11]);
					if ($f[10] eq "NA")
					{
						$somaticGenes{$gene}{af} = sprintf("%.2f", $f[9]);
					}
					else
					{
						$somaticGenes{$gene}{af} = sprintf("%.2f", $f[10]);
					}
	            }
	        }


		}


	}
}
close FILE;

my %genesOfInterest = (
	"KRAS" => "chr12,25358180,25403854",
	"TP53" => "chr17,7571720,7590868",
	"CDKN2A" => "chr9,21967751,21994490",
	"SMAD4" => "chr18,48556583,48611411",
	"ARID1A" => "chr1,27022522,27108601",
	"MAP2K4" => "chr17,11924135,12047051",
	"RNF43" => "chr17,56431038,56494931",
	"TGFBR2" => "chr3,30647994,30735633",
	"KDM6A" => "chrX,44732423,44971845",
);


$consequence = "";


print "\n";
print "Gene\tVariant\tRead Freq\tCopy Number\tMAF\tAdditional Information\n";

for my $gene (qw/KRAS TP53 CDKN2A SMAD4 MAP2K4 ARID1A RNF43 TGFBR2 KDM6A/)
{
	$iter = $tabix->query(split(/,/, $genesOfInterest{$gene}));

	if (defined $iter->{"_"})       # happens if the contig isn't in the bed file
	{
		while ($l = $tabix->read($iter))		# not handling a gene split over multiple segments!
		{
			@f = split(/\t/, $l);
			$somaticGenes{$gene}{cn} = sprintf("%.2f", $f[11]);
			if ($f[10] eq "NA")
			{
				$somaticGenes{$gene}{af} = sprintf("%.2f", $f[9]);
			}
			else
			{
				$somaticGenes{$gene}{af} = sprintf("%.2f", $f[10]);
			}
		}
	}

	if (exists $somaticGenes{$gene}{variant})
	{
		for $consequence (sort keys %{ $somaticGenes{$gene}{variant} })
		{
			print "$gene\t$consequence\t$somaticGenes{$gene}{variant}{$consequence}{freq}\t$somaticGenes{$gene}{cn}\t$somaticGenes{$gene}{af}\t$somaticGenes{$gene}{variant}{$consequence}{notes}\n";
		#	delete $somaticGenes{$gene}{variant}{$consequence};
		}
	}
	else
	{
		print "$gene\t\t\t$somaticGenes{$gene}{cn}\t$somaticGenes{$gene}{af}\t\n";
	}
}

print "\n";
print "Gene\tVariant\tRead Freq\tCopy Number\tMAF\tAdditional Information\n";

for $gene (sort keys %somaticGenes)
{
	if (exists $somaticGenes{$gene}{variant})
	{
		for $consequence (sort keys %{ $somaticGenes{$gene}{variant} })
		{
			print "$gene\t$consequence\t$somaticGenes{$gene}{variant}{$consequence}{freq}\t$somaticGenes{$gene}{cn}\t$somaticGenes{$gene}{af}\t$somaticGenes{$gene}{variant}{$consequence}{notes}\n";
		}
	}
}



## draw history context plots

my %histVals;
my %histCols;

my %histColours = (
	"snv_count" => "#40409699",
	"indel_count" => "#57A3AD99",
	"sv_count" => "#DEA73A99",
	"neo_antigens" => "#D9212099",
	"compass" => "#000000",
);

my %histRange = (
	"snv_count" => "c(0,20000)",
	"indel_count" => "c(0,2000)",
	"sv_count" => "c(0,200)",
	"neo_antigens" => "c(0,100)",
);

for my $type (qw/snv_count indel_count sv_count neo_antigens/)
{
	@{ $histVals{$type} } = valueVector(\%history, \@historySamps, $type);

	for (my $i = 0; $i < scalar(@{ $histVals{$type} }); $i++)
	{
		push(@{ $histCols{$type} }, $histColours{$type});
	}
	push (@{ $histCols{$type} }, $histColours{compass});
}

my $snvPR = getPercentRank($snvCount, \@{ $histVals{snv_count} });
my $indelPR = getPercentRank($indelCount, \@{ $histVals{indel_count} });
my $svPR = getPercentRank($svCount, \@{ $histVals{sv_count} });
my $neoPR = getPercentRank($neoCount, \@{ $histVals{neo_antigens} });



push (@{ $histVals{snv_count} }, $snvCount);
push (@{ $histVals{indel_count} }, $indelCount);
push (@{ $histVals{sv_count} }, $svCount);
push (@{ $histVals{neo_antigens} }, $neoCount);

`mkdir $htmlDir/files/${compass}_${pcsi}`;
my $rfile;
open ($rfile, ">$htmlDir/files/${compass}_${pcsi}/histPlots.R") or die "Couldn't open >$htmlDir/files/${compass}_${pcsi}/histPlots.R\n";

for my $type (qw/snv_count indel_count sv_count neo_antigens/)
{
	printArrayToR($rfile,$type,\@{ $histVals{$type} });
	printQuotedArrayToR($rfile,"${type}_col",\@{ $histCols{$type} });

	print $rfile "png(\"$htmlDir/files/${compass}_${pcsi}/${type}_hist.png\",width=240,height=100)\n";
	print $rfile "par(mar=c(1,3,1,0)+0.1, las=1,cex=0.5)\n";
	print $rfile "bp = barplot(${type}\[sort($type,index.return=T)\$ix], ylim=$histRange{$type}, col=as.vector(${type}_col[sort($type,index.return=T)\$ix]), border=NA, space=0)\n";
	print $rfile "abline(v=bp[match(length($type),sort($type,index.return=T)\$ix)])\n";
	print $rfile "dev.off()\n";

	print $rfile "\n";
}


close $rfile;
`Rscript $htmlDir/files/${compass}_${pcsi}/histPlots.R`;




### Print website

my $htmlFiles = "./files/${compass}_${pcsi}/";
my $htmlCommon = "./files/common/";

`cp $analysisPath/*/$pcsi/wgs/bwa/0.6.2/integration/${pcsi}.onePage.png $htmlDir/files/${compass}_${pcsi}/onePage.png`;
`cp $analysisPath/*/$pcsi/wgs/bwa/0.6.2/celluloid/v11/solution/contour_${pcsi}.png $htmlDir/files/${compass}_${pcsi}/contour.png`;

open (HTML, ">$htmlDir/${compass}_${pcsi}_summary.html") or die "Couldn't open >$htmlDir/${compass}_${pcsi}_summary.html\n";

print HTML "<!DOCTYPE html>
<html>
<head>
<link rel=\"stylesheet\" type=\"text/css\" href=\"$htmlCommon/report_style.css\">
</head>\n";

print HTML "<body>\n";

print HTML "<div id=\"header\">\n";
print HTML "<table border=\"0\" style=\"width:100%\"><tr style=\"vertical-align:top;\"><td style=\"width:33.33%\"><img src=\"$htmlCommon/COMPASS_logo.png\" style=\"height:150px;\"></td><td style=\"width:33.33%\"><h4 style=\"text-align:center;\">$compass</h4></td><td style=\"width:33.33%\"><h4 style=\"text-align:right;\">$currentTime</h4></td></tr></table>\n";
print HTML "</div>\n";

print HTML "<div id=\"main\">\n";
print HTML "<h1>COMPASS Genomics Summary</h1>\n";

print HTML "<table style=\"width:100%\">\n";
print HTML "<tr>\n";
print HTML "<td><b>Study ID:</b> COMPASS Trial</td>";
print HTML "<td><b>Biopsy Date:</b> $biopsyDate</td>";
print HTML "<td><b>Biopsy Tissue:</b> $tumTissue</td>";
print HTML "</tr>\n<tr>\n";
print HTML "<td><b>Patient ID:</b> $compass</td>";
print HTML "<td><b>Report Date:</b> $reportDate</td>";
print HTML "<td></td>";
print HTML "</tr>\n";
print HTML "</table>\n";

print HTML "<p>This report is for research purposes only.</p>\n";


print HTML "<h2>DNA Sequencing Quality Metrics</h2>\n";

print HTML "<table style=\"width:100%\">\n";
print HTML "<tr>\n";
print HTML "<td><b>Tumour Coverage:</b> " . sprintf("%.1f",$tumCoverage) . "x (>45x)</td>";
print HTML "<td><b>Error Rate:</b> " . sprintf("%.2f",$tumError) . "% (<1%)</td>";
print HTML "<td><b>Soft Clip Rate:</b> " . sprintf("%.2f",$tumSoft) . "% (<1%)</td>";
print HTML "</tr>\n<tr>\n";
print HTML "<td><b>Normal Coverage:</b> " . sprintf("%.1f",$norCoverage) . "x (>30x)</td>";
print HTML "<td><b>Error Rate:</b> " . sprintf("%.2f",$norError) . "% (<1%)</td>";
print HTML "<td><b>Soft Clip Rate:</b> " . sprintf("%.2f",$norSoft) . "% (<1%)</td>";
print HTML "</tr>\n";
print HTML "<tr>\n";
print HTML "<td><b>Cellularity:</b> $cellularity%</td>";
print HTML "<td><b>Inferred Sex:</b> $inferredSex ($ratioString)</td>";
print HTML "<td></td>";
print HTML "</tr>\n";
print HTML "</table>\n";





print HTML "<h2>Germline Variants</h2>\n";

print HTML "<p><b>Germline Variants:</b> $germVars<br>\n";
print HTML "<b>Germline Non-silent Variants:</b> $germNonsil<br>\n";
print HTML "<b>Germline Non-silent Variants in HPCS Genes:</b> $germNonsilGenes<br>\n";
print HTML "<b>Germline Rare Non-silent Variants in HPCS Genes:</b> $germNonsilGenesRare<br>\n";
print HTML "<b>Germline Pathogenic Variants:</b> $germPathogenic<br>\n";

print HTML "<br>*HPCS genes: POLE, POLD1, EPCAM, MLH1, MSH2, MSH6, PMS2, BRCA1, BRCA2, PALB2, ATM, APC, MUTYH, CDKN2A, STK11, PRSS1, PRSS2, TP53, SMAD4</p>\n";




print HTML "<h2>Somatic Mutations</h2>\n";

print HTML "<h3>\tMutational Load</h3>\n";

print HTML "<table style=\"width:100%\">\n";
print HTML "<tr>\n";
print HTML "<td style=\"text-align:center\"><b>$snvCount SNVs ($snvPR%PR)</b></td>";
print HTML "<td style=\"text-align:center\"><b>$indelCount Indels ($indelPR%PR)</b></td>";
print HTML "<td style=\"text-align:center\"><b>$svCount SVs ($svPR%PR)</b></td>";
print HTML "<td style=\"text-align:center\"><b>$neoCount NAs ($neoPR%PR)</b></td>";
print HTML "</tr>\n";
print HTML "<tr>\n";
print HTML "<td style=\"text-align:center\"><img src=\"$htmlFiles/snv_count_hist.png\" stype=\"width:100%\"></td>";
print HTML "<td style=\"text-align:center\"><img src=\"$htmlFiles/indel_count_hist.png\" stype=\"width:100%\"></td>";
print HTML "<td style=\"text-align:center\"><img src=\"$htmlFiles/sv_count_hist.png\" stype=\"width:100%\"></td>";
print HTML "<td style=\"text-align:center\"><img src=\"$htmlFiles/neo_antigens_hist.png\" stype=\"width:100%\"></td>";
print HTML "</tr>\n";
print HTML "</table>\n";

print HTML "<h3>Driver Gene Status</h3>\n";

print HTML "<table style=\"width:100%\">\n";
print HTML "<tr>\n";
print HTML "<th>Gene</th><th>Variant</th><th>Read Freq</th><th>Copy Number</th><th>MAF</th><th>Additional Information</th>";
print HTML "</tr>\n";

for my $gene (qw/KRAS TP53 CDKN2A SMAD4 MAP2K4 ARID1A RNF43 TGFBR2 KDM6A/)
{
	if (exists $somaticGenes{$gene}{variant})
	{
		for $consequence (sort keys %{ $somaticGenes{$gene}{variant} })
		{
			print HTML "<tr>";
			print HTML "<td>$gene</td>";
			if ($consequence =~ /^frameshift/)
			{
				print HTML "<td class=\"frameshift\">$consequence</td>";
			}
			elsif ($consequence =~ /^stopgain/)
			{
				print HTML "<td class=\"stopgain\">$consequence</td>";
			}
			else
			{
				print HTML "<td>$consequence</td>";
			}
			print HTML "<td>$somaticGenes{$gene}{variant}{$consequence}{freq}</td>";
			if ($somaticGenes{$gene}{cn} < 0.5)
			{
				print HTML "<td class=\"cn_loss\">$somaticGenes{$gene}{cn}</td>";
			}
			elsif ($somaticGenes{$gene}{cn} < 1.5)
			{
				print HTML "<td class=\"cn_loh\">$somaticGenes{$gene}{cn}</td>";
			}
			elsif ($somaticGenes{$gene}{cn} > ($ploidy*4))
			{
				print HTML "<td class=\"cn_very_gain\">$somaticGenes{$gene}{cn}</td>";
			}
			elsif ($somaticGenes{$gene}{cn} > ($ploidy + 1))
			{
				print HTML "<td class=\"cn_gain\">$somaticGenes{$gene}{cn}</td>";
			}
			else
			{
				print HTML "<td>$somaticGenes{$gene}{cn}</td>";
			}
			print HTML "<td>$somaticGenes{$gene}{af}</td>";
			print HTML "<td>$somaticGenes{$gene}{variant}{$consequence}{notes}</td>";
			print HTML "</tr>\n";
		}
	}
	else
	{
		print HTML "<tr>";
		print HTML "<td>$gene</td>";
		print HTML "<td></td>";
		print HTML "<td></td>";
		if ($somaticGenes{$gene}{cn} < 0.5)
		{
			print HTML "<td class=\"cn_loss\">$somaticGenes{$gene}{cn}</td>";
		}
		elsif ($somaticGenes{$gene}{cn} < 1.5)
		{
			print HTML "<td class=\"cn_loh\">$somaticGenes{$gene}{cn}</td>";
		}
		elsif ($somaticGenes{$gene}{cn} > ($ploidy*4))
		{
			print HTML "<td class=\"cn_very_gain\">$somaticGenes{$gene}{cn}</td>";
		}
		elsif ($somaticGenes{$gene}{cn} > ($ploidy + 1))
		{
			print HTML "<td class=\"cn_gain\">$somaticGenes{$gene}{cn}</td>";
		}
		else
		{
			print HTML "<td>$somaticGenes{$gene}{cn}</td>";
		}
		print HTML "<td>$somaticGenes{$gene}{af}</td>";
		print HTML "<td></td>";
		print HTML "</tr>\n";
	}

}
print HTML "</table>\n";

print HTML "<h3>Structural Variant Classifications</h3>\n";

print HTML "<table style=\"width:100%\">\n";
print HTML "<tr>\n";
print HTML "<td><b>Mutation Class:</b> $alexandrovClass</td>\n";
print HTML "<td><b>Waddell Class:</b> $waddellClass</td>\n";
print HTML "</tr>\n";
print HTML "</table>\n";


print HTML "<h3>Whole Genome Summary</h3>\n";
print HTML "<img src=\"$htmlFiles/onePage.png\" style=\"width:100%\">\n";

print HTML "<h3>Ploidy and Chromothripsis</h3>\n";
print HTML "<img src=\"$htmlFiles/contour.png\" style=\"width:100%\">\n";


print HTML "<h3>Full Deleterious Somatic Variant List</h3>\n";

print HTML "<table style=\"width:100%\">\n";
print HTML "<tr>\n";
print HTML "<td><b>Nonsynonymous-SNVs:</b> $nonsynCount</td>\n";
print HTML "<td><b>Frameshift-deletions:</b> $framedelCount</td>\n";
print HTML "<td><b>Deletion-BPs:</b> $delBPCount</td>\n";
print HTML "</tr>\n";
print HTML "<tr>\n";
print HTML "<td><b>Stopgain-SNVs:</b> $stopgainCount</td>\n";
print HTML "<td><b>Nonframeshift-deletions:</b> $nondelCount</td>\n";
print HTML "<td><b>Duplication-BPs:</b> $dupBPCount</td>\n";
print HTML "</tr>\n";
print HTML "<tr>\n";
print HTML "<td><b>Stoploss-SNVs:</b> $stoplossCount</td>\n";
print HTML "<td><b>Frameshift-Insertions:</b> $frameinsCount</td>\n";
print HTML "<td><b>Inversion-BPs</b> $invBPCount</td>\n";
print HTML "</tr>\n";
print HTML "<tr>\n";
print HTML "<td><b>Splicing-SNVs:</b> $splicingCount</td>\n";
print HTML "<td><b>Nonframeshift-Insertions:</b> $noninsCount</td>\n";
print HTML "<td><b>Translocation-BPs:</b> $traBPCount</td>\n";
print HTML "</tr>\n";
print HTML "<tr>\n";
print HTML "<td><b>Total Deleterious SNVs:</b> $deleteriousSNVCount</td>\n";
print HTML "<td><b>Total Deleterious Indels:</b> $deleteriousIndelCount</td>\n";
print HTML "<td><b>Total Deleterious SV Breakpoints:</b> $deleteriousSVCount</td>\n";
print HTML "</tr>\n";
print HTML "</table>\n";

print HTML "<br>\n";

print HTML "<table style=\"width:100%\">\n";
print HTML "<tr>\n";
print HTML "<th>Gene</th><th>Variant</th><th>Read Freq</th><th>Copy Number</th><th>MAF</th><th>Additional Information</th>";
print HTML "</tr>\n";

for my $gene (sort keys %somaticGenes)
{
	if (exists $somaticGenes{$gene}{variant})
	{
		for $consequence (sort keys %{ $somaticGenes{$gene}{variant} })
		{
			print HTML "<tr>";
			print HTML "<td>$gene</td>";
			if ($consequence =~ /^frameshift/)
			{
				print HTML "<td class=\"frameshift\">$consequence</td>";
			}
			elsif ($consequence =~ /^stopgain/)
			{
				print HTML "<td class=\"stopgain\">$consequence</td>";
			}
			else
			{
				print HTML "<td>$consequence</td>";
			}
			print HTML "<td>$somaticGenes{$gene}{variant}{$consequence}{freq}</td>";
			if ($somaticGenes{$gene}{cn} eq "")
			{
				print HTML "<td>$somaticGenes{$gene}{cn}</td>";
			}
			elsif ($somaticGenes{$gene}{cn} < 0.5)
			{
				print HTML "<td class=\"cn_loss\">$somaticGenes{$gene}{cn}</td>";
			}
			elsif ($somaticGenes{$gene}{cn} < 1.5)
			{
				print HTML "<td class=\"cn_loh\">$somaticGenes{$gene}{cn}</td>";
			}
			elsif ($somaticGenes{$gene}{cn} > ($ploidy * 4))
			{
				print HTML "<td class=\"cn_very_gain\">$somaticGenes{$gene}{cn}</td>";
			}
			elsif ($somaticGenes{$gene}{cn} > ($ploidy + 1))
			{
				print HTML "<td class=\"cn_gain\">$somaticGenes{$gene}{cn}</td>";
			}
			else
			{
				print HTML "<td>$somaticGenes{$gene}{cn}</td>";
			}
			print HTML "<td>$somaticGenes{$gene}{af}</td>";
			print HTML "<td>$somaticGenes{$gene}{variant}{$consequence}{notes}</td>";
			print HTML "</tr>\n";
		}
	}
}
print HTML "</table>";



print HTML "</div>\n";

print HTML "<div id=\"footer\">\n";
print HTML "<table border=\"0\" style=\"width:100%\"><tr style=\"vertical-align:bottom\"><td style=\"width:25%\"><h4>COMPASS Summary $version<br>For Research Purposes Only</h4></td>";
print HTML "<td style=\"width:25%;text-align:center;vertical-align:middle\"><img src=\"$htmlCommon/OICR_logo.jpg\"style=\"height:125px;\"></td>";
print HTML "<td style=\"width:25%;text-align:center;vertical-align:middle\"><img src=\"$htmlCommon/UHN_logo.png\" style=\"height:50px;\"></td>";
print HTML "<td style=\"width:100%\"></td></tr></table>\n";
print HTML "</div>\n";

print HTML "</body>\n</html>\n";








sub valueVector
{
    my $dataRef = shift;
    my $sampRef = shift;
    my $type = shift;

    my @vector;

    for my $samp (@{ $sampRef })
    {
        if (exists $dataRef->{$samp}{$type})
        {
            unless ($dataRef->{$samp}{$type} eq "NA")
            {
                push(@vector,$dataRef->{$samp}{$type});
            }
            else
            {
				#push(@vector, "NA");
            }
        }
        else
        {
			#push(@vector, "NA");
        }
    }

    return @vector;
}



sub printArrayToR
{
    my $fh = shift;
    my $var = shift;
    my $array = shift;

    print $fh "$var <- c($array->[0]";
    for (my $i = 1; $i < scalar(@{ $array }); $i++)
    {
        print $fh ",$array->[$i]";
    }
    print $fh ")\n";
}

sub printQuotedArrayToR
{
    my $fh = shift;
    my $var = shift;
    my $array = shift;

    print $fh "$var <- c(\"$array->[0]\"";
    for (my $i = 1; $i < scalar(@{ $array }); $i++)
    {
        print $fh ",\"$array->[$i]\"";
    }
    print $fh ")\n";
}


sub getPercentRank
{
	my $value = shift;
	my $arrayRef = shift;

	my $rank = "";
	my $lowerCount = 0;
	my $equalCount = 0;

	for (my $i = 0; $i < scalar(@$arrayRef); $i++)
	{
		if ($value > $arrayRef->[$i])
		{
			$lowerCount++;
		}
		elsif ($value == $arrayRef->[$i])
		{
			$equalCount++;
		}
	}

	$rank = ($lowerCount + (0.5 * $equalCount)) / (scalar(@$arrayRef) + 1);

	print $rank . "\n";
	return sprintf("%.0f", ($rank * 100));
}



