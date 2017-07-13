#!/usr/bin/perl

use strict;
use warnings;

use Bio::DB::Fasta;
my %fastaHandles = ("path" => "/oicr/data/reference/genomes/homo_sapiens_mc/UCSC/hg19_random/Genomic/references/fasta/");
my %referenceHash = ("chunkSize" => 10000);


my $dataFile = shift;
my $outName = shift;


my %data;
open (FILE, $dataFile) or die "Couldn't open $dataFile\n";
my $l;
my (@header, @f);
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
        for (my $i = 0; $i < scalar(@f); $i++)
        {
            $data{$header[$i]} = $f[$i];
        }
    }
}
close FILE;



my $sample = $data{tumour};

my %dataPath = (
	"somatic" => $data{ssm_file},
	"germline" => $data{sgv_file},
	"cn_states" => $data{seg_file},
	"structural" => $data{sv_file}
);

my $outputFile = "$outName.onePage.png";





my %colours = (
	"C>A" => "#332288",
	"C>G" => "#88CCEE",
	"C>T" => "#117733",	# transversion
	"T>A" => "#DDCC77",
	"T>C" => "#CC6677",	# transversion
	"T>G" => "#AA4499",
	"ins" => "#44AA99",
	"del" => "#882255",
);

my %fadeColours = (
	"C>A" => "#33228866",
	"C>G" => "#88CCEE66",
	"C>T" => "#11773366",	# transversion
	"T>A" => "#DDCC7766",
	"T>C" => "#CC667766",	# transversion
	"T>G" => "#AA449966",
	"ins" => "#44AA9966",
	"del" => "#88225566",
);

my %cnColours = (
	"1" => "#3A89C9",
	"2" => "#99C7EC",
	"3" => "#E6F5FE",
	"4" => "#FFE3AA",
	"5" => "#F5A275",
	"6" => "#D24D3E"
);

my %bpColours = (
	"DEL" => "#4477AA",
	"DUP" => "#117733",
	"INV" => "#DDCC77",
	"TRA" => "#CC6677",
#	"ITX" => "#AA4499",
);


my %genesOfInterest = ("AKT3",1,"PARP1",1,"NRAS",1,"NOTCH2",1,"NOTCH2NL",1,"PIK3C2B",1,"PTCH2",1,"PIK3CD",1,"PIK3R3",1,"DPYD",1,"GSTM1",1,"EPHA2",1,"MTHFR",1,"PLK3",1,"CDC7",1,"WNT2B",1,"WNT3A",1,"WNT4",1,"WNT9A",1,"ARID1A",1,"MCL1",1,"PTEN",1,"RET",1,"MGMT",1,"FGFR2",1,"CYP2C19",1,"CYP2C9",1,"CYP2C8",1,"CYP2E1",1,"ABCC2",1,"WNT8B",1,"TET1",1,"CYP17A1",1,"CHEK1",1,"HRAS",1,"PIK3C2A",1,"CCND1",1,"ATM",1,"WEE1",1,"RRM1",1,"GSTP1",1,"WNT11",1,"CBL",1,"PGR",1,"SLC22A6",1,"SLCO2B1",1,"WT1",1,"KRAS",1,"ErbB3",1,"MDM2",1,"CDK4",1,"PIK3C2G",1,"WNT1",1,"PTPN11",1,"PXN",1,"GLI1",1,"WNT10B",1,"WNT5B",1,"SLCO1B1",1,"SLCO1B3",1,"SOCS2",1,"BRCA2",1,"FLT1",1,"FLT3",1,"RB1",1,"ERCC5",1,"AKT1",1,"PARP2",1,"ESR2",1,"TEP1",1,"PGF",1,"IGF1R",1,"MAP2K1",1,"IDH2",1,"CYP1A1",1,"DLL4",1,"CYP1A2",1,"CDH1",1,"TSC2",1,"PLK1",1,"TUBB3",1,"ERCC4",1,"SULT1A1",1,"SOCS1",1,"SH2B",1,"BRCA1",1,"ErbB2",1,"TP53",1,"RARA",1,"TOP2A",1,"AURKB",1,"PIK3R5",1,"NF1",1,"WNT3",1,"WNT9B",1,"SOCS3",1,"STAT3",1,"BCL2",1,"DCC",1,"PIK3C3",1,"TYMS",1,"SMAD4",1,"ROCK1",1,"AKT2",1,"MAP2K2",1,"NOTCH3",1,"XRCC1",1,"AURKC",1,"ERCC1",1,"PIK3R2",1,"ERCC2",1,"STK11",1,"CYP2A6",1,"CYP2B6",1,"AXL",1,"ALK",1,"MSH2",1,"ErbB4",1,"IDH1",1,"UGT1A1",1,"GLI2",1,"ERCC3",1,"WNT10A",1,"WNT6",1,"SRC",1,"AURKA",1,"TOP1",1,"GART",1,"BCR",1,"CHEK2",1,"CYP2D6",1,"GSTT1",1,"NF2",1,"EWSR1",1,"WNT7B",1,"PIK3CA",1,"MLH1",1,"RAF1",1,"CTNNB1",1,"VHL",1,"PIK3CB",1,"FANCD2",1,"TOP2B",1,"ATR",1,"MST1R",1,"TERC",1,"RASSF1",1,"WNT5A",1,"WNT7A",1,"SLC15A2",1,"KIT",1,"KDR",1,"PDGFRA",1,"FGFR3",1,"EIF4E",1,"PLK4",1,"ABCG2",1,"UGT2B15",1,"UGT2B17",1,"UGT2B7",1,"TET2",1,"PDGFRB",1,"DHFR",1,"FLT4",1,"PIK3R1",1,"TERT",1,"PLK2",1,"WNT8A",1,"NPM1",1,"ESR1",1,"NOTCH4",1,"CCND3",1,"TUBB",1,"SLC29A1",1,"TPMT",1,"SLC22A1",1,"SLC22A2",1,"EGFR",1,"BRAF",1,"MET",1,"SHH",1,"SMO",1,"PIK3CG",1,"CDK5",1,"CYP3A4",1,"EPHB4",1,"GLI3",1,"CYP3A5",1,"WNT16",1,"WNT2",1,"ABCB1",1,"FGFR1",1,"PTK2",1,"PTK2B",1,"ANGPT1",1,"ANGPT2",1,"NAT1",1,"NAT2",1,"ABL1",1,"JAK2",1,"CDKN2A",1,"NOTCH1",1,"PTCH1",1,"TSC1",1,"ROR2",1,"AR",1,"ARAF",1);

my %chrTopOffset = (
	"chr1" => 0,
	"chr2" => 249250621,
	"chr3" => 492449994,
	"chr4" => 690472424,
	"chr5" => 881626700,
	"chr6" => 1062541960,
	"chr7" => 1233657027,
	"chr8" => 1392795690,
	"end" => 1539159712,
);
my %chrBotOffset =(
	"chr9" => 0,
	"chr10" => 141213431,
	"chr11" => 276748178,
	"chr12" => 411754694,
	"chr13" => 545606589,
	"chr14" => 660776467,
	"chr15" => 768126007,
	"chr16" => 870657399,
	"chr17" => 961012152,
	"chr18" => 1042207362,
	"chr19" => 1120284610,
	"chr20" => 1179413593,
	"chr21" => 1242439113,
	"chr22" => 1290569008,
	"chrX" => 1341873574,
	"chrY" => 1497144134,
	"end" => 1556517700,
);

my $lastL;

my %count = (
	"snv" => 0,
	"ins" => 0,
	"del" => 0,
	"nonsyn" => 0,
	"frameins" => 0,
	"framedel" => 0,
	"tv" => 0,
	"kataegis" => 0
);
my %vcfHash;
my @vcfOrder;

my $tSample;
my $nSample;

my $tCoverage;
my $nCoverage;

my $tCol;
my $nCol;

my $date = `date`;
chomp $date;

my ($chr, $pos, $id, $ref, $alt, $freq, $depth, $consequence, $gene);
my ($refBases, $altBases, $totBases);
my $info;
my $type;
my $change;
my $context;

my $lastChrSNV = "null";
my $lastPosSNV = "null";
my $prevDist;

my %comp = (
	"A>C" => "T>G",
	"A>G" => "T>C",
	"A>T" => "T>A",
	"G>A" => "C>T",
	"G>C" => "C>G",
	"G>T" => "C>A",
);

my %contextComp = (
    "AA" => "TT",
    "AC" => "GT",
    "AG" => "CT",
    "AT" => "AT",
    "CA" => "TG",
    "CC" => "GG",
    "CG" => "CG",
    "CT" => "AG",
    "GA" => "TC",
    "GC" => "GC",
    "GG" => "CC",
    "GT" => "AC",
    "TA" => "TA",
    "TC" => "GA",
    "TG" => "CA",
    "TT" => "AA",
);

# make gene of interest ticks while reading somatic vcf
my (@genePosTop, @genePosBot, @geneNameTop, @geneNameBot);

warn "Opening $dataPath{somatic}\n";

open (FILE, "<$dataPath{somatic}") or die "Couldn't open $dataPath{somatic}\n";
while ($l = <FILE>)
{
	if ($l =~ /##DCC=<analyzed_sample_id="(.*?)">/)
	{
		$tSample = $1;
	}
	elsif ($l =~ /##DCC=<matched_sample_id="(.*?)">/)
	{
		$nSample = $1;
	}
	elsif ($l =~ /##DCC=<analyzed_seq_coverage=(.*?)>/)
	{
		$tCoverage = sprintf("%0.2f", $1);
	}
	elsif ($l =~ /##DCC=<matched_seq_coverage=(.*?)>/)
	{
		$nCoverage = sprintf("%0.2f", $1);
	}
	elsif ($l =~ /^#CHROM/)
	{
		chomp $l;
		@f = split(/\t/, $l);
		if ($f[9] eq $tSample)
		{
			$tCol = 9;
			$nCol = 10;
		}
		elsif ($f[10] eq $tSample)
		{
			$tCol = 10;
			$nCol = 9;
		}
		elsif ($f[10] eq "TUMOR")
		{
			$tCol = 10;
			$nCol = 9;
		}
		else
		{
			die "Couldn't detect genotype column!\n";
		}
	}

	unless ($l =~ /^#/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		$chr = $f[0];
		$pos = $f[1];
		$id = $f[2];
		$ref = $f[3];
		$info = $f[7];
		for $alt (split(/,/, $f[4]))
		{
			push(@vcfOrder, "$chr:$pos");
			if (length($ref) == length($alt))
			{
				$type = "snv";
			
				$change = "$ref>$alt";
				if (exists $comp{$change})
				{
					$change = $comp{$change};
				}

				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{change} = $change;

				$count{snv}++;



				if (($change eq "C>T") or ($change eq "T>C"))
				{
					$count{ts}++;
				}

				if ($chr eq $lastChrSNV)
				{
					$prevDist = $pos - $lastPosSNV;		# assuming sorted input!
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{prevdist} = $prevDist;
				}


				$lastChrSNV = $chr;
				$lastPosSNV = $pos;

			}
			elsif (length($ref) > length($alt))
			{
				$type = "del";

				$count{del}++;

				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{size} = length($ref) - length($alt);
			}
			else
			{
				$type = "ins";

				$count{ins}++;

				$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{size} = length($alt) - length($ref);
			}

			$gene = "";
			$consequence = "";

			for $info (split(/;/, $f[7]))
			{
				if ($info =~ /COSMIC=(.*)/)
				{
					if ($id eq ".")
					{
						$id = $1;
					}
					else
					{
						$id .= ",$1";
					}
				}
				elsif ($info =~ /ANNOVAR=exonic,(.*)/)
				{
					$gene = $1;
				}
				elsif ($info =~ /ANNOVAR=splicing,(.*)\((.*)\)/)
				{
					$gene = $1;
					$consequence =~ "splicing,$1($2)";
				}
				elsif ($info =~ /ANNOVAR_EXONIC=(.*)/)
				{
					$consequence = $1;
					$consequence =~ s/$gene/ $gene/g;
				}
			}

			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{id} = $id;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{gene} = $gene;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} = $consequence;

			if ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^frameshift/)
			{
				$count{"frame$type"}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "fs:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "fs:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^nonsynonymous/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "ns:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "ns:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^stopgain/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "sg:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "sg:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^stoploss/)
			{
				$count{nonsyn}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "sl:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "sl:$gene");
					}
				}
			}
			elsif ($vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{consequence} =~ /^splicing/)
			{
				$count{"splice$type"}++;
				if (exists $genesOfInterest{$gene})
				{
					$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{print} = "yes";

					if (exists $chrTopOffset{$chr})
					{
						push(@genePosTop, $pos + $chrTopOffset{$chr});
						push(@geneNameTop, "sp:$gene");
					}
					elsif (exists $chrBotOffset{$chr})
					{
						push(@genePosBot, $pos + $chrBotOffset{$chr});
						push(@geneNameBot, "sp:$gene");
					}
				}
			}

			$freq = "";
			$depth = "";

	        if ($f[8] eq "GT:AD:DP:GQ:PL")      # GATK format
	        {
	            if ($f[$tCol] =~ /^.*?:(.*?),(.*?):(.*?):.*/)
	            {
	                $refBases = $1;
	                $altBases = $2;
	                $totBases = $3;
	
					$depth = $totBases;
					$freq = $altBases / $depth;
	            }
	            else
	            {
	                die "Assumed GATK format, couldn't parse frequencies from $f[$tCol]\n";
	            }
	        }
	        elsif ($f[8] eq "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")     # strelka format
	        {
	            if ($f[$tCol] =~ /^(.*?):.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
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
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$tCol]\n";
	            }
	        }
			elsif ($f[8] eq "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50")
			{
				if ($f[$tCol] =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
				{
					$totBases = $1;
					$altBases = $2;
	
					$depth = $totBases;
					$freq = $altBases / $depth;
				}
	            else
	            {
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$tCol]\n";
	            }
			}
	
	
			unless ($freq eq "")
			{
				$freq = sprintf("%.4f", $freq);
			}

			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{tfreq} = $freq;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{tdepth} = $depth;



			$freq = "";
			$depth = "";

	        if ($f[8] eq "GT:AD:DP:GQ:PL")      # GATK format
	        {
	            if ($f[$nCol] =~ /^.*?:(.*?),(.*?):(.*?):.*/)
	            {
	                $refBases = $1;
	                $altBases = $2;
	                $totBases = $3;
	
					$depth = $totBases;
					$freq = $altBases / $depth;
	            }
	            else
	            {
	                die "Assumed GATK format, couldn't parse frequencies from $f[$nCol]\n";
	            }
	        }
	        elsif ($f[8] eq "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")     # strelka format
	        {
	            if ($f[$nCol] =~ /^(.*?):.*:(.*?),.*?:(.*?),.*?:(.*?),.*?:(.*?),.*?$/)
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
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$nCol]\n";
	            }
	        }
			elsif ($f[8] eq "DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50")
			{
				if ($f[$nCol] =~ /^(.*?):.*?:.*?,.*?:(.*?),.*/)
				{
					$totBases = $1;
					$altBases = $2;
	
					$depth = $totBases;
					$freq = 0;
					if ($depth > 0)
					{
						$freq = $altBases / $depth;
					}
				}
	            else
	            {
	                die "Assumed Strelka format, couldn't parse frequencies from $f[$nCol]\n";
	            }
			}
	
	
			unless ($freq eq "")
			{
				$freq = sprintf("%.4f", $freq);
			}

			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{nfreq} = $freq;
			$vcfHash{$chr}{$pos}{$type}{"$ref>$alt"}{ndepth} = $depth;


			freeMemChunks($chr, $pos, \%referenceHash);
		}
	}
}
close FILE;


my $titvRatio = "NA";
if ($count{tv} > 0)
{
	$titvRatio = sprintf("%0.4f", ($count{snv} - $count{tv}) / $count{tv});
}







# grab het freqs from GATK

my @hetFreqsTop;
my @hetFreqsBot;

my %hetCountsTop;
my %hetCountsBot;

my $imageWindow = 786280;
my $imageY = 85;

my $germN;
my $germT;

my $freq1;
my $freq2;

warn "Opening $dataPath{germline}\n";
open (FILE, "<$dataPath{germline}") or die "Couldn't open $dataPath{germline}\n";

while ($l = <FILE>)
{
	if ($l =~ /##DCC=<analyzed_sample_id="(.*)">/)
	{
		$germT = $1;
	}
	elsif ($l =~ /##DCC=<matched_sample_id="(.*)">/)
	{
		$germN = $1;
	}
	elsif ($l =~ /#CHROM/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		if ($f[9] eq $germT)
		{
			$germT = 9;
			$germN = 10;
		}
		elsif ($f[9] eq $germN)
		{
			$germN = 9;
			$germT = 10;
		}
		else
		{
			die "Couldn't determine tumour/normal genotype columns in $dataPath{germline}\n";
		}
	}
	elsif ($l =~ /^#/)
	{
	}
	else
	{
		chomp $l;

		@f = split(/\t/, $l);

		if (length($f[3]) == length($f[4]))
		{
			if ($f[$germN] =~ /^0\/1:/)
			{
				unless ($f[$germT] =~ /^\.\/\.:/)
				{
					if ($f[$germT] =~ m/^.*?:(.*?),(.*?):(.*?):/)
					{
						$freq1 = $1/$3;
						$freq2 = $2/$3;

						$chr = $f[0];
						$pos = $f[1];

						if (exists $chrTopOffset{$chr})
						{
							$hetCountsTop{int(($chrTopOffset{$chr} + $pos) / $imageWindow)}{int($freq1 * $imageY)}++;
							$hetCountsTop{int(($chrTopOffset{$chr} + $pos) / $imageWindow)}{int($freq2 * $imageY)}++;
						}
						elsif (exists $chrBotOffset{$chr})
						{
							$hetCountsBot{int(($chrBotOffset{$chr} + $pos) / $imageWindow)}{int($freq1 * $imageY)}++;
							$hetCountsBot{int(($chrBotOffset{$chr} + $pos) / $imageWindow)}{int($freq2 * $imageY)}++;
						}

					}
				}
			}
		}
	}
}

close FILE;

for ($freq = 0; $freq <= $imageY; $freq++)
{
	for ($pos = 0; $pos <= int($chrTopOffset{end} / $imageWindow); $pos++)
	{
		if (exists $hetCountsTop{$pos}{$freq})
		{
			push(@hetFreqsTop, $hetCountsTop{$pos}{$freq});
		}
		else
		{
			push(@hetFreqsTop, 0);
		}
	}
}

for ($freq = 0; $freq <= $imageY; $freq++)
{
	for ($pos = 0; $pos <= int($chrBotOffset{end} / $imageWindow); $pos++)
	{
		if (exists $hetCountsBot{$pos}{$freq})
		{
			push(@hetFreqsBot, $hetCountsBot{$pos}{$freq});
		}
		else
		{
			push(@hetFreqsBot, 0);
		}
	}
}




# pull cn from HMMcopy
my (@cnStateTop, @cnStartTop, @cnEndTop, @cnStateBot, @cnStartBot, @cnEndBot, @cnColTop, @cnColBot);
my %cnPoints;

my $cnWindow = 10000;

my %cnStateType;

my $file = $dataPath{cn_states};
$file =~ s/\.sorted\.bed\.gz$//;
chomp $file;

warn "Opening $file\n";
open (FILE, $file) or die "Couldn't open $file\n";

my $cnMax = 12;

my %cnColHash = (
    "0 copy" => "#3D52A1",
    "1 copy" => "#77B7E5",
    "2 copy" => "#E6F5FE",
    "3-4 copy" => "#FFE3AA",
    "5-6 copy" => "#F9BD7E",
    "7-8 copy" => "#ED875E",
    "9-10 copy" => "#D24D3E",
    ">10 copy" => "#AE1C3E",
);

while ($l = <FILE>)
{
	unless ($l =~ /^chr\t/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		if ($f[0] eq "chr23")
		{
			$f[0] = "chrX";
		}
		if ($f[0] eq "chr24")
		{
			$f[0] = "chrY";
		}

		if (exists $chrTopOffset{$f[0]})
		{
			push (@cnStartTop, $f[1] + $chrTopOffset{$f[0]});
			push (@cnEndTop, $f[2] + $chrTopOffset{$f[0]});
			if ($f[11] > $cnMax)
			{
				$f[11] = $cnMax;
			}
			if ($f[11] < 0)
			{
				$f[11] = 0;
			}
			push (@cnStateTop, $f[11]);
			

			if ($f[11] < 0.5)
			{
				push (@cnColTop, $cnColHash{"0 copy"});
			}
			elsif ($f[11] < 1.5)
			{
				push (@cnColTop, $cnColHash{"1 copy"});
			}
			elsif ($f[11] < 2.5)
			{
				push (@cnColTop, $cnColHash{"2 copy"});
			}
			elsif ($f[11] < 4.5)
			{
				push (@cnColTop, $cnColHash{"3-4 copy"});
			}
			elsif ($f[11] < 6.5)
			{
				push (@cnColTop, $cnColHash{"5-6 copy"});
			}
			elsif ($f[11] < 8.5)
			{
				push (@cnColTop, $cnColHash{"7-8 copy"});
			}
			elsif ($f[11] < 10.5)
			{
				push (@cnColTop, $cnColHash{"9-10 copy"});
			}
			else
			{
				push (@cnColTop, $cnColHash{">10 copy"});
			}
		}
		elsif (exists $chrBotOffset{$f[0]})
		{
			push (@cnStartBot, $f[1] + $chrBotOffset{$f[0]});
			push (@cnEndBot, $f[2] + $chrBotOffset{$f[0]});
			if ($f[11] > $cnMax)
			{
				$f[11] = $cnMax;
			}
			if ($f[11] < 0)
			{
				$f[11] = 0;
			}
			push (@cnStateBot, $f[11]);

			if ($f[11] < 0.5)
			{
				push (@cnColBot, $cnColHash{"0 copy"});
			}
			elsif ($f[11] < 1.5)
			{
				push (@cnColBot, $cnColHash{"1 copy"});
			}
			elsif ($f[11] < 2.5)
			{
				push (@cnColBot, $cnColHash{"2 copy"});
			}
			elsif ($f[11] < 4.5)
			{
				push (@cnColBot, $cnColHash{"3-4 copy"});
			}
			elsif ($f[11] < 6.5)
			{
				push (@cnColBot, $cnColHash{"5-6 copy"});
			}
			elsif ($f[11] < 8.5)
			{
				push (@cnColBot, $cnColHash{"7-8 copy"});
			}
			elsif ($f[11] < 10.5)
			{
				push (@cnColBot, $cnColHash{"9-10 copy"});
			}
			else
			{
				push (@cnColBot, $cnColHash{">10 copy"});
			}
		}

		push(@{ $cnStateType{$f[0]}{pos} }, $f[2]);
		push(@{ $cnStateType{$f[0]}{state} }, $f[11]);

	}
}
close FILE;





# get SV from crest
my (@breakX1, @breakX2, @breakY1, @breakY2, @breakCol);
my ($leftChr,$leftPos,$rightChr,$rightPos);
my $bpWidth = 1500000000;

my (@topBreakX1, @topBreakX2, @topBreakY1, @topBreakY2, @topBreakCol);
my (@botBreakX1, @botBreakX2, @botBreakY1, @botBreakY2, @botBreakCol);

warn "Opening $dataPath{structural}\n";
open (FILE, "$dataPath{structural}") or die "Couldn't open $dataPath{structural}\n";

while ($l = <FILE>)
{
	chomp $l;

	unless ($l =~ /^chrom1/)
	{
		$count{sv}++;
		@f = split(/\t/, $l);
	
		$leftChr = $f[0];
		$leftPos = $f[1];
		$rightChr = $f[2];
		$rightPos = $f[3];
		$type = $f[4];
	
		if ((exists $chrTopOffset{$leftChr}) and (exists $chrTopOffset{$rightChr}))
		{
			$leftPos = (($leftPos + $chrTopOffset{$leftChr}) / $chrTopOffset{end}) * $bpWidth;
			$rightPos = (($rightPos + $chrTopOffset{$rightChr}) / $chrTopOffset{end}) * $bpWidth;
	
			push(@breakX1,$leftPos);
			push(@breakY1, 1);
	
			push(@breakX2, ($leftPos + $rightPos) / 2);
			push(@breakY2, 0.5);
			
			push(@breakX1, ($leftPos + $rightPos) / 2);
			push(@breakY1, 0.5);
	
			push(@breakX2, $rightPos);
			push(@breakY2, 1);
	
			push(@breakCol,$bpColours{$type});
			push(@breakCol,$bpColours{$type});

			$leftPos = $f[1] + $chrTopOffset{$leftChr};
			$rightPos = $f[3] + $chrTopOffset{$rightChr};

			unless ($type eq "TRA")
			{
				push(@topBreakX1, $leftPos);
				push(@topBreakY1, 0);
	
				push(@topBreakX2, ($leftPos + $rightPos) / 2);
				push(@topBreakY2, 1);

				push(@topBreakX1, ($leftPos + $rightPos) / 2);
				push(@topBreakY1, 1);

				push(@topBreakX2, $rightPos);
				push(@topBreakY2, 0);

				push(@topBreakCol,$bpColours{$type});
				push(@topBreakCol,$bpColours{$type});
			}
			else # TRANSLOCATION
			{
				push(@topBreakX1, $leftPos);
				push(@topBreakY1, 0);
				push(@topBreakX2, $leftPos);
				push(@topBreakY2, 1);

				push(@topBreakX1, $rightPos);
				push(@topBreakY1, 0);
				push(@topBreakX2, $rightPos);
				push(@topBreakY2, 1);

				push(@topBreakCol,$bpColours{$type});
				push(@topBreakCol,$bpColours{$type});
			}
		}
		elsif ((exists $chrBotOffset{$leftChr}) and (exists $chrBotOffset{$rightChr}))
		{
			$leftPos = (($leftPos + $chrBotOffset{$leftChr}) / $chrBotOffset{end}) * $bpWidth;
			$rightPos = (($rightPos + $chrBotOffset{$rightChr}) / $chrBotOffset{end}) * $bpWidth;
	
			push(@breakX1,$leftPos);
			push(@breakY1, 0);
	
			push(@breakX2, ($leftPos + $rightPos) / 2);
			push(@breakY2, 0.5);
			
			push(@breakX1, ($leftPos + $rightPos) / 2);
			push(@breakY1, 0.5);
	
			push(@breakX2, $rightPos);
			push(@breakY2, 0);
	
			push(@breakCol,$bpColours{$type});
			push(@breakCol,$bpColours{$type});

			$leftPos = $f[1] + $chrBotOffset{$leftChr};
			$rightPos = $f[3] + $chrBotOffset{$rightChr};
			unless ($type eq "TRA")
			{
				push(@botBreakX1, $leftPos);
				push(@botBreakY1, 0);
	
				push(@botBreakX2, ($leftPos + $rightPos) / 2);
				push(@botBreakY2, 1);

				push(@botBreakX1, ($leftPos + $rightPos) / 2);
				push(@botBreakY1, 1);

				push(@botBreakX2, $rightPos);
				push(@botBreakY2, 0);

				push(@botBreakCol,$bpColours{$type});
				push(@botBreakCol,$bpColours{$type});
			}
			else # TRANSLOCATION
			{
				push(@botBreakX1, $leftPos);
				push(@botBreakY1, 0);
				push(@botBreakX2, $leftPos);
				push(@botBreakY2, 1);

				push(@botBreakX1, $rightPos);
				push(@botBreakY1, 0);
				push(@botBreakX2, $rightPos);
				push(@botBreakY2, 1);

				push(@botBreakCol,$bpColours{$type});
				push(@botBreakCol,$bpColours{$type});
			}
		}
		elsif ((exists $chrTopOffset{$leftChr}) and (exists $chrBotOffset{$rightChr}))
		{
			$leftPos = (($leftPos + $chrTopOffset{$leftChr}) / $chrTopOffset{end}) * $bpWidth;
			$rightPos = (($rightPos + $chrBotOffset{$rightChr}) / $chrBotOffset{end}) * $bpWidth;
	
			push(@breakX1,$leftPos);
			push(@breakY1, 1);
	
			push(@breakX2, ($leftPos + $rightPos) / 2);
			push(@breakY2, 0.5);
			
			push(@breakX1, ($leftPos + $rightPos) / 2);
			push(@breakY1, 0.5);
	
			push(@breakX2, $rightPos);
			push(@breakY2, 0);
	
			push(@breakCol,$bpColours{$type});
			push(@breakCol,$bpColours{$type});


			$leftPos = $f[1] + $chrTopOffset{$leftChr};
			$rightPos = $f[3] + $chrBotOffset{$rightChr};
			# assuming TRA
			push(@topBreakX1, $leftPos);
			push(@topBreakY1, 0);
			push(@topBreakX2, $leftPos);
			push(@topBreakY2, 1);

			push(@botBreakX1, $rightPos);
			push(@botBreakY1, 0);
			push(@botBreakX2, $rightPos);
			push(@botBreakY2, 1);

			push(@topBreakCol,$bpColours{$type});
			push(@botBreakCol,$bpColours{$type});

		}
		elsif ((exists $chrBotOffset{$leftChr}) and (exists $chrTopOffset{$rightChr}))
		{
			$leftPos = (($leftPos + $chrBotOffset{$leftChr}) / $chrBotOffset{end}) * $bpWidth;
			$rightPos = (($rightPos + $chrTopOffset{$rightChr}) / $chrTopOffset{end}) * $bpWidth;
	
			push(@breakX1,$leftPos);
			push(@breakY1, 0);
	
			push(@breakX2, ($leftPos + $rightPos) / 2);
			push(@breakY2, 0.5);
			
			push(@breakX1, ($leftPos + $rightPos) / 2);
			push(@breakY1, 0.5);
			
			push(@breakX2, $rightPos);
			push(@breakY2, 1);
	
			push(@breakCol,$bpColours{$type});
			push(@breakCol,$bpColours{$type});

			$leftPos = $f[1] + $chrBotOffset{$leftChr};
			$rightPos = $f[3] + $chrTopOffset{$rightChr};
			# assuming TRA
			push(@topBreakX1, $rightPos);
			push(@topBreakY1, 0);
			push(@topBreakX2, $rightPos);
			push(@topBreakY2, 1);

			push(@botBreakX1, $leftPos);
			push(@botBreakY1, 0);
			push(@botBreakX2, $leftPos);
			push(@botBreakY2, 1);

			push(@topBreakCol,$bpColours{$type});
			push(@botBreakCol,$bpColours{$type});
		}
	}

}
close FILE;


# make rainfall plot
my @positionsTop;
my @positionsBot;

my @distancesTop;
my @distancesBot;

my @coloursTop;
my @coloursBot;

my $rainMaxTop = 0;
my $rainMaxBot = 0;

for $chr (keys %vcfHash)
{
	for $pos (keys %{ $vcfHash{$chr} })
	{
		for $alt (keys %{ $vcfHash{$chr}{$pos}{snv} })
		{
			if (exists $chrTopOffset{$chr})
			{
				if (exists $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist})
				{
					push(@positionsTop, $pos + $chrTopOffset{$chr});
					push(@distancesTop, $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist});
					push(@coloursTop, $colours{$vcfHash{$chr}{$pos}{snv}{$alt}{change}});

					if ($vcfHash{$chr}{$pos}{snv}{$alt}{prevdist} > $rainMaxTop)
					{
						$rainMaxTop = $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist};
					}
				}
			}
			elsif (exists $chrBotOffset{$chr})
			{
				if (exists $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist})
				{
					push(@positionsBot, $pos + $chrBotOffset{$chr});
					push(@distancesBot, $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist});
					push(@coloursBot, $colours{$vcfHash{$chr}{$pos}{snv}{$alt}{change}});

					if ($vcfHash{$chr}{$pos}{snv}{$alt}{prevdist} > $rainMaxBot)
					{
						$rainMaxBot = $vcfHash{$chr}{$pos}{snv}{$alt}{prevdist};
					}
				}
			}
		}
	}
}

my @verticalsTop;
my @verticalsBot;
for $chr (keys %chrTopOffset)
{
	push(@verticalsTop, $chrTopOffset{$chr});
}
for $chr (keys %chrBotOffset)
{
	push(@verticalsBot, $chrBotOffset{$chr});
}

my @chrListTop = qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 end/;
my @atLabelsTop;
my @chrLabelsTop;
for (my $i = 0; $i < (scalar @chrListTop - 1); $i++)
{
	push(@atLabelsTop, ($chrTopOffset{$chrListTop[$i]} + $chrTopOffset{$chrListTop[$i + 1]}) / 2);

	push(@chrLabelsTop, $chrListTop[$i]);
	$chrLabelsTop[$i] =~ s/chr//;
}

my @chrListBot = qw/chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY end/;
my @atLabelsBot;
my @chrLabelsBot;
for (my $i = 0; $i < (scalar @chrListBot - 1); $i++)
{
	push(@atLabelsBot, ($chrBotOffset{$chrListBot[$i]} + $chrBotOffset{$chrListBot[$i + 1]}) / 2);

	push(@chrLabelsBot, $chrListBot[$i]);
	$chrLabelsBot[$i] =~ s/chr//;
}



# read karotypes
$file = "/u/rdenroche/hg19_bands.txt";
open (FILE, $file) or die "Couldn't open $file\n";

my (@karoStartTop, @karoStartBot, @karoEndTop, @karoEndBot, @karoColTop, @karoColBot);

my %centroPos;
my $centroStart;

my %karoColour = (
	"acen" => "#CC6677",
	"gpos25" => "grey75",
	"gpos50" => "grey50",
	"gpos75" => "grey25",
	"gpos100" => "grey0"
);

while ($l = <FILE>)
{
	unless ($l =~ /^#/)
	{
		chomp $l;
		@f = split(/\t/, $l);

		if (exists $chrTopOffset{$f[0]})
		{
			push(@karoStartTop, $f[1] + $chrTopOffset{$f[0]});
			push(@karoEndTop, $f[2] + $chrTopOffset{$f[0]});
			if (exists $karoColour{$f[4]})
			{
				push(@karoColTop, $karoColour{$f[4]});
			}
			else
			{
				push(@karoColTop, "white");
			}

            if ($f[4] eq "acen")
            {
                if ($f[3] =~ /p/)
                {
                    $centroStart = $f[1];
                }
                else
                {
                    $centroPos{$f[0]} = (($centroStart + $f[2]) / 2) + $chrTopOffset{$f[0]};
                }
            }

		}
		elsif (exists $chrBotOffset{$f[0]})
		{
			push(@karoStartBot, $f[1] + $chrBotOffset{$f[0]});
			push(@karoEndBot, $f[2] + $chrBotOffset{$f[0]});
			if (exists $karoColour{$f[4]})
			{
				push(@karoColBot, $karoColour{$f[4]});
			}
			else
			{
				push(@karoColBot, "white");
			}

            if ($f[4] eq "acen")
            {
                if ($f[3] =~ /p/)
                {
                    $centroStart = $f[1];
                }
                else
                {
                    $centroPos{$f[0]} = (($centroStart + $f[2]) / 2) + $chrBotOffset{$f[0]};
                }
            }

		}


	}
}
close FILE;











warn "Writing to $outputFile.R\n";
open (RFILE, ">$outputFile.R") or die "Couldn't open >$outputFile.R\n";

print RFILE "cat(\"R: Loading variables\\n\")\n";

print RFILE "xvalsTop <- c($positionsTop[0]";
for (my $i = 1; $i < scalar @positionsTop; $i++)
{
	print RFILE ", \"$positionsTop[$i]\"";
}
print RFILE ")\n";

print RFILE "yvalsTop <- c($distancesTop[0]";
for (my $i = 1; $i < scalar @distancesTop; $i++)
{
	print RFILE ", \"$distancesTop[$i]\"";
}
print RFILE ")\n";

print RFILE "colsTop <- c(\"$coloursTop[0]\"";
for (my $i = 1; $i < scalar @coloursTop; $i++)
{
	print RFILE ", \"$coloursTop[$i]\"";
}
print RFILE ")\n";

print RFILE "vertsTop <- c($verticalsTop[0]";
for (my $i = 1; $i < scalar @verticalsTop; $i++)
{
	print RFILE ", $verticalsTop[$i]";
}
print RFILE ")\n";

print RFILE "hetVertsTop <- c($verticalsTop[0]";
for (my $i = 1; $i < scalar @verticalsTop; $i++)
{
	print RFILE ", " . int($verticalsTop[$i] / $imageWindow);
}
print RFILE ")\n";

print RFILE "atLabelsTop <- c($atLabelsTop[0]";
for (my $i = 1; $i < scalar @atLabelsTop; $i++)
{
	print RFILE ", $atLabelsTop[$i]";
}
print RFILE ")\n";

print RFILE "chrLabelsTop <- c(\"$chrLabelsTop[0]\"";
for (my $i = 1; $i < scalar @chrLabelsTop; $i++)
{
	print RFILE ", \"$chrLabelsTop[$i]\"";
}
print RFILE ")\n";

print RFILE "chrLabelsTopS <- c(\"chr$chrLabelsTop[0]\"";
for (my $i = 1; $i < scalar @chrLabelsTop; $i++)
{
	print RFILE ", \"chr$chrLabelsTop[$i]\"";
}
print RFILE ")\n";



print RFILE "xvalsBot <- c($positionsBot[0]";
for (my $i = 1; $i < scalar @positionsBot; $i++)
{
	print RFILE ", \"$positionsBot[$i]\"";
}
print RFILE ")\n";

print RFILE "yvalsBot <- c($distancesBot[0]";
for (my $i = 1; $i < scalar @distancesBot; $i++)
{
	print RFILE ", \"$distancesBot[$i]\"";
}
print RFILE ")\n";

print RFILE "colsBot <- c(\"$coloursBot[0]\"";
for (my $i = 1; $i < scalar @coloursBot; $i++)
{
	print RFILE ", \"$coloursBot[$i]\"";
}
print RFILE ")\n";

print RFILE "vertsBot <- c($verticalsBot[0]";
for (my $i = 1; $i < scalar @verticalsBot; $i++)
{
	print RFILE ", $verticalsBot[$i]";
}
print RFILE ")\n";

print RFILE "hetVertsBot <- c($verticalsBot[0]";
for (my $i = 1; $i < scalar @verticalsBot; $i++)
{
	print RFILE ", " . int(($verticalsBot[$i] / $imageWindow));
}
print RFILE ")\n";

print RFILE "atLabelsBot <- c($atLabelsBot[0]";
for (my $i = 1; $i < scalar @atLabelsBot; $i++)
{
	print RFILE ", $atLabelsBot[$i]";
}
print RFILE ")\n";

print RFILE "chrLabelsBot <- c(\"$chrLabelsBot[0]\"";
for (my $i = 1; $i < scalar @chrLabelsBot; $i++)
{
	print RFILE ", \"$chrLabelsBot[$i]\"";
}
print RFILE ")\n";

print RFILE "chrLabelsBotS <- c(\"chr$chrLabelsBot[0]\"";
for (my $i = 1; $i < scalar @chrLabelsBot; $i++)
{
	print RFILE ", \"chr$chrLabelsBot[$i]\"";
}
print RFILE ")\n";




# het freq plot R values
print RFILE "hetFreqsTop <- c($hetFreqsTop[0]";
for (my $i = 1; $i < scalar @hetFreqsTop; $i++)
{
	print RFILE ", $hetFreqsTop[$i]";
}
print RFILE ")\n";

print RFILE "hetFreqsBot <- c($hetFreqsBot[0]";
for (my $i = 1; $i < scalar @hetFreqsBot; $i++)
{
	print RFILE ", $hetFreqsBot[$i]";
}
print RFILE ")\n";


print RFILE "hetTopX <- c(0";
for (my $i = $imageWindow; $i <= int($chrTopOffset{end} / $imageWindow) * $imageWindow; $i += $imageWindow)
{
	print RFILE ", $i";
}
print RFILE ")\n";

print RFILE "hetBotX <- c(0";
for (my $i = $imageWindow; $i <= int($chrBotOffset{end} / $imageWindow) * $imageWindow; $i += $imageWindow)
{
	print RFILE ", $i";
}
print RFILE ")\n";

print RFILE "hetY <- c(0";
for (my $i = 1; $i <= $imageY; $i++)
{
	print RFILE ", " . ($i / $imageY);
}
print RFILE ")\n";


#print RFILE "hetPosTop <- c($hetPosTop[0]";
#for (my $i = 1; $i < scalar @hetPosTop; $i++)
#{
#	print RFILE ", $hetPosTop[$i]";
#}
#print RFILE ")\n";
#
#print RFILE "hetPosBot <- c($hetPosBot[0]";
#for (my $i = 1; $i < scalar @hetPosBot; $i++)
#{
#	print RFILE ", $hetPosBot[$i]";
#}
#print RFILE ")\n";



# cn state plot R values

print RFILE "cnStateTop <- c($cnStateTop[0]";
for (my $i = 1; $i < scalar @cnStateTop; $i++)
{
	print RFILE ", $cnStateTop[$i]";
}
print RFILE ")\n";

print RFILE "cnColTop <- c(\"$cnColTop[0]\"";
for (my $i = 1; $i < scalar @cnColTop; $i++)
{
	print RFILE ", \"$cnColTop[$i]\"";
}
print RFILE ")\n";

print RFILE "cnStartTop <- c($cnStartTop[0]";
for (my $i = 1; $i < scalar @cnStartTop; $i++)
{
	print RFILE ", $cnStartTop[$i]";
}
print RFILE ")\n";

print RFILE "cnEndTop <- c($cnEndTop[0]";
for (my $i = 1; $i < scalar @cnEndTop; $i++)
{
	print RFILE ", $cnEndTop[$i]";
}
print RFILE ")\n";

print RFILE "cnStateBot <- c($cnStateBot[0]";
for (my $i = 1; $i < scalar @cnStateBot; $i++)
{
	print RFILE ", $cnStateBot[$i]";
}
print RFILE ")\n";

print RFILE "cnColBot <- c(\"$cnColBot[0]\"";
for (my $i = 1; $i < scalar @cnColBot; $i++)
{
	print RFILE ", \"$cnColBot[$i]\"";
}
print RFILE ")\n";

print RFILE "cnStartBot <- c($cnStartBot[0]";
for (my $i = 1; $i < scalar @cnStartBot; $i++)
{
	print RFILE ", $cnStartBot[$i]";
}
print RFILE ")\n";

print RFILE "cnEndBot <- c($cnEndBot[0]";
for (my $i = 1; $i < scalar @cnEndBot; $i++)
{
	print RFILE ", $cnEndBot[$i]";
}
print RFILE ")\n";



# karo bands
#my (@karoStartTop, @karoStartBot, @karoEndTop, @karoEndBot, @karoColTop, @karoColBot);

print RFILE "karoStartTop <- c($karoStartTop[0]";
for (my $i = 1; $i < scalar @karoStartTop; $i++)
{
	print RFILE ", $karoStartTop[$i]";
}
print RFILE ")\n";

print RFILE "karoEndTop <- c($karoEndTop[0]";
for (my $i = 1; $i < scalar @karoEndTop; $i++)
{
	print RFILE ", $karoEndTop[$i]";
}
print RFILE ")\n";

print RFILE "karoColTop <- c(\"$karoColTop[0]\"";
for (my $i = 1; $i < scalar @karoColTop; $i++)
{
	print RFILE ", \"$karoColTop[$i]\"";
}
print RFILE ")\n";

print RFILE "karoStartBot <- c($karoStartBot[0]";
for (my $i = 1; $i < scalar @karoStartBot; $i++)
{
	print RFILE ", $karoStartBot[$i]";
}
print RFILE ")\n";

print RFILE "karoEndBot <- c($karoEndBot[0]";
for (my $i = 1; $i < scalar @karoEndBot; $i++)
{
	print RFILE ", $karoEndBot[$i]";
}
print RFILE ")\n";

print RFILE "karoColBot <- c(\"$karoColBot[0]\"";
for (my $i = 1; $i < scalar @karoColBot; $i++)
{
	print RFILE ", \"$karoColBot[$i]\"";
}
print RFILE ")\n";

my $breakWidth = 7000000;
my $gapWidth = 5000000;

my @karoBreakTop;
my @karoBreakBot;

my @chrBreakTop;
my @chrBreakBot;

for $chr (qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8/)
{
    push (@karoBreakTop, $centroPos{$chr});
}
for $chr (qw/chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY/)
{
    push (@karoBreakBot, $centroPos{$chr});
}
for $chr (qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 end/)
{
    push (@chrBreakTop, $chrTopOffset{$chr});
}
for $chr (qw/chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY end/)
{
    push (@chrBreakBot, $chrBotOffset{$chr});
}


# triangles
print RFILE "karoBreakTop <- c($karoBreakTop[0]," . ($karoBreakTop[0] - $breakWidth) . "," . ($karoBreakTop[0] + $breakWidth) . ",NA,$karoBreakTop[0]," . ($karoBreakTop[0] - $breakWidth) . "," . ($karoBreakTop[0] + $breakWidth) . ",NA";
for (my $i = 1; $i < scalar @karoBreakTop; $i++)
{
	print RFILE ",$karoBreakTop[$i]," . ($karoBreakTop[$i] - $breakWidth) . "," . ($karoBreakTop[$i] + $breakWidth) . ",NA,$karoBreakTop[$i]," . ($karoBreakTop[$i] - $breakWidth) . "," . ($karoBreakTop[$i] + $breakWidth) . ",NA";
}
print RFILE ")\n";

print RFILE "karoBreakTopY <- c(0.5,1,1,NA,0.5,0,0,NA";
for (my $i = 1; $i < scalar @karoBreakTop; $i++)
{
	print RFILE ",0.5,1,1,NA,0.5,0,0,NA";
}
print RFILE ")\n";

print RFILE "karoBreakBot <- c($karoBreakBot[0]," . ($karoBreakBot[0] - $breakWidth) . "," . ($karoBreakBot[0] + $breakWidth) . ",NA,$karoBreakBot[0]," . ($karoBreakBot[0] - $breakWidth) . "," . ($karoBreakBot[0] + $breakWidth) . ",NA";
for (my $i = 1; $i < scalar @karoBreakBot; $i++)
{
	print RFILE ",$karoBreakBot[$i]," . ($karoBreakBot[$i] - $breakWidth) . "," . ($karoBreakBot[$i] + $breakWidth) . ",NA,$karoBreakBot[$i]," . ($karoBreakBot[$i] - $breakWidth) . "," . ($karoBreakBot[$i] + $breakWidth) . ",NA";
}
print RFILE ")\n";

print RFILE "karoBreakBotY <- c(0.5,1,1,NA,0.5,0,0,NA";
for (my $i = 1; $i < scalar @karoBreakBot; $i++)
{
	print RFILE ",0.5,1,1,NA,0.5,0,0,NA";
}
print RFILE ")\n";


# chrom spaces
print RFILE "chromSpaceStartTop <- c(0";
for ($chr = 1; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakTop[$chr] - $gapWidth);
}
print RFILE ")\n";

print RFILE "chromSpaceEndTop <- c($gapWidth";
for ($chr = 1; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakTop[$chr] + $gapWidth);
}
print RFILE ")\n";


print RFILE "chromSpaceStartBot <- c(0";
for ($chr = 1; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakBot[$chr] - $gapWidth);
}
print RFILE ")\n";

print RFILE "chromSpaceEndBot <- c($gapWidth";
for ($chr = 1; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakBot[$chr] + $gapWidth);
}
print RFILE ")\n";




# border
my ($start, $end, $centStart, $centEnd, $top, $bot);

print RFILE "karoBorderTopX <- c(";
for ($chr = 0; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    $start = $chrBreakTop[$chr] + $gapWidth;
    $centStart = $karoBreakTop[$chr] - $breakWidth;
    $centEnd = $karoBreakTop[$chr] + $breakWidth;
    $end = $chrBreakTop[$chr + 1] - $gapWidth;

    print RFILE "$start, $centStart, $centEnd, $end, $end, $centEnd, $centStart, $start, $start, NA, ";
}
print RFILE "NA)\n";

print RFILE "karoBorderTopY <- c(";
for ($chr = 0; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    $top = 1;
    $bot = 0;

    print RFILE "$top, $top, $bot, $bot, $top, $top, $bot, $bot, $top, NA, ";
}
print RFILE "NA)\n";


print RFILE "karoBorderBotX <- c(";
for ($chr = 0; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    $start = $chrBreakBot[$chr] + $gapWidth;
    $centStart = $karoBreakBot[$chr] - $breakWidth;
    $centEnd = $karoBreakBot[$chr] + $breakWidth;
    $end = $chrBreakBot[$chr + 1] - $gapWidth;

    print RFILE "$start, $centStart, $centEnd, $end, $end, $centEnd, $centStart, $start, $start, NA, ";
}
print RFILE "NA)\n";

print RFILE "karoBorderBotY <- c(";
for ($chr = 0; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    $top = 1;
    $bot = 0;

    print RFILE "$top, $top, $bot, $bot, $top, $top, $bot, $bot, $top, NA, ";
}
print RFILE "NA)\n";


my $breakWidthS = 700000;
my $gapWidthS = 500000;
# chrom spaces for singles
print RFILE "chromSpaceStartTopS <- c(0";
for ($chr = 1; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakTop[$chr] - $gapWidthS);
}
print RFILE ")\n";

print RFILE "chromSpaceEndTopS <- c($gapWidthS";
for ($chr = 1; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakTop[$chr] + $gapWidthS);
}
print RFILE ")\n";


print RFILE "chromSpaceStartBotS <- c(0";
for ($chr = 1; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakBot[$chr] - $gapWidthS);
}
print RFILE ")\n";

print RFILE "chromSpaceEndBotS <- c($gapWidth";
for ($chr = 1; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    print RFILE ", " . ($chrBreakBot[$chr] + $gapWidthS);
}
print RFILE ")\n";




# border for singles
print RFILE "karoBorderTopXS <- c(";
for ($chr = 0; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    $start = $chrBreakTop[$chr] + $gapWidthS;
    $centStart = $karoBreakTop[$chr] - $breakWidth;
    $centEnd = $karoBreakTop[$chr] + $breakWidth;
    $end = $chrBreakTop[$chr + 1] - $gapWidthS;

    print RFILE "$start, $centStart, $centEnd, $end, $end, $centEnd, $centStart, $start, $start, NA, ";
}
print RFILE "NA)\n";

print RFILE "karoBorderTopYS <- c(";
for ($chr = 0; $chr < scalar(@chrBreakTop) - 1; $chr++)
{
    $top = 1;
    $bot = 0;

    print RFILE "$top, $top, $bot, $bot, $top, $top, $bot, $bot, $top, NA, ";
}
print RFILE "NA)\n";


print RFILE "karoBorderBotXS <- c(";
for ($chr = 0; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    $start = $chrBreakBot[$chr] + $gapWidthS;
    $centStart = $karoBreakBot[$chr] - $breakWidth;
    $centEnd = $karoBreakBot[$chr] + $breakWidth;
    $end = $chrBreakBot[$chr + 1] - $gapWidthS;

    print RFILE "$start, $centStart, $centEnd, $end, $end, $centEnd, $centStart, $start, $start, NA, ";
}
print RFILE "NA)\n";

print RFILE "karoBorderBotYS <- c(";
for ($chr = 0; $chr < scalar(@chrBreakBot) - 1; $chr++)
{
    $top = 1;
    $bot = 0;

    print RFILE "$top, $top, $bot, $bot, $top, $top, $bot, $bot, $top, NA, ";
}
print RFILE "NA)\n";




# print sv

print RFILE "breakX1 <- c($breakX1[0]";
for (my $i = 1; $i < scalar @breakX1; $i++)
{
	print RFILE ", $breakX1[$i]";
}
print RFILE ")\n";

print RFILE "breakY1 <- c($breakY1[0]";
for (my $i = 1; $i < scalar @breakY1; $i++)
{
	print RFILE ", $breakY1[$i]";
}
print RFILE ")\n";

print RFILE "breakX2 <- c($breakX2[0]";
for (my $i = 1; $i < scalar @breakX2; $i++)
{
	print RFILE ", $breakX2[$i]";
}
print RFILE ")\n";

print RFILE "breakY2 <- c($breakY2[0]";
for (my $i = 1; $i < scalar @breakY2; $i++)
{
	print RFILE ", $breakY2[$i]";
}
print RFILE ")\n";

print RFILE "breakCol <- c(\"$breakCol[0]\"";
for (my $i = 1; $i < scalar @breakCol; $i++)
{
	print RFILE ", \"$breakCol[$i]\"";
}
print RFILE ")\n";

print RFILE "topBreakX1 <- c($topBreakX1[0]";
for (my $i = 1; $i < scalar @topBreakX1; $i++)
{
	print RFILE ", $topBreakX1[$i]";
}
print RFILE ")\n";

print RFILE "topBreakY1 <- c($topBreakY1[0]";
for (my $i = 1; $i < scalar @topBreakY1; $i++)
{
	print RFILE ", $topBreakY1[$i]";
}
print RFILE ")\n";

print RFILE "topBreakX2 <- c($topBreakX2[0]";
for (my $i = 1; $i < scalar @topBreakX2; $i++)
{
	print RFILE ", $topBreakX2[$i]";
}
print RFILE ")\n";

print RFILE "topBreakY2 <- c($topBreakY2[0]";
for (my $i = 1; $i < scalar @topBreakY2; $i++)
{
	print RFILE ", $topBreakY2[$i]";
}
print RFILE ")\n";

print RFILE "topBreakCol <- c(\"$topBreakCol[0]\"";
for (my $i = 1; $i < scalar @topBreakCol; $i++)
{
	print RFILE ", \"$topBreakCol[$i]\"";
}
print RFILE ")\n";

print RFILE "botBreakX1 <- c($botBreakX1[0]";
for (my $i = 1; $i < scalar @botBreakX1; $i++)
{
	print RFILE ", $botBreakX1[$i]";
}
print RFILE ")\n";

print RFILE "botBreakY1 <- c($botBreakY1[0]";
for (my $i = 1; $i < scalar @botBreakY1; $i++)
{
	print RFILE ", $botBreakY1[$i]";
}
print RFILE ")\n";

print RFILE "botBreakX2 <- c($botBreakX2[0]";
for (my $i = 1; $i < scalar @botBreakX2; $i++)
{
	print RFILE ", $botBreakX2[$i]";
}
print RFILE ")\n";

print RFILE "botBreakY2 <- c($botBreakY2[0]";
for (my $i = 1; $i < scalar @botBreakY2; $i++)
{
	print RFILE ", $botBreakY2[$i]";
}
print RFILE ")\n";

print RFILE "botBreakCol <- c(\"$botBreakCol[0]\"";
for (my $i = 1; $i < scalar @botBreakCol; $i++)
{
	print RFILE ", \"$botBreakCol[$i]\"";
}
print RFILE ")\n";


# genes of interest ticks
my $geneBump = 18000000;
if (scalar @geneNameTop > 0)
{
	print RFILE "geneNameTop <- c(\"$geneNameTop[0]\"";
	for (my $i = 1; $i < scalar @geneNameTop; $i++)
	{
		print RFILE ", \"$geneNameTop[$i]\"";
	}
	print RFILE ")\n";
	
	print RFILE "genePosTop <- c($genePosTop[0]";
	for (my $i = 1; $i < scalar @genePosTop; $i++)
	{
		if ($genePosTop[$i - 1] + $geneBump > $genePosTop[$i])
		{
			$genePosTop[$i] += $geneBump - ($genePosTop[$i] - $genePosTop[$i - 1]);
		}
		print RFILE ", \"$genePosTop[$i]\"";
	}
	print RFILE ")\n";
}

if (scalar @geneNameBot > 0)
{
	print RFILE "geneNameBot <- c(\"$geneNameBot[0]\"";
	for (my $i = 1; $i < scalar @geneNameBot; $i++)
	{
		print RFILE ", \"$geneNameBot[$i]\"";
	}
	print RFILE ")\n";
	
	print RFILE "genePosBot <- c($genePosBot[0]";
	for (my $i = 1; $i < scalar @genePosBot; $i++)
	{
		if ($genePosBot[$i - 1] + $geneBump > $genePosBot[$i])
		{
			$genePosBot[$i] += $geneBump - ($genePosBot[$i] - $genePosBot[$i - 1]);
		}
		print RFILE ", \"$genePosBot[$i]\"";
	}
	print RFILE ")\n";
}





my $headerString = "   $tSample (${tCoverage}x) // $nSample (${nCoverage}x)     $count{snv} somatic SNVs, " . ($count{ins} + $count{del}) . " somatic indels, $count{sv} structrual variants";
my $footerString = "   $outputFile $date";



print RFILE "cat(\"R: Plotting...\\n\")\n";

print RFILE "png(filename = \"$outName.whole_genome.png\", width = 2150, height = 1650)\n";

print RFILE "par(oma=c(2.5,0,2.5,0),mgp=c(1,0,-1), cex=2)\n";

# top rainfall
print RFILE "par(fig=c(0,1,0.86,1), mar=c(0, 3, 0, 0) + 0.1)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(1,$rainMaxTop), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"Somatic SNV\\nRainfall\",)\n";
print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
print RFILE "points(xvalsTop, yvalsTop, col=colsTop, pch = 20)\n";
print RFILE "axis(side=3, at=atLabelsTop, labels=chrLabelsTop, tick=FALSE)\n";
if (scalar @geneNameTop > 0)
{
	print RFILE "text(genePosTop,1, labels=geneNameTop, pos=4, srt=45, cex=0.75)\n";
}

# top het freqs
print RFILE "par(fig=c(0,1,0.72,0.86), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(0, 1), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"Germline Het\\nFreq in Tumour\",)\n";
print RFILE "image(hetTopX, hetY, matrix(hetFreqsTop," . (int($chrTopOffset{end} / $imageWindow) + 1) . ",$imageY + 1), add=TRUE, col=c(grey.colors(5,start=1,end=0.5),grey.colors(35,start=0.5,end=0)))\n";
print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
#print RFILE "points(hetPosTop, hetFreqsTop, col=\"black\", pch=\".\")\n";

# top cnv
print RFILE "par(fig=c(0,1,0.58,0.72), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(0,$cnMax), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", main=\"\", xlab=\"\", ylab=\"Copy Number\\nVariation\",)\n";
print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
#print RFILE "points(cnPosTop, cnValTop, col=cnColTop, pch=20)\n";
print RFILE "segments(cnStartTop, cnStateTop, cnEndTop, cnStateTop, col=cnColTop, lwd=20)\n";
print RFILE "segments(cnStartTop, cnStateTop, cnEndTop, cnStateTop, col=\"black\", lwd=2)\n";

#top karotype
print RFILE "par(fig=c(0,1,0.555,0.58), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrTopOffset{end}), c(0,1),type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\")\n";
print RFILE "rect(karoStartTop, 0, karoEndTop, 1, col=karoColTop, border=NA)\n";
print RFILE "rect(chromSpaceStartTop, 0, chromSpaceEndTop, 1, col=\"white\", border=NA)\n";
print RFILE "polygon(karoBreakTop,karoBreakTopY,col=\"white\",border=NA)\n";
print RFILE "lines(karoBorderTopX,karoBorderTopY,col=\"black\")\n";

# structural
print RFILE "par(fig=c(0,1,0.445,0.555), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$bpWidth), c(0,1), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"Structural\\nBreak Points\")\n";
print RFILE "segments(breakX1,breakY1,breakX2,breakY2,col=breakCol,lwd=2)\n";

# bottom karotype
print RFILE "par(fig=c(0,1,0.42,0.445), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(0,1), type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\")\n";
print RFILE "rect(karoStartBot, 0, karoEndBot, 1, col=karoColBot, border=NA)\n";
print RFILE "rect(chromSpaceStartBot, 0, chromSpaceEndBot, 1, col=\"white\", border=NA)\n";
print RFILE "polygon(karoBreakBot,karoBreakBotY,col=\"white\",border=NA)\n";
print RFILE "lines(karoBorderBotX,karoBorderBotY,col=\"black\")\n";

# bottom cnv
print RFILE "par(fig=c(0,1,0.28,0.42), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(0,$cnMax), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", main=\"\", xlab=\"\", ylab=\"Copy Number\\nVariation\",)\n";
print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
#print RFILE "points(cnPosBot, cnValBot, col=cnColBot, pch=20)\n";
print RFILE "segments(cnStartBot, cnStateBot, cnEndBot, cnStateBot, col=cnColBot, lwd=20)\n";
print RFILE "segments(cnStartBot, cnStateBot, cnEndBot, cnStateBot, col=\"black\", lwd=2)\n";

# bottom het freqs
print RFILE "par(fig=c(0,1,0.14,0.28), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(0, 1), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"Germline Het\\nFreq in Tumour\",)\n";
print RFILE "image(hetBotX, hetY, matrix(hetFreqsBot," . (int($chrBotOffset{end} / $imageWindow) + 1) . ",$imageY + 1), add=TRUE, col=c(grey.colors(5,start=1,end=0.5),grey.colors(35,start=0.5,end=0)))\n";
print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
#print RFILE "points(hetPosBot, hetFreqsBot, col=\"black\", pch=\".\")\n";

# bottom rainfall
print RFILE "par(fig=c(0,1,0.0,0.14), mar=c(0, 3, 0, 0) + 0.1, new=TRUE)\n";
print RFILE "plot(c(0,$chrBotOffset{end}), c(1,$rainMaxBot), cex.lab=0.75, cex.axis=0.5, type = \"n\", bty=\"n\", xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"Somatic SNV\\nRainfall\",)\n";
print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
print RFILE "points(xvalsBot, yvalsBot, col=colsBot, pch=20)\n";
print RFILE "axis(side=1, at=atLabelsBot, labels=chrLabelsBot, tick=FALSE)\n";
if (scalar @geneNameBot > 0)
{
	print RFILE "text(genePosBot,1, labels=geneNameBot, pos=4, srt=45, cex=0.75)\n";
}




print RFILE "mtext(\"$headerString\", NORTH<-3, line=1, adj=0, outer=TRUE, cex=2)\n";
print RFILE "mtext(\"$footerString\", SOUTH<-1, line=1, adj=0, outer=TRUE, cex=1.5)\n";

print RFILE "dev.off()\n";






my %chromTopPanels = (
	"chr1" => "c(0,249250621)",
	"chr2" => "c(249250621,492449994)",
	"chr3" => "c(492449994,690472424)",
	"chr4" => "c(690472424,881626700)",
	"chr5" => "c(881626700,1062541960)",
	"chr6" => "c(1062541960,1233657027)",
	"chr7-8" => "c(1233657027,1539159712)"
);

my %chromBotPanels = (
	"chr9-10" => "c(0,276748178)",
	"chr11-12" => "c(276748178,545606589)",
	"chr13-14" => "c(545606589,768126007)",
	"chr15-17" => "c(768126007,1042207362)",
	"chr18-22" => "c(1042207362,1341873574)",
	"chrX-Y" => "c(1341873574,1556517700)"
);

my $axisText = "axis(side=1, line=0.25,";

print RFILE "singleChromLabs = c(0, NA, 20, NA, 40, NA, 60, NA, 80, NA, 100, NA, 120, NA, 140, NA, 160, NA, 180, NA, 200, NA, 220, NA, 240, NA, 260, NA, 280, NA, 300, NA, 320, NA, 340, NA, 360)\n";

my %chromTopAxis = (
	"chr1" => "$axisText at=c(0:24*10000000 + 0), labels=c(singleChromLabs[1:25]))\n",
	"chr2" => "$axisText at=c(0:24*10000000 + 249250621), labels=c(singleChromLabs[1:25]))\n",
	"chr3" => "$axisText at=c(0:19*10000000 + 492449994), labels=c(singleChromLabs[1:20]))\n",
	"chr4" => "$axisText at=c(0:19*10000000 + 690472424), labels=c(singleChromLabs[1:20]))\n",
	"chr5" => "$axisText at=c(0:18*10000000 + 881626700), labels=c(singleChromLabs[1:19]))\n",
	"chr6" => "$axisText at=c(0:17*10000000 + 1062541960), labels=c(singleChromLabs[1:18]))\n",
	"chr7-8" => "$axisText at=c(0:15*10000000 + 1233657027), labels=c(singleChromLabs[1:16]))\n$axisText at=c(0:14*10000000 + 1392795690), labels=c(singleChromLabs[1:15]))\n",
);

my %chromBotAxis = (
	"chr9-10" => "$axisText at=c(0:14*10000000 + 0), labels=c(singleChromLabs[1:15]))\n$axisText at=c(0:14*10000000 + 141213431), labels=c(singleChromLabs[1:15]))\n",
	"chr11-12" => "$axisText at=c(0:13*10000000 + 276748178), labels=c(singleChromLabs[1:14]))\n$axisText at=c(0:13*10000000 + 411754694), labels=c(singleChromLabs[1:14]))\n",
	"chr13-14" => "$axisText at=c(0:11*10000000 + 545606589), labels=c(singleChromLabs[1:12]))\n$axisText at=c(0:10*10000000 + 660776467), labels=c(singleChromLabs[1:11]))\n",
	"chr15-17" => "$axisText at=c(0:10*10000000 + 768126007), labels=c(singleChromLabs[1:11]))\n$axisText at=c(0:9*10000000 + 870657399), labels=c(singleChromLabs[1:10]))\n$axisText at=c(0:8*10000000 + 961012152), labels=c(singleChromLabs[1:9]))\n",
	"chr18-22" => "$axisText at=c(0:7*10000000 + 1042207362), labels=c(singleChromLabs[1:8]))\n$axisText at=c(0:5*10000000 + 1120284610), labels=c(singleChromLabs[1:6]))\n$axisText at=c(0:6*10000000 + 1179413593), labels=c(singleChromLabs[1:7]))\n$axisText at=c(0:4*10000000 + 1242439113), labels=c(singleChromLabs[1:5]))\n$axisText at=c(0:5*10000000 + 1290569008), labels=c(singleChromLabs[1:6]))\n",
	"chrX-Y" => "$axisText at=c(0:15*10000000 + 1341873574), labels=c(singleChromLabs[1:16]))\n$axisText at=c(0:5*10000000 + 1497144134), labels=c(singleChromLabs[1:6]))\n",
);



for my $range (keys %chromTopPanels)
{
	print RFILE "png(filename = \"$outName.whole_$range.png\", width = 1000, height = 300)\n";

	# top rainfall
	print RFILE "par(fig=c(0,1,0.78,1), mar=c(0, 10, 0, 1) + 0.1)\n";
	print RFILE "plot($chromTopPanels{$range}, c(1,$rainMaxTop),  type = \"n\", bty=\"n\", axes=F, xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
	print RFILE "points(xvalsTop, yvalsTop, col=colsTop, pch = 20, cex=2)\n";
	print RFILE "axis(side=2, cex.axis=.8, line=0.25, las=1)\n";
	print RFILE "mtext(\"Somatic SNV\\nRainfall\", side=2, line=3.5, las=1)\n";
	if (scalar @geneNameTop > 0)
	{
		print RFILE "text(genePosTop,1, labels=geneNameTop, pos=4, srt=45, cex=0.8)\n";
	}
	
	# top het freqs
	print RFILE "par(fig=c(0,1,0.55,0.77), mar=c(0, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0, 1), axes=F, type = \"n\", bty=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "image(hetTopX, hetY, matrix(hetFreqsTop," . (int($chrTopOffset{end} / $imageWindow) + 1) . ",$imageY + 1), add=TRUE, col=c(grey.colors(5,start=1,end=0.5),grey.colors(35,start=0.5,end=0)))\n";
	print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
	print RFILE "axis(side=2, cex.axis=.8, line=0.25, las=1)\n";
	print RFILE "mtext(\"Germline Het\\nFreq in Tumour\", side=2, line=3.5, las=1)\n";
	#print RFILE "points(hetPosTop, hetFreqsTop, col=\"black\", pch=\".\")\n";
	
	# top cnv
	print RFILE "par(fig=c(0,1,0.32,0.54), mar=c(0, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0,$cnMax), type = \"n\", bty=\"n\", axes=F, main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
	#print RFILE "points(cnPosTop, cnValTop, col=cnColTop, pch=20)\n";
	print RFILE "segments(cnStartTop, cnStateTop, cnEndTop, cnStateTop, col=cnColTop, lwd=12)\n";
	print RFILE "segments(cnStartTop, cnStateTop, cnEndTop, cnStateTop, col=\"black\", lwd=2)\n";
	print RFILE "axis(side=2, cex.axis=.8, line=0.25, las=1)\n";
	print RFILE "mtext(\"Copy Number\\nVariation\", side=2, line=3.5, las=1)\n";
	
	# structural
	print RFILE "par(fig=c(0,1,0.23,0.32), mar=c(0, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0,1), cex.lab=1, cex.axis=0.8, type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "segments(topBreakX1,topBreakY1,topBreakX2,topBreakY2,col=topBreakCol,lwd=2)\n";
	print RFILE "mtext(\"Structural\\nBreak Points\", side=2, line=3.5, las=1)\n";

	#top karotype
	print RFILE "par(fig=c(0,1,0.0,0.23), mar=c(3.5, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0,1),type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\")\n";
	print RFILE "rect(karoStartTop, 0, karoEndTop, 1, col=karoColTop, border=NA)\n";
	print RFILE "rect(chromSpaceStartTopS, 0, chromSpaceEndTopS, 1, col=\"white\", border=NA)\n";
	print RFILE "polygon(karoBreakTop,karoBreakTopY,col=\"white\",border=NA)\n";
	print RFILE "lines(karoBorderTopXS,karoBorderTopYS,col=\"black\")\n";
	
	print RFILE "$chromTopAxis{$range}\n";
	print RFILE "axis(side=1, at=atLabelsTop, labels=chrLabelsTopS, tick=FALSE, line = 1.5, cex.axis=1.5)\n";

	print RFILE "dev.off()\n";

}


for my $range (keys %chromBotPanels)
{
	print RFILE "png(filename = \"$outName.whole_$range.png\", width = 1000, height = 300)\n";

	# top rainfall
	print RFILE "par(fig=c(0,1,0.78,1), mar=c(0, 10, 0, 1) + 0.1)\n";
	print RFILE "plot($chromBotPanels{$range}, c(1,$rainMaxBot),  type = \"n\", bty=\"n\", axes=F, xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
	print RFILE "points(xvalsBot, yvalsBot, col=colsBot, pch = 20, cex=2)\n";
	print RFILE "axis(side=2, cex.axis=.8, line=0.25, las=1)\n";
	print RFILE "mtext(\"Somatic SNV\\nRainfall\", side=2, line=3.5, las=1)\n";
	if (scalar @geneNameBot > 0)
	{
		print RFILE "text(genePosBot,1, labels=geneNameBot, pos=4, srt=45, cex=0.8)\n";
	}
	
	# top het freqs
	print RFILE "par(fig=c(0,1,0.55,0.77), mar=c(0, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0, 1), axes=F, type = \"n\", bty=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "image(hetBotX, hetY, matrix(hetFreqsBot," . (int($chrBotOffset{end} / $imageWindow) + 1) . ",$imageY + 1), add=TRUE, col=c(grey.colors(5,start=1,end=0.5),grey.colors(35,start=0.5,end=0)))\n";
	print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
	print RFILE "axis(side=2, cex.axis=.8, line=0.25, las=1)\n";
	print RFILE "mtext(\"Germline Het\\nFreq in Tumour\", side=2, line=3.5, las=1)\n";
	#print RFILE "points(hetPosBot, hetFreqsBot, col=\"black\", pch=\".\")\n";
	
	# top cnv
	print RFILE "par(fig=c(0,1,0.32,0.54), mar=c(0, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0,$cnMax), type = \"n\", bty=\"n\", axes=F, main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
	#print RFILE "points(cnPosBot, cnValBot, col=cnColBot, pch=20)\n";
	print RFILE "segments(cnStartBot, cnStateBot, cnEndBot, cnStateBot, col=cnColBot, lwd=12)\n";
	print RFILE "segments(cnStartBot, cnStateBot, cnEndBot, cnStateBot, col=\"black\", lwd=2)\n";
	print RFILE "axis(side=2, cex.axis=.8, line=0.25, las=1)\n";
	print RFILE "mtext(\"Copy Number\\nVariation\", side=2, line=3.5, las=1)\n";
	
	# structural
	print RFILE "par(fig=c(0,1,0.23,0.32), mar=c(0, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0,1), cex.lab=1, cex.axis=0.8, type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "segments(botBreakX1,botBreakY1,botBreakX2,botBreakY2,col=botBreakCol,lwd=2)\n";
	print RFILE "mtext(\"Structural\\nBreak Points\", side=2, line=3.5, las=1)\n";

	#top karotype
	print RFILE "par(fig=c(0,1,0.0,0.23), mar=c(3.5, 10, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0,1),type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\")\n";
	print RFILE "rect(karoStartBot, 0, karoEndBot, 1, col=karoColBot, border=NA)\n";
	print RFILE "rect(chromSpaceStartBotS, 0, chromSpaceEndBotS, 1, col=\"white\", border=NA)\n";
	print RFILE "polygon(karoBreakBot,karoBreakBotY,col=\"white\",border=NA)\n";
	print RFILE "lines(karoBorderBotXS,karoBorderBotYS,col=\"black\")\n";
	
	print RFILE "$chromBotAxis{$range}\n";
	print RFILE "axis(side=1, at=atLabelsBot, labels=chrLabelsBotS, tick=FALSE, line = 1.5, cex.axis=1.5)\n";

	print RFILE "dev.off()\n";

}




# single chroms with no axis

for my $range (keys %chromTopPanels)
{
	print RFILE "png(filename = \"$outName.whole_${range}_no-axis.png\", width = 1000, height = 300)\n";

	# top rainfall
	print RFILE "par(fig=c(0,1,0.78,1), mar=c(0, 1, 0, 1) + 0.1)\n";
	print RFILE "plot($chromTopPanels{$range}, c(1,$rainMaxTop),  type = \"n\", bty=\"n\", axes=F, xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
	print RFILE "points(xvalsTop, yvalsTop, col=colsTop, pch = 20, cex=2)\n";
	if (scalar @geneNameTop > 0)
	{
		print RFILE "text(genePosTop,1, labels=geneNameTop, pos=4, srt=45, cex=0.8)\n";
	}
	
	# top het freqs
	print RFILE "par(fig=c(0,1,0.55,0.77), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0, 1), axes=F, type = \"n\", bty=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "image(hetTopX, hetY, matrix(hetFreqsTop," . (int($chrTopOffset{end} / $imageWindow) + 1) . ",$imageY + 1), add=TRUE, col=c(grey.colors(5,start=1,end=0.5),grey.colors(35,start=0.5,end=0)))\n";
	print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
	#print RFILE "points(hetPosTop, hetFreqsTop, col=\"black\", pch=\".\")\n";
	
	# top cnv
	print RFILE "par(fig=c(0,1,0.32,0.54), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0,$cnMax), type = \"n\", bty=\"n\", axes=F, main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsTop, col=\"darkgrey\")\n";
	#print RFILE "points(cnPosTop, cnValTop, col=cnColTop, pch=20)\n";
	print RFILE "segments(cnStartTop, cnStateTop, cnEndTop, cnStateTop, col=cnColTop, lwd=12)\n";
	print RFILE "segments(cnStartTop, cnStateTop, cnEndTop, cnStateTop, col=\"black\", lwd=2)\n";
	
	# structural
	print RFILE "par(fig=c(0,1,0.23,0.32), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0,1), cex.lab=1, cex.axis=0.8, type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "segments(topBreakX1,topBreakY1,topBreakX2,topBreakY2,col=topBreakCol,lwd=2)\n";

	#top karotype
	print RFILE "par(fig=c(0,1,0.0,0.23), mar=c(3.5, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromTopPanels{$range}, c(0,1),type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\")\n";
	print RFILE "rect(karoStartTop, 0, karoEndTop, 1, col=karoColTop, border=NA)\n";
	print RFILE "rect(chromSpaceStartTopS, 0, chromSpaceEndTopS, 1, col=\"white\", border=NA)\n";
	print RFILE "polygon(karoBreakTop,karoBreakTopY,col=\"white\",border=NA)\n";
	print RFILE "lines(karoBorderTopXS,karoBorderTopYS,col=\"black\")\n";
	
	print RFILE "$chromTopAxis{$range}\n";
	print RFILE "axis(side=1, at=atLabelsTop, labels=chrLabelsTopS, tick=FALSE, line = 1.5, cex.axis=1.5)\n";

	print RFILE "dev.off()\n";

}


for my $range (keys %chromBotPanels)
{
	print RFILE "png(filename = \"$outName.whole_${range}_no-axis.png\", width = 1000, height = 300)\n";

	# top rainfall
	print RFILE "par(fig=c(0,1,0.78,1), mar=c(0, 1, 0, 1) + 0.1)\n";
	print RFILE "plot($chromBotPanels{$range}, c(1,$rainMaxBot),  type = \"n\", bty=\"n\", axes=F, xaxt=\"n\", log = \"y\", main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
	print RFILE "points(xvalsBot, yvalsBot, col=colsBot, pch = 20, cex=2)\n";
	if (scalar @geneNameBot > 0)
	{
		print RFILE "text(genePosBot,1, labels=geneNameBot, pos=4, srt=45, cex=0.8)\n";
	}
	
	# top het freqs
	print RFILE "par(fig=c(0,1,0.55,0.77), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0, 1), axes=F, type = \"n\", bty=\"n\", ylim=c(0, 1), main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "image(hetBotX, hetY, matrix(hetFreqsBot," . (int($chrBotOffset{end} / $imageWindow) + 1) . ",$imageY + 1), add=TRUE, col=c(grey.colors(5,start=1,end=0.5),grey.colors(35,start=0.5,end=0)))\n";
	print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
	#print RFILE "points(hetPosBot, hetFreqsBot, col=\"black\", pch=\".\")\n";
	
	# top cnv
	print RFILE "par(fig=c(0,1,0.32,0.54), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0,$cnMax), type = \"n\", bty=\"n\", axes=F, main=\"\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "abline(v=vertsBot, col=\"darkgrey\")\n";
	#print RFILE "points(cnPosBot, cnValBot, col=cnColBot, pch=20)\n";
	print RFILE "segments(cnStartBot, cnStateBot, cnEndBot, cnStateBot, col=cnColBot, lwd=12)\n";
	print RFILE "segments(cnStartBot, cnStateBot, cnEndBot, cnStateBot, col=\"black\", lwd=2)\n";
	
	# structural
	print RFILE "par(fig=c(0,1,0.23,0.32), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0,1), cex.lab=1, cex.axis=0.8, type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\", las=1)\n";
	print RFILE "segments(botBreakX1,botBreakY1,botBreakX2,botBreakY2,col=botBreakCol,lwd=2)\n";

	#top karotype
	print RFILE "par(fig=c(0,1,0.0,0.23), mar=c(3.5, 1, 0, 1) + 0.1, new=TRUE)\n";
	print RFILE "plot($chromBotPanels{$range}, c(0,1),type = \"n\", bty=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", xaxs=\"i\")\n";
	print RFILE "rect(karoStartBot, 0, karoEndBot, 1, col=karoColBot, border=NA)\n";
	print RFILE "rect(chromSpaceStartBotS, 0, chromSpaceEndBotS, 1, col=\"white\", border=NA)\n";
	print RFILE "polygon(karoBreakBot,karoBreakBotY,col=\"white\",border=NA)\n";
	print RFILE "lines(karoBorderBotXS,karoBorderBotYS,col=\"black\")\n";
	
	print RFILE "$chromBotAxis{$range}\n";
	print RFILE "axis(side=1, at=atLabelsBot, labels=chrLabelsBotS, tick=FALSE, line = 1.5, cex.axis=1.5)\n";

	print RFILE "dev.off()\n";

}
















close RFILE;
`Rscript $outputFile.R`;





sub mean
{
	my $sum = 0;
	my $count = 0;
	my $array = $_[0];
	for my $val (@{$array})
	{
		$sum += $val;
		$count++;
	}
	return $sum/$count;
}





# getBase returns the base at a specific position in the reference.  It will initialize fasta handles and pull new chunks of reference if necessary
# input is chromosome, position, reference hash and fasta handles
# output is a single base (and the reference hash and fasta handles may be modified)
sub getBase
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];
    my $fastaHandles = $_[3];

    my $chunkStart = int(($pos - 1) / $reference->{"chunkSize"}) * $reference->{"chunkSize"} + 1;       # +1 because the first base in the reference is 1, $pos - 1 so that multiples of chunk size resolve to the correct chunk
    my $chunkEnd = $chunkStart + $reference->{"chunkSize"} - 1;

    unless (exists $reference->{$chr}{$chunkStart}{$pos})       # if the position isn't in our hash, we need to get a new chunk from the reference
    {
        unless (exists ($fastaHandles->{$chr}))     # create a handle for the chromosome fasta, if it doesn't exist
        {
#           warn "Creating fasta handle for $chr\n";
            $fastaHandles->{$chr} = Bio::DB::Fasta->new("$fastaHandles{path}/$chr.fa");
        }

#       warn "Pulling $chr:$chunkStart-$chunkEnd from fasta\n";
        my $newChunk = uc($fastaHandles{$chr}->seq($chr, $chunkStart, $chunkEnd));
        my $i = $chunkStart;
        for my $base (split("", $newChunk))
        {
            $reference->{$chr}{$chunkStart}{$i} = $base;
            $i++;
        }
    }
#   warn "returning $reference->{$chr}{$chunkStart}{$pos}\n";
    return $reference->{$chr}{$chunkStart}{$pos};
}

# getRange returns a string of bases from the reference in the specified range by calling getBase
# input is chromosome, start pos, end pos, reference hash and fasta handles
# output is a string of bases
sub getRange
{
    my $chr = $_[0];
    my $start = $_[1];
    my $end = $_[2];
    my $reference = $_[3];
    my $fastaHandles = $_[4];

    my $seq = "";

    for (my $p = $start; $p <= $end; $p++)
    {
        $seq .= getBase($chr, $p, $reference, $fastaHandles);
#       warn "Got base: $chr:$p\t$seq\n";
    }

    return $seq;
}


# freeMemChunks tests if the next indel to be processed is on a different chromosome or more than a chunk away from the reference sequences currently in memory
#   if there exist chunks that we won't need again (assuming the input is sorted) the chunks will be removed from the reference hash
# input is the chromosome and position of the current indel, and the reference hash
# there is no output
sub freeMemChunks
{
    my $chr = $_[0];
    my $pos = $_[1];
    my $reference = $_[2];

    # delete chunks from non-current chromosomes
    for my $refChr (keys %$reference)
    {
        if (($refChr ne $chr) and ($refChr ne "chunkSize"))
        {
#           warn "deleting all chunks for $refChr.\n";
            delete $reference->{$refChr};
        }
    }

    # delete chunks if they are more than 1.5 chunks away from the current indel
    # 1.5 so that we are at least in the middle of the current chunk before dropping the previous one
    for my $chunkPos (keys %{ $reference->{$chr} })
    {
        if ($chunkPos < ($pos - (1.5 * $reference->{"chunkSize"})))
        {
#           warn "deleting $chr:$chunkPos chunk.\n";
            delete $reference->{$chr}{$chunkPos};
        }
    }

    return;
}


