#!/usr/bin/perl

use strict;
use warnings;
#use Data::Dumper;

my $dataFile = shift;
my $fileName = shift;

my %hasKRAS = (		# these samples were confirmed to be KRAS positive by visual inspection
	"PCSI_0171_Pa_P_526" => 1,
	"PCSI_0280_Pa_P_526" => 1,
	"PCSI_0300_Pa_P_526" => 1,
	"PCSI_0338_Pa_P_526" => 1,
);

my %hitRank = (
	"nonsense" => 1,
	"breakpoint" => 2,
	"missense" => 3,
	"splicing" => 4,
	"promoter" => 5,
	"loss" => 6
);

my %groupCol = (
	"Typical" => "\"#332288\"",
	"DSBR" => "\"#88CCEE\"",
	"MMR" => "\"#117733\"",
	"Sig8" => "\"#999933\"",
	"Sig17" => "\"#661100\"",
);

my $l;
my @f;
my @header;
my %data;
my %values;
my $sample;
my @samples;

my %chrOffset = (
"chr1" => 0,
"chr2" => 249250621,
"chr3" => 492449994,
"chr4" => 690472424,
"chr5" => 881626700,
"chr6" => 1062541960,
"chr7" => 1233657027,
"chr8" => 1392795690,
"chr9" => 1539159712,
"chr10" => 1680373143,
"chr11" => 1815907890,
"chr12" => 1950914406,
"chr13" => 2084766301,
"chr14" => 2199936179,
"chr15" => 2307285719,
"chr16" => 2409817111,
"chr17" => 2500171864,
"chr18" => 2581367074,
"chr19" => 2659444322,
"chr20" => 2718573305,
"chr21" => 2781598825,
"chr22" => 2829728720,
"chrX" => 2881033286,
"chrY" => 3036303846,
"end" => 3095677412,
);


open (FILE, $dataFile) or die "Couldn't open $dataFile\n";

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
			$data{$sample}{$header[$i]} = $f[$i];
		}
		push(@samples,$sample);
	}
}
close FILE;

for $sample (keys %hasKRAS)
{
	if (exists $data{$sample})
	{
		$data{$sample}{'KRAS:missense'}++;
	}
}

my @genes;
my %seenGene;

for my $h (@header)
{
	if ($h =~ /(.*):(.*)/)
	{
		unless (($1 eq "DEL") or ($1 eq "DUP") or ($1 eq "INV"))
		{
			unless (exists $seenGene{$1})
			{
				push (@genes, $1);
				$seenGene{$1}++;
			}
		}
	}
}


my $numGenes = scalar(@genes);
my $numSamples = scalar(@samples);

my $xPad = 0.1;
my $yPad = 0.11;

my %markPos = (
    "nonsense" => 0.2,
    "missense" => 0.3,
    "noncoding" => 0.4,
    "splicing" => 0.5,
    "split" => 0.6
);

my %markPt = (
    "nonsense" => 22,
    "missense" => 21,
    "noncoding" => 23,
    "splicing" => 24,
    "split" => 25,
	"loss" => 1
);

my $blockSplit = 0.7;

my (@xLeft, @xRight, @yTop, @yBottom, @yTopFirst, @yTopSecond, @yBottomFirst, @yBottomSecond,@yTopCopy,@yBottomCopy);

my (@senseCol, @copyCol);
my (@xMarks, %yMarks, %ptMarks);

my $misColour = "#117733BB";
my $nonColour = "#AA4499BB";
my $lossColour = "#33228899";
#my $lossColour = "#3C53A1BB";
#my $lossColour = "#4477AA";
my $noCol = "lightgrey";

my %oncoCol = (
	"nonsense" => $nonColour,
	"split" => $nonColour,
	"missense" => $misColour,
	"splicing" => $misColour,
	"noncoding" => $misColour,
	"loss" => $lossColour,
);

my %summaryCount;

my ($val,$r,$g,$blue);

my ($firstHitFound, $secondHitFound);
my (@firstHitCol, @secondHitCol);
my $firstHitSize = (1.1 - $yPad*2) / 3;
my $secondHitSize = (1.1 - $yPad*2) / 3;
my $copySize = (1 - $yPad*2) / 3;

my (@yMarksFirst, @yMarksSecond, @ptMarksFirst, @ptMarksSecond);

my (%firstHitCounts, %secondHitCounts);

my (@yTopSumFirst,@yBottomSumFirst,@yTopSumSecond,@yBottomSumSecond,@yTopSumCopy,@yBottomSumCopy);
my (@xLeftSumFirst,@xRightSumFirst,@xLeftSumSecond,@xRightSumSecond,@xLeftSumCopy,@xRightSumCopy);

my $cnMin = 0;
my $cnMax = 6;

for (my $i = 0; $i < $numSamples; $i++)
{

    for (my $j = 0; $j < $numGenes; $j++)
    {
        push(@xLeft, $j + (($j/($numGenes - 1)) * $xPad));
        push(@xRight, ($j + 1) - ((($numGenes - 1 - $j)/($numGenes - 1)) * $xPad));

		push(@yTopFirst, ($numSamples - $i) - $yPad*($i/($numSamples-1)));
		push(@yBottomFirst, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $firstHitSize);
		push(@yTopSecond, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $firstHitSize);
		push(@yBottomSecond, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $firstHitSize - $secondHitSize);
		push(@yTopCopy, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $firstHitSize - $secondHitSize - ($yPad / 2));
		push(@yBottomCopy, (($numSamples - $i) - 1) - $yPad*($i/($numSamples-1)) + $yPad);

		push(@xMarks, $j + 0.5);
		push(@yMarksFirst, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - ($firstHitSize/2));
		push(@yMarksSecond, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $firstHitSize - ($secondHitSize/2));

		$firstHitFound = 0;
		$secondHitFound = 0;


		for my $mut (qw/nonsense split missense splicing noncoding/)
		{
			if (exists $data{$samples[$i]}{"$genes[$j]:$mut"})
			{
				if ($data{$samples[$i]}{"$genes[$j]:$mut"} > 0)
				{
					unless ($secondHitFound == 1)
					{
						if ($data{$samples[$i]}{"$genes[$j]:$mut"} > 1)
						{
							if ($firstHitFound == 0)
							{
								$firstHitFound = 1;
								$secondHitFound = 1;

								$summaryCount{first}{$genes[$j]}{$mut}++;
								$summaryCount{second}{$genes[$j]}{$mut}++;
	
								push(@firstHitCol, $oncoCol{$mut});
								push(@secondHitCol, $oncoCol{$mut});
	
								push(@ptMarksFirst, $markPt{$mut});
								push(@ptMarksSecond, $markPt{$mut});
							}
							else
							{
								$secondHitFound = 1;
								$summaryCount{second}{$genes[$j]}{$mut}++;
								push(@secondHitCol, $oncoCol{$mut});
								push(@ptMarksSecond, $markPt{$mut});
							}
						}
						else
						{
							if ($firstHitFound == 0)
							{
								$firstHitFound = 1;
								$summaryCount{first}{$genes[$j]}{$mut}++;
								push(@firstHitCol, $oncoCol{$mut});
								push(@ptMarksFirst, $markPt{$mut});
							}
							else
							{
								$secondHitFound = 1;
								$summaryCount{second}{$genes[$j]}{$mut}++;
								push(@secondHitCol, $oncoCol{$mut});
								push(@ptMarksSecond, $markPt{$mut});
							}
						}
					}
				}
			}
		}



		unless ($secondHitFound == 1)
		{
			if ($data{$samples[$i]}{"$genes[$j]:cn"} <= 0.5)
			{
				if ($firstHitFound == 0)
				{
					$firstHitFound = 1;
					$secondHitFound = 1;

					$summaryCount{first}{$genes[$j]}{loss}++;
					$summaryCount{second}{$genes[$j]}{loss}++;

					push(@firstHitCol, $oncoCol{loss});
					push(@secondHitCol, $oncoCol{loss});

					push(@ptMarksFirst, $markPt{loss});
					push(@ptMarksSecond, $markPt{loss});
				}
				else
				{
					$secondHitFound = 1;
					$summaryCount{second}{$genes[$j]}{loss}++;
					push(@secondHitCol, $oncoCol{loss});
					push(@ptMarksSecond, $markPt{loss});
				}
			}
			elsif($data{$samples[$i]}{"$genes[$j]:cn"} <= 1.5)
			{
				if ($firstHitFound == 0)
				{
					$firstHitFound = 1;
					$summaryCount{first}{$genes[$j]}{loss}++;
					push(@firstHitCol, $oncoCol{loss});
					push(@ptMarksFirst, $markPt{loss});
				}
				else
				{
					$secondHitFound = 1;
					$summaryCount{second}{$genes[$j]}{loss}++;
					push(@secondHitCol, $oncoCol{loss});
					push(@ptMarksSecond, $markPt{loss});
				}
			}
		}


		if ($firstHitFound == 0)
		{
			push(@firstHitCol, $noCol);
			push(@ptMarksFirst, "NA");
		}
		if ($secondHitFound == 0)
		{
			push(@secondHitCol, $noCol);
			push(@ptMarksSecond, "NA");
		}

		$firstHitCounts{$genes[$j]} += $firstHitFound;
		$secondHitCounts{$genes[$j]}+= $secondHitFound;


        if ($data{$samples[$i]}{"$genes[$j]:cn"} <= $cnMin)
        {
			$data{$samples[$i]}{"$genes[$j]:cn"} = $cnMin;
        }
        elsif ($data{$samples[$i]}{"$genes[$j]:cn"} >= $cnMax)
        {
			$data{$samples[$i]}{"$genes[$j]:cn"} = $cnMax;
        }

        $val = ($data{$samples[$i]}{"$genes[$j]:cn"} - $cnMin) / ($cnMax - $cnMin);
        $r = sprintf("%02x",255 * ( 0.237 - 2.13*$val + 26.92*($val**2) - 65.5*($val**3) + 63.5*($val**4) - 22.36*($val**5)));
        $g = sprintf("%02x",255 * ( ( (0.572 + 1.524*$val - 1.811*($val**2)) / (1 - 0.291*$val + 0.1574*($val**2)) )**2));
        $blue = sprintf("%02x",255 * ( 1/(1.579 - 4.03*$val + 12.92*($val**2) - 31.4*($val**3) + 48.6*($val**4) - 23.36*($val**5))));

        push(@copyCol, "#${r}${g}${blue}BB");

		push(@{ $summaryCount{cn_val}{$genes[$j]} }, $data{$samples[$i]}{"$genes[$j]:cn"});
		push(@{ $summaryCount{cn_col}{$genes[$j]} }, "#${r}${g}${blue}BB");

	}
}

my $count;
for my $g (keys %firstHitCounts)
{
	$count = $firstHitCounts{$g};
	$firstHitCounts{$g} = sprintf("%.0f", $count / $numSamples * 100);
}
for my $g (keys %secondHitCounts)
{
	$count = $secondHitCounts{$g};
	$secondHitCounts{$g} = sprintf("%.0f", $count / $numSamples * 100);
}







my $rfile;
open ($rfile, ">$fileName-oncoPrint.Rcode") or die "Couldn't open $fileName-oncoPrint.Rcode\n";

# trim underscore from donor names
my @trimDonors;
my $trim;

for my $donor (@samples)
{
	$trim = $donor;
#	$trim =~ s/^...._//;
#	$trim =~ s/^(....).*/$1/;
	push(@trimDonors, $trim);
}
printQuotedArrayToR($rfile,"donors",\@trimDonors);

my @donorsAt;
for (my $i = 0; $i < scalar(@samples); $i++)
{
	push(@donorsAt, $i + $i * 0.2 + 0.7);
}
printArrayToR($rfile,"donorsAt",\@donorsAt);



my %valVectors;
my %colVectors;

my %minVal = (
	"snv_count" => 0,
	"indel_count" => 0,
	"sv_count" => 0,
	"neo_antigens" => 0,
	"snv_ca" => 0,
	"snv_cg" => 0,
	"snv_ct" => 0,
	"snv_ta" => 0,
	"snv_tc" => 0,
	"snv_tg" => 0,
	"del_1_count" => 0,
	"del_4_count" => 0,
	"ins_1_count" => 0,
	"ins_4_count" => 0,
	"sv_del_count" => 0,
	"sv_dup_count" => 0,
	"sv_inv_count" => 0,
	"sv_tra_count" => 0,
);

my %maxVal = (
	"snv_count" => 20000,
	"indel_count" => 2000,
	"sv_count" => 200,
	"neo_antigens" => 100,
	"snv_ca" => 1,
	"snv_cg" => 1,
	"snv_ct" => 1,
	"snv_ta" => 1,
	"snv_tc" => 1,
	"snv_tg" => 1,
	"del_1_count" => 1,
	"del_4_count" => 1,
	"ins_1_count" => 1,
	"ins_4_count" => 1,
	"sv_del_count" => 1,
	"sv_dup_count" => 1,
	"sv_inv_count" => 1,
	"sv_tra_count" => 1,
);

my $type;

for $type (qw/snv_count indel_count sv_count snv_ca snv_cg snv_ct snv_ta snv_tc snv_tg del_1_count del_4_count ins_1_count ins_4_count sv_del_count sv_dup_count sv_inv_count sv_tra_count/)
{
	@{ $valVectors{$type} } = valueScaleVector(\%data, \@samples, $type, $minVal{$type}, $maxVal{$type});
}

$type = "neo_antigens";
@{ $valVectors{$type} } = onesVector(\%data, \@samples, $type);
@{ $colVectors{$type} } = colourScaleVector(\%data, \@samples, $type, $minVal{$type}, $maxVal{$type});

$type = "waddell";
my %waddellCols = (
	"stable" => "#94BE56",
	"locally rearranged" => "#F7C619",
	"scattered" => "#B09AC8",
	"unstable" => "#8F342B",
);
@{ $valVectors{$type} } = onesVector(\%data, \@samples, $type);
@{ $colVectors{$type} } = colourHashVector(\%data, \@samples, $type, \%waddellCols);

$type = "moffitt";
my %moffittCols = (
	"basal-like" => "orange",
	"unknown" => "grey",
	"classic" => "purple",
);
@{ $valVectors{$type} } = onesVector(\%data, \@samples, $type);
@{ $colVectors{$type} } = colourHashVector(\%data, \@samples, $type, \%moffittCols);

$type = "collisson";
my %collissonCols = (
	"classical" => "purple",
	"exocrine-like" => "yellow",
	"quasi-mesenchymal" => "navyblue",
);
@{ $valVectors{$type} } = onesVector(\%data, \@samples, $type);
@{ $colVectors{$type} } = colourHashVector(\%data, \@samples, $type, \%collissonCols);

for my $type (qw/snv_count indel_count sv_count snv_ca snv_cg snv_ct snv_ta snv_tc snv_tg del_1_count del_4_count ins_1_count ins_4_count sv_del_count sv_dup_count sv_inv_count sv_tra_count/)
{
	printArrayToR($rfile,$type,\@{ $valVectors{$type} });
}

for my $type (qw/neo_antigens waddell moffitt collisson/)
{
	printArrayToR($rfile,$type,\@{ $valVectors{$type} });
	printQuotedArrayToR($rfile,"${type}_col",\@{ $colVectors{$type} });
}






# onco vals here
# print rect coords

printArrayToR($rfile,"xLeft",\@xLeft);
printArrayToR($rfile,"xRight",\@xRight);

printArrayToR($rfile,"yTopFirst",\@yTopFirst);
printArrayToR($rfile,"yBottomFirst",\@yBottomFirst);

printArrayToR($rfile,"yTopSecond",\@yTopSecond);
printArrayToR($rfile,"yBottomSecond",\@yBottomSecond);

printArrayToR($rfile,"yTopCopy",\@yTopCopy);
printArrayToR($rfile,"yBottomCopy",\@yBottomCopy);

# print colours
printQuotedArrayToR($rfile,"colFirst",\@firstHitCol);
printQuotedArrayToR($rfile,"colSecond",\@secondHitCol);
printQuotedArrayToR($rfile,"colCopy",\@copyCol);

# print labels
my @geneLabels;
for (my $i = 0; $i < scalar(@genes); $i++)
{
#	push(@geneLabels, "$genes[$i]\\n$firstHitCounts{$genes[$i]}% ($secondHitCounts{$genes[$i]}%)");
	push(@geneLabels, "$genes[$i]");
}
printQuotedArrayToR($rfile,"labGene",\@geneLabels);


# print dots
printArrayToR($rfile,"xMarks",\@xMarks);

printArrayToR($rfile,"yMarksFirst",\@yMarksFirst);
printArrayToR($rfile,"ptMarksFirst",\@ptMarksFirst);

printArrayToR($rfile,"yMarksSecond",\@yMarksSecond);
printArrayToR($rfile,"ptMarksSecond",\@ptMarksSecond);




# collect cn segments

my %cnCols = (
	"0 copy" => "#3D52A1",
	"1 copy" => "#77B7E5",
	"2 copy" => "#E6F5FE",
	"3-4 copy" => "#FFE3AA",
	"5-6 copy" => "#F9BD7E",
	"7-8 copy" => "#ED875E",
	"9-10 copy" => "#D24D3E",
	">10 copy" => "#AE1C3E",
);
my $svBarSize = (1.1 - $yPad*2) / 3;
my $samp;
my $lines;
my ($chr, $start, $end, $cn, $stripTop, $stripBot);
my (@cnLeft, @cnRight, @cnTop, @cnBot, @cnCol);
my (@chrTop, @chrBot, @chrPos);
for (my $i = 0; $i < $numSamples; $i++)
{
	$samp = $samples[$i];

	unless ($data{$samp}{seg_file} eq "NA")
	{
		$lines = `zcat $data{$samp}{seg_file}`;
		chomp $lines;

		$stripTop = ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $svBarSize - ($yPad / 2);
		$stripBot = (($numSamples - $i) - 1) - $yPad*($i/($numSamples-1)) + $yPad;
		for $l (split(/\n/, $lines))
		{
			unless ($l =~ /^#chrom/)
			{
				@f = split(/\t/, $l);
				$chr = $f[0];
				$start = $f[1];
				$end = $f[2];
				$cn = $f[11];

				push(@cnLeft, int($start + $chrOffset{$chr}) / 1000);
				push(@cnRight, int($end + $chrOffset{$chr}) / 1000);

				push(@cnTop, $stripTop);
				push(@cnBot, $stripBot);

				if ($cn < 0.5)
				{
					push(@cnCol, $cnCols{"0 copy"});
				}
				elsif ($cn < 1.5)
				{
					push(@cnCol, $cnCols{"1 copy"});
				}
				elsif ($cn < 2.5)
				{
					push(@cnCol, $cnCols{"2 copy"});
				}
				elsif ($cn < 4.5)
				{
					push(@cnCol, $cnCols{"3-4 copy"});
				}
				elsif ($cn < 6.5)
				{
					push(@cnCol, $cnCols{"5-6 copy"});
				}
				elsif ($cn < 8.5)
				{
					push(@cnCol, $cnCols{"7-8 copy"});
				}
				elsif ($cn < 10.5)
				{
					push(@cnCol, $cnCols{"9-10 copy"});
				}
				else
				{
					push(@cnCol, $cnCols{">10 copy"});
				}

			}
		}

		# push lines for chromosome dividers
		push(@chrTop, $stripTop);
		push(@chrBot, $stripBot);
		push(@chrPos, 0);

		for $chr (qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY end/)
		{
			push(@chrTop, $stripTop);
			push(@chrBot, $stripBot);
			push(@chrPos, int($chrOffset{$chr} / 1000));
		}
	}
}


printArrayToR($rfile,"cnLeft",\@cnLeft);
printArrayToR($rfile,"cnRight",\@cnRight);
printArrayToR($rfile,"cnTop",\@cnTop);
printArrayToR($rfile,"cnBot",\@cnBot);
printQuotedArrayToR($rfile,"cnCol",\@cnCol);

printArrayToR($rfile,"chrTop",\@chrTop);
printArrayToR($rfile,"chrBot",\@chrBot);
printArrayToR($rfile,"chrPos",\@chrPos);



# bin svs for cn strip

my %svBins;
my $binSize = 25000000;

my $svBinMax = 20;
my $svBinHeight;
my $pos;
my (@binLeft, @binRight, @binTop, @binBot);

for (my $i = 0; $i < $numSamples; $i++)
{
	$samp = $samples[$i];
	%svBins = ();


	unless ($data{$samp}{sv_file} eq "NA")
	{
		open FILE, $data{$samp}{sv_file} or die "Couldn't open $data{$samp}{sv_file}\n";
	
		$l = <FILE>; # skip header
		while ($l = <FILE>)
		{
			@f = split(/\t/, $l);

			$chr = $f[0];
			$pos = $f[1];
			$pos += $chrOffset{$chr};
			$pos = int($pos/$binSize);
			$svBins{$pos}++;

			$chr = $f[2];
			$pos = $f[3];
			$pos += $chrOffset{$chr};
			$pos = int($pos/$binSize);
			$svBins{$pos}++;
		}
	}
	close FILE;

	for (my $bin = 0; $bin < int($chrOffset{end}/$binSize); $bin++)
	{
		if (exists $svBins{$bin})
		{
			push(@binLeft, ($bin * $binSize) / 1000);
			push(@binRight, (($bin * $binSize) + $binSize) / 1000);

			push(@binBot, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $svBarSize);

			if ($svBins{$bin} > $svBinMax)
			{
				$svBins{$bin} = $svBinMax;
			}

			$svBinHeight = $svBarSize - (($svBins{$bin} / $svBinMax) * $svBarSize);

			push(@binTop, ($numSamples - $i) - $yPad*($i/($numSamples-1)) - $svBinHeight);
		}
	}

}


printArrayToR($rfile,"svBinLeft",\@binLeft);
printArrayToR($rfile,"svBinRight",\@binRight);
printArrayToR($rfile,"svBinTop",\@binTop);
printArrayToR($rfile,"svBinBot",\@binBot);


# chrom labeles for cnv strip
my %chrLabPos = (
"1" => 124625310.5,
"2" => 370850307.5,
"3" => 591461209,
"4" => 786049562,
"5" => 972084330,
"6" => 1148099493.5,
"7" => 1313226358.5,
"8" => 1465977701,
"9" => 1609766427.5,
"10" => 1748140516.5,
"11" => 1883411148,
"12" => 2017840353.5,
"13" => 2142351240,
"14" => 2253610949,
"15" => 2358551415,
"16" => 2454994487.5,
"17" => 2540769469,
"18" => 2620405698,
"19" => 2689008813.5,
"20" => 2750086065,
"21" => 2805663772.5,
"22" => 2855381003,
"X" => 2958668566,
"Y" => 3065990629,
);

my @stripLab;
my @stripLabPos;

for $chr (sort keys %chrLabPos)
{
	push(@stripLab, $chr);
	push(@stripLabPos, int($chrLabPos{$chr} / 1000));
}

printQuotedArrayToR($rfile,"stripLab",\@stripLab);
printArrayToR($rfile,"stripLabPos",\@stripLabPos);

# main plot


# sort order
print $rfile "sort_order <- rev(c(1:$numSamples))\n";

# 

my $oncoLeft = 0;
my $oncoRight = 0.15;
my $snvLeft = 0.155;
my $snvRight = 0.255;
my $indelLeft = 0.255;
my $indelRight = 0.355;
my $svLeft = 0.355;
my $svRight = 0.455;
my $cnLeft = 0.46;
my $cnRight = 0.94;
my $neoLeft = 0.945;
my $neoRight = 0.955;
my $wadLeft = 0.955;
my $wadRight = 0.965;
my $mofLeft = 0.965;
my $mofRight = 0.975;
my $colLeft = 0.975;
my $colRight = 0.985;

my $axisCex = 0.8;
my $geneCex = 0.6;
my $barLabelOffset = 1;
my $topMar = 3.5;
my $botMar = 1;

#print $rfile "png(filename = \"$fileName-oncoPrint.png\", width = 700, height = 1000)\n";
#print $rfile "png(filename = \"$fileName-oncoPrint.png\", width = 1000, height = 150 + 50*$numSamples)\n";
print $rfile "png(filename = \"$fileName-oncoPrint.png\", width = 1800, height = 150 + 50*$numSamples)\n";

print $rfile "par(cex=2)\n";

# onco
print $rfile "par(fig=c($oncoLeft,$oncoRight,0,1), mar=c($botMar, 4.5, $topMar,0) + 0.1)\n";
print $rfile "plot(c(0,$numGenes), c(0,$numSamples), type=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", bty=\"n\",xaxs=\"i\",yaxs=\"i\")\n";
print $rfile "rect(xLeft, yBottomFirst, xRight, yTopFirst, col=c(t(matrix(colFirst,$numSamples,$numGenes,byrow=TRUE)[rev(sort_order),1:$numGenes])), border=\"NA\")\n";
print $rfile "rect(xLeft, yBottomSecond, xRight, yTopSecond, col=c(t(matrix(colSecond,$numSamples,$numGenes,byrow=TRUE)[rev(sort_order),1:$numGenes])), border=\"NA\")\n";
print $rfile "rect(xLeft, yBottomCopy, xRight, yTopCopy, col=c(t(matrix(colCopy,$numSamples,$numGenes,byrow=TRUE)[rev(sort_order),1:$numGenes])), border=\"NA\")\n";
#print $rfile "axis(2, at=($numGenes:1)-.5, labels=labGene, las=2, tick=FALSE, cex.axis=1)\n";
print $rfile "axis(2, at=($numSamples:1)-.5, labels=rev(donors[sort_order]), las=2, tick=FALSE, cex.axis=$axisCex,line=-0.8)\n";
print $rfile "axis(3, at=($numGenes:1)-.5, labels=rev(labGene), las=2, tick=FALSE, cex.axis=$geneCex,line=-0.8)\n";
print $rfile "points(xMarks,c(t(matrix(yMarksFirst,$numSamples,$numGenes,byrow=TRUE)[sort_order,1:$numGenes])), pch=c(t(matrix(ptMarksFirst,$numSamples,$numGenes,byrow=TRUE)[sort_order,1:$numGenes])), cex=0.5, bg=\"black\")\n";
print $rfile "points(xMarks,c(t(matrix(yMarksSecond,$numSamples,$numGenes,byrow=TRUE)[sort_order,1:$numGenes])), pch=c(t(matrix(ptMarksSecond,$numSamples,$numGenes,byrow=TRUE)[sort_order,1:$numGenes])), cex=0.5, bg=\"black\")\n";

print $rfile "zeros = rep(\"NA\",$numSamples)\n";
# snv count
print $rfile "par(fig=c($snvLeft,$snvRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
my $snvSubCol = "\"#E8601C\",\"#F6C141\",\"#CAE0AB\",\"#4EB265\",\"#5289C7\",\"#B178A6\"";
print $rfile "barplot(matrix(rbind(snv_ca[sort_order],snv_cg[sort_order],snv_ct[sort_order],snv_ta[sort_order],snv_tc[sort_order],snv_tg[sort_order],zeros,zeros,zeros,zeros,zeros,zeros),ncol=$numSamples*2), col=c($snvSubCol), ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\",horiz=TRUE,width=c(.7,.3),space=c(0.4,0.2), border=\"NA\")\n";
print $rfile "barplot(matrix(rbind(zeros,zeros,zeros,zeros,zeros,zeros,snv_count[sort_order],zeros,zeros,zeros,zeros,zeros),ncol=$numSamples*2), col=c(\"black\"), ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\",horiz=TRUE,width=c(.7,.3),add=T,space=c(0.4,0.2), border=\"NA\")\n";
print $rfile "axis(3, at=c(0,0.5,1), labels=c(\"\",\"\",\"\"),line=0.2)\n";
print $rfile "mtext(\"SNVs\",side=3, line=1)\n";

# indel count
print $rfile "par(fig=c($indelLeft,$indelRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
my $indelSubCol = "\"#DC050C\",\"#F7EE55\",\"#7BAFDE\",\"#D6C1DE\"";
print $rfile "barplot(matrix(rbind(del_1_count[sort_order],del_4_count[sort_order],ins_1_count[sort_order],ins_4_count[sort_order],zeros,zeros,zeros,zeros),ncol=$numSamples*2), col=c($indelSubCol), ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\",horiz=TRUE,width=c(.7,.3),space=c(0.4,0.2), border=\"NA\")\n";
print $rfile "barplot(matrix(rbind(zeros,zeros,zeros,zeros,indel_count[sort_order],zeros,zeros,zeros),ncol=$numSamples*2), col=c(\"black\"), ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\",horiz=TRUE,width=c(.7,.3),add=T,space=c(0.4,0.2), border=\"NA\")\n";
print $rfile "axis(3, at=c(0,0.5,1), labels=c(\"\",\"\",\"\"),line=0.2)\n";
print $rfile "mtext(\"Indels\",side=3,line=1)\n";


# sv count
print $rfile "par(fig=c($svLeft,$svRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
my $svSubCol = "\"#F1932D\",\"#90C987\",\"#1965B0\",\"#882E72\"";
print $rfile "barplot(matrix(rbind(sv_del_count[sort_order],sv_dup_count[sort_order],sv_inv_count[sort_order],sv_tra_count[sort_order],zeros,zeros,zeros,zeros),ncol=$numSamples*2), col=c($svSubCol), ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\",horiz=TRUE,width=c(.7,.3),space=c(0.4,0.2), border=\"NA\")\n";
print $rfile "barplot(matrix(rbind(zeros,zeros,zeros,zeros,sv_count[sort_order],zeros,zeros,zeros),ncol=$numSamples*2), col=c(\"black\"), ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\",horiz=TRUE,width=c(.7,.3),add=T,space=c(0.4,0.2), border=\"NA\")\n";
print $rfile "axis(3, at=c(0,0.5,1), labels=c(\"\",\"\",\"\"),line=0.2)\n";
print $rfile "mtext(\"SVs\",side=3,line=1)\n";


# neo-antigen count
print $rfile "par(fig=c($neoLeft,$neoRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
print $rfile "barplot(neo_antigens[sort_order], col=neo_antigens_col[sort_order], ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\", horiz=TRUE)\n";
print $rfile "axis(3, at=.5, labels=\"Neo-antigens\", las=2, tick=FALSE, cex.axis=$geneCex,line=-0.8)\n";

# waddell
print $rfile "par(fig=c($wadLeft,$wadRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
print $rfile "barplot(waddell[sort_order], col=waddell_col[sort_order], ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\", horiz=TRUE)\n";
print $rfile "axis(3, at=.5, labels=\"Waddell Class\", las=2, tick=FALSE, cex.axis=$geneCex,line=-0.8)\n";

# moffitt
print $rfile "par(fig=c($mofLeft,$mofRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
print $rfile "barplot(moffitt[sort_order], col=moffitt_col[sort_order], ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\", horiz=TRUE)\n";
print $rfile "axis(3, at=.5, labels=\"Moffitt Class\", las=2, tick=FALSE, cex.axis=$geneCex,line=-0.8)\n";

# collisson
print $rfile "par(fig=c($colLeft,$colRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
print $rfile "barplot(collisson[sort_order], col=collisson_col[sort_order], ylab=\"\", axes=FALSE, xlim=c(0,1),yaxs=\"i\", horiz=TRUE)\n";
print $rfile "axis(3, at=.5, labels=\"Collisson Class\", las=2, tick=FALSE, cex.axis=$geneCex,line=-0.8)\n";


# copy number strip
print $rfile "par(fig=c($cnLeft,$cnRight,0,1), mar=c($botMar, 0, $topMar, 0) + 0.1, new=TRUE)\n";
print $rfile "plot(c(0,$chrOffset{end}/1000), c(0,$numSamples), type=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", bty=\"n\",xaxs=\"i\",yaxs=\"i\")\n";
print $rfile "rect(cnLeft, cnBot, cnRight, cnTop, col=cnCol, border=\"NA\")\n";
print $rfile "segments(chrPos,chrTop,chrPos,chrBot,col=\"black\")\n";

print $rfile "rect(svBinLeft, svBinBot, svBinRight, svBinTop, col=\"black\", border=\"NA\")\n";

print $rfile "mtext(\"CNV and SV across the Genome\",side=3,line=1)\n";
print $rfile "axis(3, at=stripLabPos, labels=stripLab, las=2, tick=FALSE, cex.axis=$geneCex,line=-0.8)\n";

# summary plots



# snv count
#print $rfile "par(fig=c($snvBarsSum), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
#print $rfile "barplot(snv_count[rev(order(snv_count))], col=snv_count_col[rev(order(snv_count))], ylab=\"\", cex.lab=0.5, cex.axis=0.5, axes=FALSE, border=NA, space=0, ylim=c(0,$maxCounts{snv_count}),xaxs=\"i\")\n";

# indel count
#print $rfile "par(fig=c($indelBarsSum), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
#print $rfile "barplot(indel_count[rev(order(indel_count))], col=indel_count_col[rev(order(indel_count))], ylab=\"\", cex.lab=0.5, cex.axis=0.5, axes=FALSE, border=NA, space=0, ylim=c(0,$maxCounts{indel_count}),xaxs=\"i\")\n";

# sv count
#print $rfile "par(fig=c($svBarsSum), mar=c(0, 1, 0, 1) + 0.1, new=TRUE)\n";
#print $rfile "barplot(sv_count[rev(order(sv_count))], col=sv_count_col[rev(order(sv_count))], ylab=\"\", cex.lab=0.5, cex.axis=0.5, axes=FALSE, border=NA, space=0, ylim=c(0,$maxCounts{sv_count}),xaxs=\"i\")\n";

# onco side
#print $rfile "par(fig=c($oncoPlotSum), mar=c(0, 1, 3, 1) + 0.1, new=TRUE)\n";
#print $rfile "plot(c(0,1), c(0,$numGenes), type=\"n\", xaxt=\"n\", yaxt=\"n\", xlab=\"\", ylab=\"\", bty=\"n\",xaxs=\"i\",yaxs=\"i\")\n";
#print $rfile "rect(x_left_sum_first, y_bottom_sum_first, x_right_sum_first, y_top_sum_first, col=sum_mut_col, border=\"NA\")\n";
#print $rfile "rect(x_left_sum_second, y_bottom_sum_second, x_right_sum_second, y_top_sum_second, col=sum_mut_col, border=\"NA\")\n";
#print $rfile "rect(x_left_sum_copy, y_bottom_sum_copy, x_right_sum_copy, y_top_sum_copy, col=sum_copy_col, border=\"NA\")\n";
#print $rfile "abline(v=c(0,.25,.5,.75,1),lty=2)\n";
#print $rfile "axis(3,at=c(0,.25,.5,.75,1),cex.axis=0.6)\n";





close $rfile;
`Rscript $fileName-oncoPrint.Rcode`;



sub valueVector
{
	my $dataRef = shift;
	my $sampRef = shift;
	my $type = shift;
	my $min = shift;
	my $max = shift;

	my @vector;

	for my $samp (@{ $sampRef })
	{
		if (exists $dataRef->{$samp}{$type})
		{
			unless ($dataRef->{$samp}{$type} eq "NA")
			{
				if ($dataRef->{$samp}{$type} >= $max)
				{
					push(@vector, $max);
				}
				elsif ($dataRef->{$samp}{$type} <= $min)
				{
					push(@vector, $min);
				}
				else
				{
					push(@vector,$dataRef->{$samp}{$type});
				}
			}
			else
			{
				push(@vector, "NA");
			}
		}
		else
		{
			push(@vector, "NA");
		}
	}

	return @vector;
}

sub valueScaleVector
{
	my $dataRef = shift;
	my $sampRef = shift;
	my $type = shift;
	my $min = shift;
	my $max = shift;

	my @vector;

	for my $samp (@{ $sampRef })
	{
		if (exists $dataRef->{$samp}{$type})
		{
			unless ($dataRef->{$samp}{$type} eq "NA")
			{
				if ($dataRef->{$samp}{$type} >= $max)
				{
					push(@vector, ($max - $min) / ($max - $min));
				}
				elsif ($dataRef->{$samp}{$type} <= $min)
				{
					push(@vector, ($min - $min) / ($max - $min));
				}
				else
				{
					push(@vector,($dataRef->{$samp}{$type} - $min) / ($max - $min));
				}
			}
			else
			{
				push(@vector, "NA");
			}
		}
		else
		{
			push(@vector, "NA");
		}
	}

	return @vector;
}

sub colourScaleVector
{
	my $dataRef = shift;
	my $sampRef = shift;
	my $type = shift;
	my $min = shift;
	my $max = shift;

	my @vector;

	for my $samp (@{ $sampRef })
	{
		if (exists $dataRef->{$samp}{$type})
		{
			unless ($dataRef->{$samp}{$type} eq "NA")
			{
				if ($dataRef->{$samp}{$type} >= $max)
				{
					push(@vector, valToCol(($max - $min) / ($max - $min)));
				}
				elsif ($dataRef->{$samp}{$type} <= $min)
				{
					push(@vector, valToCol(($min - $min) / ($min - $min)));
				}
				else
				{
					push(@vector, valToCol(($dataRef->{$samp}{$type} - $min) / ($max - $min)));
				}
			}
			else
			{
				push(@vector, "NA");
			}
		}
		else
		{
			push(@vector, "NA");
		}
	}

	return @vector;
}


sub onesVector
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
				push(@vector, 1);
			}
			else
			{
				push(@vector, "NA");
			}
		}
		else
		{
			push(@vector, "NA");
		}
	}

	return @vector;
}


sub colourHashVector
{
	my $dataRef = shift;
	my $sampRef = shift;
	my $type = shift;
	my $colRef = shift;

	my @vector;

	for my $samp (@{ $sampRef })
	{
		if (exists $dataRef->{$samp}{$type})
		{
			unless ($dataRef->{$samp}{$type} eq "NA")
			{
				push(@vector, $colRef->{$dataRef->{$samp}{$type}});
			}
			else
			{
				push(@vector, "NA");
			}
		}
		else
		{
			push(@vector, "NA");
		}
	}

	return @vector;
}

sub valToCol
{
	my $val = shift;

    my $r = sprintf("%02x",255 * ( 0.237 - 2.13*$val + 26.92*($val**2) - 65.5*($val**3) + 63.5*($val**4) - 22.36*($val**5)));
    my $g = sprintf("%02x",255 * ( ( (0.572 + 1.524*$val - 1.811*($val**2)) / (1 - 0.291*$val + 0.1574*($val**2)) )**2));
    my $blue = sprintf("%02x",255 * ( 1/(1.579 - 4.03*$val + 12.92*($val**2) - 31.4*($val**3) + 48.6*($val**4) - 23.36*($val**5))));

    return "#${r}${g}${blue}";
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







