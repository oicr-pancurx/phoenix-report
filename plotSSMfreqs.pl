#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $dataFile = shift;
my $fileName = shift;

my $plotType = "ssm_freq";

my $l;
my @header;
my @f;
my %data;
my %values;
my $sample;
my @samples;

open (FILE, $dataFile) or die "Couldn't open $dataFile\n";

$l = <FILE>;
@header = split(/,/, $l);
$l = <FILE>;
my $i = 0;
for my $val (split(/,/, $l))
{
	$data{$header[$i]} = $val;
	$i++;
}
close FILE;

my %colours = (
    "del_1" => "#DC050C",
    "del_4" => "#F7EE55",
    "ins_1" => "#7BAFDE",
    "ins_4" => "#D6C1DE",
    "CA" => "#E8601C",
    "CG" => "#F6C141",
    "CT" => "#CAE0AB",
    "TA" => "#4EB265",
    "TC" => "#5289C7",
    "TG" => "#B178A6",
);

my %comp = (
	"AC" => "TG",
	"AG" => "TC",
	"AT" => "TA",
	"GA" => "CT",
	"GC" => "CG",
	"GT" => "CA"
);

my ($chr, $pos, $id, $ref, $alt, $freq, $depth, $totBases, $altBases, $change, $bin);
my $tCol = 10;

my %freqBins;

my $binSize = 2;
my $maxBin = 100;

open (FILE, $data{ssm_file}) or die "Couldn't open $data{ssm_file}\n";

while ($l = <FILE>)
{
	chomp $l;
	if ($l =~ /^#/)
	{
	}
	else
	{
        @f = split(/\t/, $l);

        $chr = $f[0];
        $pos = $f[1];
        $id = $f[2];
        $ref = $f[3];
        for $alt (split(/,/, $f[4]))
        {
			$change = "";
            if (length($ref) == length($alt))
            {
                $change = "$ref$alt";
                if (exists $comp{$change})
                {
                    $change = $comp{$change};
                }
			}
            elsif (length($ref) > length($alt))
            {
				if (length($ref) - length($alt) > 3)
				{
                	$change = "del_4";
				}
				else
				{
					$change = "del_1";
				}
            }
            else
            {
				if (length($alt) - length($ref) > 3)
				{
                	$change = "ins_4";
				}
				else
				{
					$change = "ins_1";
				}
            }
			if ($f[8] eq "DP:FDP:SDP:SUBDP:AU:CU:GU:TU")     # strelka format
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
					die "Couldn't parse snv freq\n";
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

			$bin = int(($freq * 100) / $binSize) * $binSize;
			if ($bin == $maxBin)
			{
				$bin = $maxBin - $binSize;
			}

			$freqBins{$change}{$bin}++;


		}
	}
}

my $numBins = $maxBin / $binSize;
for $change (qw/del_1 del_4 ins_1 ins_4 CA CG CT TA TC TG/)
{
	for ($bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		unless (exists $freqBins{$change}{$bin})
		{
			$freqBins{$change}{$bin} = 0;
		}
	}
}

my @labels;
for (my $bin = 0; $bin <= $maxBin; $bin += $binSize * 5)
{
	push (@labels, "$bin%");
}
my @ticks;
for (my $i = 0; $i <= ($maxBin/$binSize); $i += 5)
{
	push (@ticks, $i + 0.2*$i + 0.1);
}

my @snvVals;
my @snvCols;
for $change (qw/CA CG CT TA TC TG/)
{
	push (@snvCols, $colours{$change});
	for (my $bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		push(@snvVals, $freqBins{$change}{$bin});
	}
}

my @indelVals;
my @indelCols;
for $change (qw/del_1 del_4 ins_1 ins_4/)
{
	push (@indelCols, $colours{$change});
	for (my $bin = 0; $bin < $maxBin; $bin += $binSize)
	{
		push(@indelVals, $freqBins{$change}{$bin});
	}
}


my $rfile;
open ($rfile, ">$fileName-$plotType.Rcode") or die "Couldn't open $fileName-$plotType.Rcode\n";

printArrayToR($rfile, "snvVals", \@snvVals);
printArrayToR($rfile, "indelVals", \@indelVals);
printArrayToR($rfile, "ticks", \@ticks);

printQuotedArrayToR($rfile, "snvCols", \@snvCols);
printQuotedArrayToR($rfile, "indelCols", \@indelCols);
printQuotedArrayToR($rfile, "labels", \@labels);


print $rfile "png(filename = \"$fileName-snv_vaf.png\", width = 1000, height = 250)\n";
print $rfile "par(mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(snvVals,6,$numBins,byrow=T), col=snvCols, main=\"$data{tumour} Somatic SNV VAF Distribution\", xlab=\"Variant Allele Frequency\", ylab=\"SNV Count\", border=NA)\n";
print $rfile "axis(side=1, at=ticks, lab=labels)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-snv_vaf-no_sample.png\", width = 1000, height = 200)\n";
print $rfile "par(mar=c(4,5,1,0) + 0.1)\n";
print $rfile "barplot(matrix(snvVals,6,$numBins,byrow=T), col=snvCols, main=\"\", xlab=\"SNV Variant Allele Frequency\", ylab=\"SNV Count\", border=NA)\n";
print $rfile "axis(side=1, at=ticks, lab=labels)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-indel_vaf.png\", width = 1000, height = 250)\n";
print $rfile "par(mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(indelVals,4,$numBins,byrow=T), col=indelCols, main=\"$data{tumour} Somatic Indel Variant Allele Frequency\", xlab=\"Frequency\", ylab=\"Count\", border=NA)\n";
print $rfile "axis(side=1, at=ticks, lab=labels)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-indel_vaf-no_sample.png\", width = 1000, height = 200)\n";
print $rfile "par(mar=c(4,5,1,0) + 0.1)\n";
print $rfile "barplot(matrix(indelVals,4,$numBins,byrow=T), col=indelCols, main=\"\", xlab=\"Indel Variant Allele Frequency\", ylab=\"Indel Count\", border=NA)\n";
print $rfile "axis(side=1, at=ticks, lab=labels)\n";
print $rfile "dev.off()\n";


close $rfile;
`Rscript $fileName-$plotType.Rcode`;



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







