#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $dataFile = shift;
my $fileName = shift;
my $plotType = "snv_context";

my $l;
my @header;
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
    "ca" => "#E8601C",
    "cg" => "#F6C141",
    "ct" => "#CAE0AB", # transversion
    "ta" => "#4EB265",
    "tc" => "#5289C7", # transversion
    "tg" => "#B178A6",
);

my %fadeColours = (
    "ca" => "#E8601C66",
    "cg" => "#F6C14166",
    "ct" => "#CAE0AB66", # transversion
    "ta" => "#4EB26566",
    "tc" => "#5289C766", # transversion
    "tg" => "#B178A666",
);

my @changes = qw/ca cg ct ta tc tg/;
my @contexts = qw/aa ac ag at ca cc cg ct ga gc gg gt ta tc tg tt/;

my @changeHeights;
my @contextHeights;

my @changeCol;
my @contextCol;

my @changeLabel = qw/C>A C>G C>T T>A T>C T>G/;
my @contextLabel = qw/ANA ANC ANG ANT CNA CNC CNG CNT GNA GNC GNG GNT TNA TNC TNG TNT ANA ANC ANG ANT CNA CNC CNG CNT GNA GNC GNG GNT TNA TNC TNG TNT ANA ANC ANG ANT CNA CNC CNG CNT GNA GNC GNG GNT TNA TNC TNG TNT ANA ANC ANG ANT CNA CNC CNG CNT GNA GNC GNG GNT TNA TNC TNG TNT ANA ANC ANG ANT CNA CNC CNG CNT GNA GNC GNG GNT TNA TNC TNG TNT ANA ANC ANG ANT CNA CNC CNG CNT GNA GNC GNG GNT TNA TNC TNG TNT/;

my $changeMax = 0;
my $contextMax = 0;


for my $change (@changes)
{
	if ($data{"snv_$change"} > $changeMax)
	{
		$changeMax = $data{"snv_$change"};
	}
	for my $context (@contexts)
	{
		push(@contextHeights, $data{"snv_${change}_$context"} / $data{snv_count});
		if ($data{"snv_${change}_$context"} > $contextMax)
		{
			$contextMax = $data{"snv_${change}_$context"};
		}
		push(@contextCol, $colours{$change});
	}
}

for my $change (@changes)
{
	for my $context (@contexts)
	{
		push(@changeHeights, (($data{"snv_$change"} / $changeMax) * $contextMax) / $data{snv_count});
		push(@changeCol, $fadeColours{$change});
	}
}






my $rfile;
open ($rfile, ">$fileName-$plotType.Rcode") or die "Couldn't open $fileName-$plotType.Rcode\n";

printArrayToR($rfile, "change", \@changeHeights);
printArrayToR($rfile, "context", \@contextHeights);

printQuotedArrayToR($rfile, "change_cols", \@changeCol);
printQuotedArrayToR($rfile, "context_cols", \@contextCol);

printQuotedArrayToR($rfile, "change_labs", \@changeLabel);
printQuotedArrayToR($rfile, "context_labs", \@contextLabel);

print $rfile "png(filename = \"$fileName-$plotType-800x600.png\", width = 800, height = 600)\n";
print $rfile "par(mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(change,16,6), col=change_cols, main=\"\", xlab=\"\", ylab=\"\", border=NA, beside=TRUE, axes=FALSE)\n";
print $rfile "barplot(matrix(context,16,6), col=context_cols, main=\"\", xlab=\"\", ylab=\"SNV Trinucleotide Context Proportion\", border=NA, beside=TRUE, add=TRUE, names.arg=context_labs, cex.names=0.6, las=3, cex.lab=1.5)\n";
print $rfile "axis(1, at=c(9,26,43,60,77,94), labels=change_labs, tick=FALSE, line=2, cex.axis=1.5)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-800x600-label.png\", width = 800, height = 600)\n";
print $rfile "par(mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(change,16,6), col=change_cols, main=\"\", xlab=\"\", ylab=\"\", border=NA, beside=TRUE, axes=FALSE)\n";
print $rfile "barplot(matrix(context,16,6), col=context_cols, main=\"$data{tumour}\", xlab=\"\", ylab=\"SNV Trinucleotide Context Proportion\", border=NA, beside=TRUE, add=TRUE, names.arg=context_labs, cex.names=0.6, las=3, cex.lab=1.5, cex.main=2)\n";
print $rfile "axis(1, at=c(9,26,43,60,77,94), labels=change_labs, tick=FALSE, line=2, cex.axis=1.5)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-800x600-stretch.png\", width = 800, height = 600)\n";
print $rfile "par(mar=c(0,0,0,0) + 0.1)\n";
print $rfile "barplot(matrix(change,16,6), col=change_cols, main=\"\", xlab=\"\", ylab=\"\", border=NA, beside=TRUE, axes=FALSE)\n";
print $rfile "barplot(matrix(context,16,6), col=context_cols, main=\"\", xlab=\"\", ylab=\"\", border=NA, beside=TRUE, add=TRUE, cex.names=0.6, las=3, cex.lab=1.5, cex.main=2, axes=F)\n";
#print $rfile "axis(1, at=c(9,26,43,60,77,94), labels=change_labs, tick=FALSE, line=2, cex.axis=1.5)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-400x300.png\", width = 400, height = 300)\n";
print $rfile "par(mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(change,16,6), col=change_cols, main=\"\", xlab=\"\", ylab=\"\", border=NA, beside=TRUE, axes=FALSE)\n";
print $rfile "barplot(matrix(context,16,6), col=context_cols, main=\"\", xlab=\"\", ylab=\"SNV Context Proportion\", border=NA, beside=TRUE, add=TRUE, cex.names=0.6, las=3, cex.lab=1.3)\n";
print $rfile "axis(1, at=c(9,26,43,60,77,94), labels=change_labs, tick=FALSE, line=0, cex.axis=1.3)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-400x300-label.png\", width = 400, height = 300)\n";
print $rfile "par(mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(change,16,6), col=change_cols, main=\"\", xlab=\"\", ylab=\"\", border=NA, beside=TRUE, axes=FALSE)\n";
print $rfile "barplot(matrix(context,16,6), col=context_cols, main=\"$data{tumour}\", xlab=\"\", ylab=\"SNV Context Proportion\", border=NA, beside=TRUE, add=TRUE, cex.names=0.6, las=3, cex.lab=1.3, cex.main=1.5)\n";
print $rfile "axis(1, at=c(9,26,43,60,77,94), labels=change_labs, tick=FALSE, line=0, cex.axis=1.3)\n";
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







