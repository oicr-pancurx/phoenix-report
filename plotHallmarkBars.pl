#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $dataFile = shift;
my $fileName = shift;
my $plotType = "hallmark_bars";

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


my $dsbrColour = "#332288";
my $mmrColour = "#117733";
my $bgColour = "#DDDDDD";




my $rfile;
open ($rfile, ">$fileName-$plotType.Rcode") or die "Couldn't open $fileName-$plotType.Rcode\n";


print $rfile "png(filename = \"$fileName-DSBR_score_bar.png\", width = 1000, height = 30)\n";
print $rfile "par(mar=c(0.8,0,0,0) + 0.1)\n";
print $rfile "barplot($data{dsbr_score},xlim=c(0,10),horiz=T,xaxt=\"n\",yaxt=\"n\",bty=\"n\",col=\"$dsbrColour\",border=F)\n";
print $rfile "abline(v=c(0,10))\n";
print $rfile "abline(v=c(5),lty=\"dashed\",col=\"grey\")\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-MMR_score_bar.png\", width = 1000, height = 30)\n";
print $rfile "par(mar=c(0.8,0,0,0) + 0.1)\n";
print $rfile "barplot($data{mmr_score},xlim=c(0,4),horiz=T,xaxt=\"n\",yaxt=\"n\",bty=\"n\",col=\"$mmrColour\",border=F)\n";
print $rfile "abline(v=c(0,4))\n";
print $rfile "abline(v=c(3),lty=\"dashed\",col=\"grey\")\n";
print $rfile "dev.off()\n";


print $rfile "par(bg = \"#DDDDDD\")\n";
print $rfile "png(filename = \"$fileName-DSBR_score_bar_simple.png\", width = 100, height = 30)\n";
print $rfile "par(mar=c(0,0,0,0))\n";
print $rfile "barplot($data{dsbr_score},xlim=c(0,10),horiz=T,xaxt=\"n\",yaxt=\"n\",xaxs=\"i\",yaxs=\"i\",bty=\"n\",col=\"$dsbrColour\",border=\"$dsbrColour\")\n";
print $rfile "barplot(c(10),xlim=c(0,10),horiz=T,xaxt=\"n\",yaxt=\"n\",xaxs=\"i\",yaxs=\"i\",bty=\"n\",col=\"$bgColour\",border=\"$bgColour\",add=T)\n";
print $rfile "barplot($data{dsbr_score},xlim=c(0,10),horiz=T,xaxt=\"n\",yaxt=\"n\",xaxs=\"i\",yaxs=\"i\",bty=\"n\",col=\"$dsbrColour\",border=\"$dsbrColour\",add=T)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-MMR_score_bar_simple.png\", width = 100, height = 30)\n";
print $rfile "par(mar=c(0,0,0,0))\n";
print $rfile "barplot($data{mmr_score},xlim=c(0,4),horiz=T,xaxt=\"n\",yaxt=\"n\",xaxs=\"i\",yaxs=\"i\",bty=\"n\",col=\"$mmrColour\",border=\"$mmrColour\")\n";
print $rfile "barplot(c(4),xlim=c(0,4),horiz=T,xaxt=\"n\",yaxt=\"n\",xaxs=\"i\",yaxs=\"i\",bty=\"n\",col=\"$bgColour\",border=\"$bgColour\",add=T)\n";
print $rfile "barplot($data{mmr_score},xlim=c(0,4),horiz=T,xaxt=\"n\",yaxt=\"n\",xaxs=\"i\",yaxs=\"i\",bty=\"n\",col=\"$mmrColour\",border=\"$mmrColour\",add=T)\n";
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







