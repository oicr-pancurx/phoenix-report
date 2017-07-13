#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $dataFile = shift;
my $fileName = shift;
my $plotType = "indel_sv_bins";

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
    "del_1" => "#DC050C",
    "del_4" => "#F7EE55",
    "ins_1" => "#7BAFDE", # transversion
    "ins_4" => "#D6C1DE",
    "DEL" => "#F1932D", # transversion
    "DUP" => "#90C987",
    "INV" => "#1965B0",
    "TRA" => "#882E72",
);


my @svTypes = qw/DEL DUP INV/;
my @svBins = qw/100 1k 10k 100k 1m 10m >10m/;

my @indels1;
my @indels4;
my @svs;

push (@indels1, $data{del_1_count});
push (@indels1, $data{ins_1_count});
push (@indels1, "NA");
push (@indels1, "NA");

push (@indels4, "NA");
push (@indels4, "NA");
push (@indels4, $data{del_4_count});
push (@indels4, $data{ins_4_count});

for my $bin (@svBins)
{
	for my $type (@svTypes)
	{
		push(@svs, $data{"$type:$bin"});
	}
	push(@svs, "NA");
}

push(@svs, 0);
push(@svs, 0);
push(@svs, 0);
push(@svs, $data{sv_tra_count});

my @indels1Col = ($colours{del_1},$colours{ins_1});
my @indels4Col = ($colours{del_4},$colours{ins_4});
my @svsCol = ($colours{DEL}, $colours{DUP}, $colours{INV}, $colours{TRA});

my $rfile;
open ($rfile, ">$fileName-$plotType.Rcode") or die "Couldn't open $fileName-$plotType.Rcode\n";

printArrayToR($rfile, "indels1", \@indels1);
printArrayToR($rfile, "indels4", \@indels4);
printArrayToR($rfile, "svs", \@svs);

printQuotedArrayToR($rfile, "indels1_col", \@indels1Col);
printQuotedArrayToR($rfile, "indels4_col", \@indels4Col);
printQuotedArrayToR($rfile, "svs_col", \@svsCol);


print $rfile "png(filename = \"$fileName-$plotType-800x600.png\", width = 800, height = 600)\n";
print $rfile "par(fig=c(0,0.25,0,1),mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(indels1,nrow=2),xlab=\"Indel Sizes\",names.arg=c(\"1-3\",\">4\"),ylim=c(0,1000),col=indels1_col,border=NA, space=0.0625, cex.lab=1.5, ylab=\"Indel Count\")\n";
print $rfile "barplot(matrix(indels4,nrow=2),ylim=c(0,1000),col=indels4_col,add=T,border=NA, space=0.0625,axes=F)\n";
print $rfile "par(fig=c(0.23,1,0,1), new = TRUE,mar=c(4, 0, 4, 4) + 0.1)\n";
print $rfile "barplot(matrix(svs,nrow=4),axes=F,names.arg=c(\"100\", \"1k\", \"10k\", \"100k\", \"1m\",\"10m\",\">10m\", \"TRA\"),ylim=c(0,150),col=c(\"#F1932D\",\"#90C987\",\"#1965B0\",\"#882E72\"),border=NA, space=0.0625, xlab=\"SV Sizes\", cex.lab=1.5)\n";
print $rfile "axis(4, line=-1)\nmtext(\"SV Count\",4,line=2,cex=1.5)\n";
print $rfile "par(new=T,fig=c(0,1,0,1),mar=c(4,5,4,4)+0.1)\n";
#print $rfile "plot(c(0,0),type=\"n\",main=\"$data{tumour}\",axes=F,ylab=\"\",xlab=\"\",pty=NA,cex.main=2)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-800x600-label.png\", width = 800, height = 600)\n";
print $rfile "par(fig=c(0,0.25,0,1),mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(indels1,nrow=2),xlab=\"Indel Sizes\",names.arg=c(\"1-3\",\">4\"),ylim=c(0,1000),col=indels1_col,border=NA, space=0.0625, cex.lab=1.5, ylab=\"Indel Count\")\n";
print $rfile "barplot(matrix(indels4,nrow=2),ylim=c(0,1000),col=indels4_col,add=T,border=NA, space=0.0625,axes=F)\n";
print $rfile "par(fig=c(0.23,1,0,1), new = TRUE,mar=c(4, 0, 4, 4) + 0.1)\n";
print $rfile "barplot(matrix(svs,nrow=4),axes=F,names.arg=c(\"100\", \"1k\", \"10k\", \"100k\", \"1m\",\"10m\",\">10m\", \"TRA\"),ylim=c(0,150),col=c(\"#F1932D\",\"#90C987\",\"#1965B0\",\"#882E72\"),border=NA, space=0.0625, xlab=\"SV Sizes\", cex.lab=1.5)\n";
print $rfile "axis(4, line=-1)\nmtext(\"SV Count\",4,line=2,cex=1.5)\n";
print $rfile "par(new=T,fig=c(0,1,0,1),mar=c(4,5,4,4)+0.1)\n";
print $rfile "plot(c(0,0),type=\"n\",main=\"$data{tumour}\",axes=F,ylab=\"\",xlab=\"\",pty=NA,cex.main=2)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-800x600-stretch.png\", width = 800, height = 600)\n";
print $rfile "par(fig=c(0,0.20,0,1),mar=c(0,0,0,0) + 0.1)\n";
print $rfile "barplot(matrix(indels1,nrow=2),xlab=\"\",ylim=c(0,1000),col=indels1_col,border=NA, space=0.0625, cex.lab=1.5, ylab=\"\", axes=F)\n";
print $rfile "barplot(matrix(indels4,nrow=2),ylim=c(0,1000),col=indels4_col,add=T,border=NA, space=0.0625,axes=F)\n";
print $rfile "par(fig=c(0.18,1,0,1), new = TRUE,mar=c(0, 0, 0, 0) + 0.1)\n";
print $rfile "barplot(matrix(svs,nrow=4),axes=F,ylim=c(0,150),col=c(\"#F1932D\",\"#90C987\",\"#1965B0\",\"#882E72\"),border=NA, space=0.0625, xlab=\"\", cex.lab=1.5)\n";
#print $rfile "axis(4, line=-1)\nmtext(\"SV Count\",4,line=2,cex=1.5)\n";
#print $rfile "par(new=T,fig=c(0,1,0,1),mar=c(4,5,4,4)+0.1)\n";
#print $rfile "plot(c(0,0),type=\"n\",main=\"$data{tumour}\",axes=F,ylab=\"\",xlab=\"\",pty=NA,cex.main=2)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-400x300.png\", width = 400, height = 300)\n";
print $rfile "par(fig=c(0,0.32,0,1),mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(indels1,nrow=2),xlab=\"Indel Sizes\",names.arg=c(\"1-3\",\">4\"),ylim=c(0,1000),col=indels1_col,border=NA, space=0.0625, cex.lab=1.3, ylab=\"Indel Count\", cex.names=0.8, las=2)\n";
print $rfile "barplot(matrix(indels4,nrow=2),ylim=c(0,1000),col=indels4_col,add=T,border=NA, space=0.0625,axes=F)\n";
print $rfile "par(fig=c(0.30,1,0,1), new = TRUE,mar=c(4, 0, 4, 4.5) + 0.1)\n";
print $rfile "barplot(matrix(svs,nrow=4),axes=F,names.arg=c(\"100\", \"1k\", \"10k\", \"100k\", \"1m\",\"10m\",\">10m\", \"TRA\"),ylim=c(0,150),col=c(\"#F1932D\",\"#90C987\",\"#1965B0\",\"#882E72\"),border=NA, space=0.0625, xlab=\"SV Sizes\", cex.lab=1.3, cex.names=0.9, las=2)\n";
print $rfile "axis(4, line=-0.4, las=2)\nmtext(\"SV Count\",4,line=2.5,cex=1.3)\n";
print $rfile "par(new=T,fig=c(0,1,0,1),mar=c(4,5,4,4.5)+0.1)\n";
#print $rfile "plot(c(0,0),type=\"n\",main=\"$data{tumour}\",axes=F,ylab=\"\",xlab=\"\",pty=NA,cex.main=1.5)\n";
print $rfile "dev.off()\n";

print $rfile "png(filename = \"$fileName-$plotType-400x300-label.png\", width = 400, height = 300)\n";
print $rfile "par(fig=c(0,0.32,0,1),mar=c(4,5,4,0) + 0.1)\n";
print $rfile "barplot(matrix(indels1,nrow=2),xlab=\"Indel Sizes\",names.arg=c(\"1-3\",\">4\"),ylim=c(0,1000),col=indels1_col,border=NA, space=0.0625, cex.lab=1.3, ylab=\"Indel Count\", cex.names=0.8, las=2)\n";
print $rfile "barplot(matrix(indels4,nrow=2),ylim=c(0,1000),col=indels4_col,add=T,border=NA, space=0.0625,axes=F)\n";
print $rfile "par(fig=c(0.30,1,0,1), new = TRUE,mar=c(4, 0, 4, 4.5) + 0.1)\n";
print $rfile "barplot(matrix(svs,nrow=4),axes=F,names.arg=c(\"100\", \"1k\", \"10k\", \"100k\", \"1m\",\"10m\",\">10m\", \"TRA\"),ylim=c(0,150),col=c(\"#F1932D\",\"#90C987\",\"#1965B0\",\"#882E72\"),border=NA, space=0.0625, xlab=\"SV Sizes\", cex.lab=1.3, cex.names=0.9, las=2)\n";
print $rfile "axis(4, line=-0.4, las=2)\nmtext(\"SV Count\",4,line=2.5,cex=1.3)\n";
print $rfile "par(new=T,fig=c(0,1,0,1),mar=c(4,5,4,4.5)+0.1)\n";
print $rfile "plot(c(0,0),type=\"n\",main=\"$data{tumour}\",axes=F,ylab=\"\",xlab=\"\",pty=NA,cex.main=1.5)\n";
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







