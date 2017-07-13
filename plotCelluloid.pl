#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $dataFile = shift;
my $fileName = shift;
my $plotType = "celluloid_contour";

my $l;
my @header;
my %data;
my %values;
my $sample;

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


# /.mounts/labs/PCSI/analysis/oicr/PCSI0292/PCSI_0292_Pa_P_526/wgs/bwa/0.6.2/celluloid/v11.2/solution/segments_PCSI_0292_Pa_P_526.txt.sorted.bed.gz
my $image = $data{seg_file};

$image =~ s/segments.*//;
$image .= "contour_$data{tumour}.png";

`cp $image $fileName-$plotType.png`;

