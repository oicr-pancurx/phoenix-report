#!/usr/bin/perl

use strict;
use warnings;


my $varFile = $ARGV[0];


my $l;
my @header;
my @row;
my %var;

my %germClass = (
	"germline snp" => 1,
	"germline indel" => 1
);

my %germType = (
	"frameshift" => 1,
	"nonframeshift" => 1,
	"stopgain" => 1,
	"stoploss" => 1,
	"splicing" => 1,
	"nonsynonymous" => 1,
);

my %somType = (
	"frameshift" => 1,
	"nonframeshift" => 1,
	"stopgain" => 1,
	"stoploss" => 1,
	"splicing" => 1,
	"nonsynonymous" => 1,
	"deletion breakpoint" => 1,
	"duplication breakpoint" => 1,
	"inversion breakpoint" => 1,
	"translocation breakpoint" => 1,
	"deletion potential fusion" => 1,
	"duplication potential fusion" => 1,
	"inversion potential fusion" => 1,
	"translocation potential fusion" => 1,
	"homozygous deletion" => 2,
);

my %genes;

open (FILE, $varFile) or die "Couldn't open $varFile\n";

$l = <FILE>;
chomp $l;
for my $f (split(/,/, $l))
{
	push (@header, $f);
}

while ($l = <FILE>)
{
	%var = ();
	chomp $l;
	@row = split(/,/, $l);
	for (my $i = 0; $i < scalar(@row); $i++)
	{
		$var{$header[$i]} = $row[$i];
	}

	if (exists $germClass{$var{mutation_class}})
	{
		if (exists $germType{$var{mutation_type}})
		{
			unless ($var{rarity} eq "common")
			{
				$genes{$var{gene}}{count} += $germType{$var{mutation_type}};
				push( @{ $genes{$var{gene}}{hits} }, $l );
			}
		}
	}
	else	# assuming somatic
	{
		if (exists $somType{$var{mutation_type}})
		{
			$genes{$var{gene}}{count} += $somType{$var{mutation_type}};
			push( @{ $genes{$var{gene}}{hits} }, $l );
		}
	}

	if ($var{ab_counts} =~ /0\./)
	{
		$genes{$var{gene}}{count} += 1;
	}
}


for my $g (sort keys %genes)
{
	if ($genes{$g}{count} > 1)
	{
		for $l (@{ $genes{$g}{hits} })
		{
			print $l . "\n";
		}
	}
}



