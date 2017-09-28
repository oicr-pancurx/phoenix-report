#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $dataFile = shift;
my $fileName = shift;
my $plotType = "histogram";

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

my $snvCount = $data{snv_count};
my $indelCount = $data{indel_count};
my $svCount = $data{sv_count};
my $neoCount = $data{neo_antigens};

my %histVals;
my %histCols;


my $historyFile = "/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/pdac_20170926_matrix.csv";
## read file for population numbers
my %history;
my @historySamps;

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

my %histRange400 = (
    "snv_count" => "c(0,30000)",
    "indel_count" => "c(0,3000)",
    "sv_count" => "c(0,300)",
    "neo_antigens" => "c(0,200)",
);

my %denseMin = (
	"snv_count" => 500,
	"indel_count" => 50,
	"sv_count"=> 10,
	"neo_antigens" => 10,
);

my %denseMax = (
	"snv_count" => 50000,
	"indel_count" => 5000,
	"sv_count"=> 500,
	"neo_antigens" => 200,
);

my %valuePos = (
	"snv_count" => $snvCount,
	"indel_count" => $indelCount,
	"sv_count"=> $svCount,
	"neo_antigens" => $neoCount,
);


my %boxMax = (
    "snv_count" => 20000,
    "indel_count" => 2000,
    "sv_count" => 400,
    "neo_antigens" => 200,
);

my %histBreaks = (
	"snv_count" => 500,
	"indel_count" => 5000,
	"sv_count"=> 50,
	"neo_antigens" => 500,
);



for my $type (qw/snv_count indel_count sv_count neo_antigens/)
{
    @{ $histVals{$type} } = valueVector(\%history, \@historySamps, $type);

    for (my $i = 0; $i < scalar(@{ $histVals{$type} }); $i++)
    {
        push(@{ $histCols{$type} }, $histColours{$type});
    }
    push (@{ $histCols{$type} }, $histColours{compass});


	if ($valuePos{$type} > $boxMax{$type})
	{
		$boxMax{$type} = int($valuePos{$type}*1.05);
	}
}

my $snvPR = getPercentRank($snvCount, \@{ $histVals{snv_count} });
my $indelPR = getPercentRank($indelCount, \@{ $histVals{indel_count} });
my $svPR = getPercentRank($svCount, \@{ $histVals{sv_count} });
my $neoPR = getPercentRank($neoCount, \@{ $histVals{neo_antigens} });

my $boxRange;

push (@{ $histVals{snv_count} }, $snvCount);
push (@{ $histVals{indel_count} }, $indelCount);
push (@{ $histVals{sv_count} }, $svCount);
push (@{ $histVals{neo_antigens} }, $neoCount);


my $rfile;
open ($rfile, ">$fileName-$plotType.Rcode") or die "Couldn't open $fileName-$plotType.Rcode\n";

my %titleHash = (
	"snv_count" => "$snvCount SNVs ($snvPR% PR)",
	"indel_count" => "$indelCount Indels ($indelPR% PR)",
	"sv_count" => "$svCount SVs ($svPR% PR)",
	"neo_antigens" => "$neoCount Neo-antigens ($neoPR% PR)",
);

my %titleDense = (
	"snv_count" => "$snvCount SNVs ($snvPR% PR)",
	"indel_count" => "$indelCount Indels ($indelPR% PR)",
	"sv_count" => "$svCount SVs ($svPR% PR)",
	"neo_antigens" => "$neoCount Neo-antigens ($neoPR% PR)",
);
my $numSamps = scalar(@historySamps);

my $vertValue;

for my $type (qw/snv_count indel_count sv_count neo_antigens/)
{
    printArrayToR($rfile,$type,\@{ $histVals{$type} });
    printQuotedArrayToR($rfile,"${type}_col",\@{ $histCols{$type} });

    print $rfile "png(\"$fileName-$plotType-${type}-240.png\",width=240,height=100)\n";
    print $rfile "par(mar=c(1,3,1,0)+0.1, las=1,cex=0.5)\n";
    print $rfile "bp = barplot(${type}\[sort($type,index.return=T)\$ix], ylim=$histRange{$type}, col=as.vector(${type}_col[sort($type,index.return=T)\$ix]), border=NA, space=0)\n";
    print $rfile "abline(v=bp[match(length($type),sort($type,index.return=T)\$ix)])\n";
    print $rfile "dev.off()\n";

    print $rfile "\n";



    print $rfile "png(\"$fileName-$plotType-${type}-400.png\",width=400,height=300)\n";
    print $rfile "par(mar=c(4,5,4,0)+0.1, las=1,cex=1)\n";
    print $rfile "bp = barplot(${type}\[sort($type,index.return=T)\$ix], ylim=$histRange400{$type}, col=as.vector(${type}_col[sort($type,index.return=T)\$ix]), border=NA, space=0, main=\"$titleHash{$type}\", xlab=\"\", ylab=\"\", cex.main=2)\n";
	print $rfile "mtext(text=\"Series of $numSamps PDAC Tumours\", side=1, line=1)\n";
	print $rfile "mtext(text=\"Mutational Load\", side=2, line=4,las=0)\n";
    print $rfile "abline(v=bp[match(length($type),sort($type,index.return=T)\$ix)])\n";
    print $rfile "dev.off()\n";

    print $rfile "\n";

	$vertValue = $data{$type};
	if ($vertValue > $denseMax{$type})
	{
		$vertValue = $denseMax{$type};
	}
	if ($vertValue < $denseMin{$type})
	{
		$vertValue = $denseMin{$type};
	}

	# plot actual histograms

    print $rfile "png(\"$fileName-density-${type}-240.png\",width=240,height=100)\n";
    print $rfile "par(mar=c(2,1,1,0)+0.1, las=1,cex=0.6)\n";
	print $rfile "d <- density(${type})\n";
    print $rfile "plot(d, xlim=c($denseMin{$type},$denseMax{$type}), log=\"x\", main=\"\", type=\"n\", yaxt=\"n\",xaxt=\"n\",bty=\"n\")\n";
	print $rfile "polygon(x=c(d\$x,0.1), y=c(d\$y,0), col=\"$histColours{$type}\",border=F)\n";
    print $rfile "abline(v=$vertValue)\n";
	print $rfile "axis(side=1,at=c(1,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000),labels=c(1,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000))\n";
    print $rfile "dev.off()\n";

    print $rfile "\n";


    print $rfile "png(\"$fileName-density-${type}-400.png\",width=400,height=300)\n";
    print $rfile "par(mar=c(4,1,4,0)+0.1, las=1,cex=1)\n";
	print $rfile "d <- density(${type})\n";
    print $rfile "plot(d, xlim=c($denseMin{$type},$denseMax{$type}), log=\"\",type=\"n\", yaxt=\"n\",xaxt=\"n\",bty=\"n\", main=\"$titleDense{$type}\", cex.main=2, xlab=\"\")\n";
	print $rfile "polygon(x=c(d\$x,0.1), y=c(d\$y,0), col=\"$histColours{$type}\",border=F)\n";
	print $rfile "mtext(text=\"Mutational load in $numSamps PDAC Tumours\", side=1, line=3,cex=1.5)\n";
    print $rfile "abline(v=$vertValue)\n";
	print $rfile "axis(side=1,at=c(1,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000),labels=c(1,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000))\n";
    print $rfile "dev.off()\n";

    print $rfile "\n";


	# plot hist + boxplots
	$boxRange = "c(0,$boxMax{$type})";
	print $rfile "png(\"$fileName-histbox-${type}-240.png\",width=240,height=100)\n";
	print $rfile "par(mar=c(0,1,0,1)+0.1,fig=c(0,1,0.5,1))\n";
	print $rfile "hist($type, xlim=$boxRange,breaks=$histBreaks{$type},axes=F,xaxs=\"i\", yaxs=\"i\",ylab=NA,xlab=NA,main=NA, col=\"$histColours{$type}\",border=NA)\n";
	print $rfile "par(mar=c(2,1,0,1)+0.1,fig=c(0,1,0,0.54),new=T)\n";
	print $rfile "boxplot($type, horizontal=T, ylim=$boxRange, col=\"$histColours{$type}\", xaxs=\"i\", yaxs=\"i\", frame=F,width=1, cex.axis=0.8,pch=20,cex=0.8)\n";
	print $rfile "points($valuePos{$type},1,cex=2,pch=4)\n";
	print $rfile "par(mar=c(2,1,0,1)+0.1,fig=c(0,1,0,1),new=T)\n";
	print $rfile "plot(c(0,1),xlim=$boxRange,type=\"n\",axes=F,xaxs=\"i\", yaxs=\"i\")\n";
	print $rfile "abline(v=$valuePos{$type},lwd=2)\n";
	print $rfile "dev.off()\n";
	print $rfile "\n";



}

close $rfile;
`Rscript $fileName-$plotType.Rcode`;



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

    return sprintf("%.0f", ($rank * 100));
}






