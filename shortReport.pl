#!/usr/bin/perl

use strict;
use warnings;
use Time::Local;

my $summaryFile = shift;
my $variantsFile = shift;
my $outPath = shift;

unless (defined $outPath)
{
	$outPath = ".";
}

my (%data);
$data{out_path} = $outPath;

my $commonLocation = "/.mounts/labs/PCSI/reports/phoenix_common";

my $l;
my (@f, @header);
$data{version} = "v0.2";

# read summary file
open (FILE, $summaryFile) or die "Couldn't open $summaryFile\n";

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


# read variants file
my (%genes, %row);
my $i;
open (FILE, $variantsFile) or die "Couldn't open $variantsFile\n";

while ($l = <FILE>)
{
    chomp $l;
    if ($l =~ /^donor/)
    {
        @header = split (/,/, $l);
    }
    else
    {
        %row = ();
        $i = 0;
        for my $val (split(/,/, $l))
        {
            $row{$header[$i]} = $val;
            $i++;
        }
		for my $h (qw/ploidy copy_number ab_counts gene_position maf_mean maf_p/)
		{
			$genes{$row{gene}}{$h} = $row{$h};
		}

		for my $h (qw/mutation_class mutation_type position fusion_genes base_change tumour_freq tumour_depth normal_freq normal_depth nuc_context aa_context dbsnp cosmic cadd_phred rarity 1000G_all ExAC_all ESP6500siv2_all clinvar cosmic_census_flag cosmic_census_data/)
		{
			$genes{$row{gene}}{variants}{"$row{mutation_type},$row{position},$row{base_change}"}{$h} = $row{$h};
		}
    }
}
close FILE;


$data{plot_dir} = "./plots/";
$data{common_dir} = "./common/";

`mkdir -p $outPath/plots`;
`ln -s $commonLocation $outPath/common`;

my $currentTime = localtime;


my %nonsilent =  (
	"frameshift" => 1,
	"nonframeshift" => 1,
	"nonsynonymous" => 1,
	"stoploss" => 1,
	"stopgain" => 1,
	"splicing" => 1,
	"strong amplification" => 1,
	"homozygous deletion" => 1,
	"deletion breakpoint" => 1,
	"duplication breakpoint" => 1,
	"inversion breakpoint" => 1,
	"translocation breakpoint" => 1,
	# leaving altered promoters and potential fusions for now
);


my %limsTissue = (
	"Ad" => "Adipose",
	"Ap" => "Appendix",
	"Ag" => "Adrenal Gland",
	"As" => "Ascites",
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


my $headerName = $data{tumour};

for my $id (split(/ /, $data{external_id}))
{
	if ($id =~ /^COMP/)
	{
		$headerName = $id;
	}
}



### Print website
my $html;
open ($html, ">$outPath/$data{tumour}_slide.html") or die "Couldn't open >$outPath/$data{tumour}_slide.html\n";


printHtmlHead($html, \%data);

print $html "<div>\n";
printHeader($html, \%data, $currentTime, $headerName);

printSlideOne($html, \%data, \%genes, \%limsTissue);

#printFooter($html, \%data); ##?

print $html "</div>\n";


print $html "<br><div>\n";

printHeader($html, \%data, $currentTime, $headerName);

printSlideTwo($html, \%data, \%genes);

print $html "</div>\n";


print $html "<br><div>\n";

printHeader($html, \%data, $currentTime, $headerName);

printSlideThree($html, \%data);

print $html "</div>\n";


print $html "<br><div>\n";

printHeader($html, \%data, $currentTime, $headerName);

printSlideFour($html, \%data);

print $html "</div>\n";


print $html "</body>\n";
print $html "</html>\n";




sub printHtmlHead
{
	my $html = shift;
	my $data = shift;

	print $html "<!DOCTYPE html>
<html>
<head>
<title>$data->{tumour} Genome Summary</title>

<link rel=\"stylesheet\" type=\"text/css\" href=\"$data->{common_dir}/slide_report_style.css\">

<script src=\"$data->{common_dir}/sorttable.js\"></script>
<style type=\"text/css\">
th, td {
  padding: 3px !important;
}
table
{
border-collapse:collapse;
}
/* Sortable tables */
table.sortable thead {
        background-color:#FFFFFF;
        color:#000000;
        font-weight: bold;
        cursor: default;
}
</style>

</head>\n";

print $html "<body>\n";

}


sub printHeader
{
	my $html = shift;
	my $data = shift;
	my $currentTime = shift;
	my $sample = shift;


	print $html "<table border=\"0\" style=\"width:100%\"><tr style=\"vertical-align:top;\">
	<td style=\"width:33.33%\"><h4>For Research Purposes Only</h4></td>
	<td style=\"width:33.33%\"><h4 style=\"text-align:center;\"></h4></td>
	<td style=\"width:33.33%\"><h4 style=\"text-align:right;\">$currentTime</h4></td></tr></table>\n";

	print $html "<h2>$sample Genomics Summary</h2>\n";
	


}


sub printSlideOne
{
	my $html = shift;
	my $data = shift;
	my $genes = shift;
	my $limsTissue = shift;

	my %qcIcon;

	my $passIcon = "./common/pass.png";
	my $failIcon = "./common/fail.png";

	my %cutoffs = (
		"tumour_coverage" => 45,
		"normal_coverage" => 30,
		"cellularity" => 20,
		"genotype_concordance" => 85,
		"tumour_error_rate" => 1,
		"normal_error_rate" => 1,
		"normal_soft_clip" => 1,
		"tumour_soft_clip" => 1,

	);


	my $tumTissue = "NA";
	if ($data{tumour} =~ /^...._...._(..)_/)
	{
		if (exists $limsTissue->{$1})
		{
			$tumTissue = $limsTissue->{$1};
		}
		else
		{
			$tumTissue = $1;
		}
	}

	my $norTissue = "NA";
	if ($data{normal} =~ /^...._...._(..)_/)
	{
		if (exists $limsTissue->{$1})
		{
			$norTissue = $limsTissue->{$1};
		}
		else
		{
			$norTissue = $1;
		}
	}


	# decimals to percent
	for my $test (qw/tumour_error_rate normal_error_rate tumour_soft_clip normal_soft_clip cellularity genotype_concordance/)
	{
		$data{$test} = $data{$test}*100;
	}


	# over values
	for my $test (qw/tumour_coverage normal_coverage cellularity genotype_concordance/)
	{
		if ($data{$test} >= $cutoffs{$test})
		{
			$qcIcon{$test} = $passIcon;
		}
		else
		{
			$qcIcon{$test} = $failIcon;
		}
	}

	# under values
	for my $test (qw/tumour_error_rate normal_error_rate tumour_soft_clip normal_soft_clip/)
	{
		if ($data{$test} < $cutoffs{$test})
		{
			$qcIcon{$test} = $passIcon;
		}
		else
		{
			$qcIcon{$test} = $failIcon;
		}
	}

	my %testName = (
		"tumour_coverage" => "Tumour Coverage (>$cutoffs{tumour_coverage}x)",
		"normal_coverage" => "Normal Coverage (>$cutoffs{normal_coverage}x)",
		"cellularity" => "Cellularity (>$cutoffs{cellularity}%)",
		"genotype_concordance" => "Genotype Concordance (>$cutoffs{genotype_concordance}%)",
		"tumour_error_rate" => "Tumour Error Rate (<$cutoffs{tumour_error_rate}%)",
		"normal_error_rate" => "Normal Error Rate (<$cutoffs{normal_error_rate}%)",
		"tumour_soft_clip" => "Tumour Soft Clip Rate (<$cutoffs{tumour_soft_clip}%)",
		"normal_soft_clip" => "Normal Soft Clip Rate (<$cutoffs{normal_soft_clip}%)",
	);

	my %testValue = (
		"tumour_coverage" => sprintf("%.1f",$data{tumour_coverage}) . "x",
		"normal_coverage" => sprintf("%.1f",$data{normal_coverage}) . "x",
		"cellularity" => "$data{cellularity}%",
		"genotype_concordance" => sprintf("%.01f", $data{genotype_concordance}) . "%",
		"tumour_error_rate" => sprintf("%.2f",$data{tumour_error_rate}) . "%",
		"normal_error_rate" => sprintf("%.2f",$data{normal_error_rate}) . "%",
		"tumour_soft_clip" => sprintf("%.2f",$data{tumour_soft_clip}) . "%",
		"normal_soft_clip" => sprintf("%.2f",$data{normal_soft_clip}) . "%",
	);

	print $html "<table style=\"width:49%;float:left\">\n";
	print $html "<tr><td colspan=\"2\"><b><h3>Sample Metadata</h3></td></tr><tr>\n";
	print $html "<td><b>Donor ID:</b> $data->{donor}</td>";
	print $html "<td><b>External ID:</b> $data->{external_id}</td>";
	print $html "</tr>\n<tr>\n";
	print $html "<td colspan=\"2\"><b>Tumour Sample:</b> $data->{tumour} ($tumTissue)</td>";
	print $html "</tr>\n<tr>\n";
	print $html "<td colspan=\"2\"><b>Normal Sample:</b> $data->{normal} ($norTissue)</td>";
	print $html "</tr>\n";
	print $html "</table>\n";


my %nonsilent =  (
	"frameshift" => 1,
	"nonframeshift" => 1,
	"nonsynonymous" => 1,
	"stoploss" => 1,
	"stopgain" => 1,
	"splicing" => 1,
	"strong amplification" => 1,
	"homozygous deletion" => 1,
	"deletion breakpoint" => 1,
	"duplication breakpoint" => 1,
	"inversion breakpoint" => 1,
	"translocation breakpoint" => 1,
	# leaving altered promoters and potential fusions for now
);


	# germline gene table
	my $didAprint = 0;
	my @germlineGeneOrder = qw/POLE POLD1 EPCAM MLH1 MSH2 MSH6 PMS2 BRCA1 BRCA2 PALB2 ATM APC MUTYH CDKN2A STK11 PRSS1 PRSS2 TP53 SMAD4 C12orf32/;

	print $html "<div style=\"width:49%;float:right;height:65%;overflow:auto;padding:0px;border:0px\">\n";

	print $html "<table style=\"width:100%\">\n";		# style=\"width:49%;float:right\">\n";
	print $html "<td colspan=\"3\"><b><h3>Germline Variants</h3></td></tr>\n";
	print $html "<tr><th>Gene</th><th>Variant</th><th style=\"text-align:center\">Copy Number</th></tr>\n";

	for my $g (@germlineGeneOrder)
	{
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($genes->{$g}{variants}{$v}{mutation_class} =~ /germline/)
			{
				if (exists $nonsilent{ $genes->{$g}{variants}{$v}{mutation_type} })
				{
					if ($genes->{$g}{variants}{$v}{rarity} ne "common")
					{
						printSimpleGeneRow($g, $v, $genes, $html, "");
						$didAprint = 1;
					}
				}
			}
		}
	}
	if ($didAprint == 0)
	{
		print $html "<tr><td colspan=\"3\">No variants of interest in $data{germline_snv_count} SNPs and $data{germline_indel_count} indels.</td></tr>\n";
	}


	# somatic gene table
	my @alwaysGenes = qw/KRAS TP53 CDKN2A SMAD4/;
	my @driverGenes = qw/ARID1A RNF43 TGFBR2 KDM6A/;

	print $html "<td colspan=\"3\"><b><h3>Somatic Mutations</h3></td></tr>\n";
	print $html "<tr><th>Gene</th><th>Mutation</th><th style=\"text-align:center\">Copy Number</th></tr>\n";


	my %printed;
	for my $g (@alwaysGenes)
	{
		$printed{$g} = 0;
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			unless ($genes->{$g}{variants}{$v}{mutation_class} =~ /germline/)
			{
				unless ($genes->{$g}{variants}{$v}{mutation_class} =~ /NA/)
				{
					printSimpleGeneRow($g, $v, $genes, $html, "");
					$printed{$g} = 1;
				}
			}
		}
	
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			if ($printed{$g} == 0)
			{
				printSimpleGeneRow($g, $v, $genes, $html, "gene only");
				$printed{$g} = 1;
			}
		}
	
	}

	for my $g (@driverGenes)
	{
		for my $v (sort keys %{ $genes->{$g}{variants} })
		{
			unless ($genes->{$g}{variants}{$v}{mutation_class} =~ /germline/)
			{
				unless ($genes->{$g}{variants}{$v}{mutation_class} =~ /NA/)
				{
					printSimpleGeneRow($g, $v, $genes, $html, "");
					$printed{$g} = 1;
				}
			}
		}
	}


	for my $g (sort keys %{ $genes })
	{
		unless (exists $printed{$g})
		{
			for my $v (sort keys %{ $genes->{$g}{variants} })
			{
				unless ($genes->{$g}{variants}{$v}{mutation_class} =~ /germline/)
				{
					if (exists $genes->{$g}{variants}{$v}{cosmic_census_flag})
					{
						if ($genes->{$g}{variants}{$v}{cosmic_census_flag} eq "cosmic_mutation")
						{
							printSimpleGeneRow($g, $v, $genes, $html, "");
						}
					}
				}
			}
		}
	}


	print $html "</table>\n";
	print $html "</div>\n";



	# QC Metrics

	print $html "<table style=\"width:49%;float:left\">\n";
	print $html "<tr><td colspan=\"3\"><b><h3>QC Metrics</h3></td></tr>\n";

	for my $test (qw/tumour_coverage normal_coverage tumour_error_rate normal_error_rate tumour_soft_clip normal_soft_clip cellularity genotype_concordance/)
	{
		print $html "<tr>";
		print $html "<td><b>$testName{$test}:</b></td>";
		print $html "<td style=\"text-align:left\">$testValue{$test}</td>";
		print $html "<td><img src=\"$qcIcon{$test}\" style=\"height:20px\"></td>";
		print $html "</tr>";
	}

	my $pipelineQC = "PASS";
	my $pipelineIcon = $passIcon;

	print $html "<tr>\n<td><b>Pipeline QCs:</b></td>";
	print $html "<td>$pipelineQC</td>";
	print $html "<td><img src=\"$pipelineIcon\" style=\"height:20px\"></td>";
	print $html "</table>\n";


	# oncoSlice

	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-onco_slice-1000_slide.png\" style=\"width:100%\">\n";


}

sub printSlideTwo
{
	my $html = shift;
	my $data = shift;
	my $genes = shift;

	# terribly hacky!
	my $snvText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep SNVs | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $indelText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep Indel | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $svText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep SVs | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $neoText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep Neo | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;

	$snvText =~ s/s \(/s<br>(/;
	$indelText =~ s/s \(/s<br>(/;
	$svText =~ s/s \(/s<br>(/;
	$neoText =~ s/s \(/s<br>(/;

	my $ploidy = "diploid ($data->{ploidy})";
	if ($data->{ploidy} > 2.3)
	{
		$ploidy = "polyploid ($data->{ploidy})";
	}

	my $ploidyForImg = $ploidy;
	$ploidyForImg =~ s/ \(.*//;

	my %classImg = (
		"dsbr_score" => "<img src=\"$data->{plot_dir}/$data->{tumour}-DSBR_score_bar_simple.png\" style=\"width:50px;height:18px\">",
		"mmr_score" => "<img src=\"$data->{plot_dir}/$data->{tumour}-MMR_score_bar_simple.png\" style=\"width:50px;height:18px\">",
		"cosmic_score" => "<img src=\"$data->{plot_dir}/$data->{tumour}-COSMIC_bar_simple.png\" style=\"width:50px;height:18px\">",
		"stable" => "<img src=\"./common/stable.png\" style=\"width:50px;height:18px\">",
		"locally rearranged" => "<img src=\"./common/locally_rearranged.png\" style=\"width:50px;height:18px\">",
		"scattered" => "<img src=\"./common/scattered.png\" style=\"width:50px;height:18px\">",
		"unstable" => "<img src=\"./common/unstable.png\" style=\"width:50px;height:18px\">",
		"classic" => "<img src=\"./common/classic.png\" style=\"width:50px;height:18px\">",
		"classical" => "<img src=\"./common/classical.png\" style=\"width:50px;height:18px\">",
		"basal-like" => "<img src=\"./common/basal-like.png\" style=\"width:50px;height:18px\">",
		"basal" => "<img src=\"./common/basal-like.png\" style=\"width:50px;height:18px\">",
		"quasimesenchymal" => "<img src=\"./common/quasi-mesenchymal.png\" style=\"width:50px;height:18px\">",
		"quasi-mesenchymal" => "<img src=\"./common/quasi-mesenchymal.png\" style=\"width:50px;height:18px\">",
		"exocrine" => "<img src=\"./common/exocrine-like.png\" style=\"width:50px;height:18px\">",
		"exocrine-like" => "<img src=\"./common/exocrine-like.png\" style=\"width:50px;height:18px\">",
		"progenitor" => "<img src=\"./common/progenitor.png\" style=\"width:50px;height:18px\">",
		"Pancreatic progenitor" => "<img src=\"./common/progenitor.png\" style=\"width:50px;height:18px\">",
		"ADEX" => "<img src=\"./common/adex.png\" style=\"width:50px;height:18px\">",
		"adex" => "<img src=\"./common/adex.png\" style=\"width:50px;height:18px\">",
		"immunogenic" => "<img src=\"./common/immunogenic.png\" style=\"width:50px;height:18px\">",
		"squamous" => "<img src=\"./common/squamous.png\" style=\"width:50px;height:18px\">",
		"diploid" => "<img src=\"./common/diploid.png\" style=\"width:50px;height:18px\">",
		"polyploid" => "<img src=\"./common/polyploid.png\" style=\"width:50px;height:18px\">",
		"NA" => ""
	);


	my %vals;
	my @sigs = qw/sig1 sig2 sig3 sig5 sig6 sig8 sig13 sig17 sig18 sig20 sig26/;
	for my $sig (@sigs)
	{
		if (exists $data{"csnnls_$sig"})
		{
			if ($data{"csnnls_$sig"} > 0)
			{
				$vals{$sig} = $data{"csnnls_$sig"};
			}
		}
	}
	my $cosmicString = "";
	for my $sig (sort { $vals{$b} <=> $vals{$a} } keys %vals)
	{
		$cosmicString .= "$sig,";
	}
	$cosmicString =~ s/,$//;
	$cosmicString =~ s/sig//g;



	print $html "<table style=\"width:49%;float:left\">\n";
	print $html "<tr><td colspan=\"3\"><b><h3>Somatic Mutation Load</h3></td></tr>\n";
	print $html "</table>\n";

	print $html "<table style=\"width:49%;float:right\">\n";
	print $html "<tr><td colspan=\"3\"><b><h3>Signatures and Classifications</h3></td></tr>\n";
	print $html "<tr><td><b>DSBR Hallmarks:</b></td><td>$data->{dsbr_score}/10</td><td>$classImg{dsbr_score}</td></tr>\n";
	print $html "<tr><td><b>MMR Hallmarks:</b></td><td>$data->{mmr_score}/4</td><td>$classImg{mmr_score}</td></tr>\n";
	print $html "<tr><td><b>COSMIC Signatures:</b></td><td>$cosmicString</td><td>$classImg{cosmic_score}</td><td></td></tr>\n";
	print $html "<tr><td><b>Waddell:</b></td><td>$data->{waddell}</td><td>$classImg{$data->{waddell}}</td></tr>\n";
	print $html "<tr><td><b>Collisson (RNA):</b></td><td>$data->{collisson}</td><td>$classImg{$data->{collisson}}</td></tr>\n";
	print $html "<tr><td><b>Moffitt (RNA):</b></td><td>$data->{moffitt}</td><td>$classImg{$data->{moffitt}}</td></tr>\n";
	print $html "<tr><td><b>Bailey (RNA):</b></td><td>$data->{bailey}</td><td>$classImg{$data->{bailey}}</td></tr>\n";
	print $html "<tr><td><b>Ploidy:</b></td><td>$ploidy</td><td>$classImg{$ploidyForImg}</td></tr>\n";
	print $html "<tr><td></td></tr>\n";
	print $html "<tr><td></td></tr>\n";
	print $html "<tr><td colspan=\"3\"><b><h3>CELLULOID Purity/Ploidy Solution</h3></td></tr>\n";
	print $html "<tr><td colspan=\"3\"><img src=\"$data->{plot_dir}/$data{tumour}-celluloid_contour.png\" style=\"width:100%\"></td></tr>\n";
	print $html "</table>\n";



	print $html "<table style=\"width:15%;height:68%;float:left;clear:left\">\n";
	print $html "<tr><td style=\"text-align:center\"><b>$snvText</b><br><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-snv_count-240.png\" style=\"width:100%\"></td></tr>\n";
	print $html "<tr><td style=\"text-align:center\"><b>$indelText</b><br><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-indel_count-240.png\" style=\"width:100%\"></td></tr>\n";
	print $html "<tr><td style=\"text-align:center\"><b>$svText</b><br><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-sv_count-240.png\" style=\"width:100%\"></td></tr>\n";
	print $html "<tr><td style=\"text-align:center\"><b>$neoText</b><br><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-neo_antigens-240.png\" style=\"width:100%\"></td></tr>\n";

	print $html "</table>\n";

	
	print $html "<table style=\"width:35%;height:72%;float:left\">\n";
	print $html "<tr><td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-snv_context-400x300.png\" style=\"width:100%\"></td></tr>";
	print $html "<tr><td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-indel_sv_bins-400x300.png\" style=\"width:100%\"></td></tr>";
	print $html "</table>\n";

#	print $html "<table style=\"width:34%;height:80%;float:left\">\n";



}


sub printSlideThree
{
	my $html = shift;
	my $data = shift;

	print $html "<table style=\"width:100%;height:80%\">\n";
	print $html "<tr><td style=\"valign:middle\"><img src=\"$data->{plot_dir}/$data->{tumour}-snv_vaf-no_sample.png\" style=\"width:100%\"></td></tr>\n";
	print $html "<tr><td style=\"valign:middle\"><img src=\"$data->{plot_dir}/$data->{tumour}-indel_vaf-no_sample.png\" style=\"width:100%\"><tr><td>\n";
	print $html "</table>\n";
}


sub printSlideFour
{
	my $html = shift;
	my $data = shift;

	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}.whole_genome.png\" style=\"height:80%;display:block;margin:auto\">\n";
}



sub printClassifications
{
	my $html = shift;
	my $data = shift;
	my $genes = shift;

	my $ploidy = "diploid ($data->{ploidy})";
	if ($data->{ploidy} > 2.3)
	{
		$ploidy = "polyploid ($data->{ploidy})";
	}

	my $HLAtypes = $data->{hla_types};
	$HLAtypes =~ s/\|/, /g;

	print $html "<h2>Classifications</h2>\n";

	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td><b>Alexandrov Class:</b> $data->{alexandrov_class}</td>\n";
	print $html "<td><b>Collisson Class:</b> $data->{collisson}</td>\n";
	print $html "<td><b>Ploidy:</b> $ploidy</td>\n";
	print $html "<tr>\n";

	print $html "<tr>\n";
	print $html "<td><b>Waddell Class:</b> $data->{waddell}</td>\n";
	print $html "<td><b>Moffitt Class:</b> $data->{moffitt}</td>\n";
	print $html "<td><b>Inferred Sex:</b> $data->{inferred_sex}</td>\n";
	print $html "</tr>\n";

	print $html "<tr>\n";
	print $html "<td colspan=\"3\"><b>HLA Types:</b> $HLAtypes</td>\n";
	print $html "</tr>\n";

	print $html "</table>\n";

}

sub printDsbrScore
{
	my $html = shift;
	my $data = shift;
	my $genes = shift;

	my %niceName = (
		"dsbr_snv_load" => "SNV Load >",
		"dsbr_ct_ratio" => "SNV C>T Ratio <",
		"dsbr_del4_load" => "4bp+ Deletion Load >",
		"dsbr_del4_ratio" => "4bp+ Deletion Ratio >",
		"dsbr_delsv_load" => "100-10kbp Deletion Load >",
		"dsbr_delsv_ratio" => "100-10kbp Deletion Ratio >",
		"dsbr_dup_load" => "10k-1mbp Duplication Load >",
		"dsbr_sv_load" => "Structural Variant Load >",
		"dsbr_first_hit" => "First Gene Hit(s)",
		"dsbr_second_hit" => "Second Gene Hit(s)"
	);

	my %tableCell;

	$data->{dsbr_ct_ratio} = sprintf("%.02f", $data->{dsbr_ct_ratio});
	$data->{dsbr_del4_ratio} = sprintf("%.02f", $data->{dsbr_del4_ratio});
	$data->{dsbr_delsv_ratio} = sprintf("%.02f", $data->{dsbr_delsv_ratio});


	my $yesCircle = "<img src=\"$data->{common_dir}/DSBR-circle.png\" style=\"height:40%\">";
	my $noCircle = "<img src=\"$data->{common_dir}/Empty-circle.png\" style=\"height:40%\">";

	for my $type (qw/dsbr_snv_load dsbr_del4_load dsbr_del4_ratio dsbr_delsv_load dsbr_delsv_ratio dsbr_dup_load dsbr_sv_load/)
	{
		if ($data->{$type} > $data->{"${type}_cut"})
		{
			$tableCell{$type} = "$yesCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
		}
		else
		{
			$tableCell{$type} = "$noCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
		}
	}

	for my $type (qw/dsbr_ct_ratio/)
	{
		if ($data->{$type} < $data->{"${type}_cut"})
		{
			$tableCell{$type} = "$yesCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
		}
		else
		{
			$tableCell{$type} = "$noCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
		}
	}

	for my $type (qw/dsbr_first_hit dsbr_second_hit/)
	{
		if ($data->{$type} eq "")
		{
			$tableCell{$type} = "$noCircle <b>$niceName{$type}:</b> $data->{$type}";
		}
		else
		{
			$data->{$type} =~ s/\|/, /g;
			$tableCell{$type} = "$yesCircle <b>$niceName{$type}:</b> $data->{$type}";
		}
	}


	print $html "<h3>DSBR Deficiency Hallmark Score</h3>\n";

	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-DSBR_score_bar.png\" style=\"width:100%\"></td>";

	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td>$tableCell{dsbr_snv_load}</td>\n<td>$tableCell{dsbr_delsv_load}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td>$tableCell{dsbr_ct_ratio}</td>\n<td>$tableCell{dsbr_delsv_ratio}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td>$tableCell{dsbr_del4_load}</td>\n<td>$tableCell{dsbr_dup_load}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td>$tableCell{dsbr_del4_ratio}</td>\n<td>$tableCell{dsbr_first_hit}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td>$tableCell{dsbr_sv_load}</td>\n<td>$tableCell{dsbr_second_hit}</td>\n";
	print $html "</tr>\n";

	print $html "</table>\n";

}

sub printMmrScore
{
	my $html = shift;
	my $data = shift;
	my $genes = shift;

	my %niceName = (
		"mmr_snv_load" => "SNV Load >",
		"mmr_indel_load" => "Indel Load >",
		"mmr_first_hit" => "First Gene Hit(s)",
		"mmr_second_hit" => "Second Gene Hit(s)"
	);

	my %tableCell;

	my $yesCircle = "<img src=\"$data->{common_dir}/MMR-circle.png\" style=\"height:40%\">";
	my $noCircle = "<img src=\"$data->{common_dir}/Empty-circle.png\" style=\"height:40%\">";

	for my $type (qw/mmr_snv_load mmr_indel_load/)
	{
		if ($data->{$type} > $data->{"${type}_cut"})
		{
			$tableCell{$type} = "$yesCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
		}
		else
		{
			$tableCell{$type} = "$noCircle <b>$niceName{$type} " . $data->{"${type}_cut"} . ":</b> $data->{$type}\n";
		}
	}

	for my $type (qw/mmr_first_hit mmr_second_hit/)
	{
		if ($data->{$type} eq "")
		{
			$tableCell{$type} = "$noCircle <b>$niceName{$type}:</b> $data->{$type}";
		}
		else
		{
			$data->{$type} =~ s/\|/, /g;
			$tableCell{$type} = "$yesCircle <b>$niceName{$type}:</b> $data->{$type}";
		}
	}


	print $html "<h3>MMR Deficiency Hallmark Score</h3>\n";

	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-MMR_score_bar.png\" style=\"width:100%\"></td>";

	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td>$tableCell{mmr_snv_load}</td>\n<td>$tableCell{mmr_first_hit}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td>$tableCell{mmr_indel_load}</td>\n<td>$tableCell{mmr_second_hit}</td>\n";
	print $html "</tr>\n";
	print $html "</table>\n";

}




sub printMutationLoad
{
	my $html = shift;
	my $data = shift;

	# terribly hacky!
	my $snvText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep SNVs | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $indelText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep Indel | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $svText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep SVs | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;
	my $neoText = `grep main $data->{out_path}/$data->{plot_dir}/$data->{tumour}-histogram.Rcode | grep Neo | grep bp | sed 's/.*main="//' | sed 's/",.*//'`;


	print $html "<h3>\tMutational Load</h3>\n";
	
	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td style=\"text-align:center\"><b>$snvText</b></td>";
	print $html "<td style=\"text-align:center\"><b>$indelText</b></td>";
	print $html "<td style=\"text-align:center\"><b>$svText</b></td>";
	print $html "<td style=\"text-align:center\"><b>$neoText</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{snv_count} SNVs (%PR)</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{indel_count} Indels (%PR)</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{sv_count} SVs (%PR)</b></td>";
#	print $html "<td style=\"text-align:center\"><b>$data->{neo_antigens} Neo-antigens (%PR)</b></td>";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-snv_count-240.png\" style=\"width:100%\"></td>";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-indel_count-240.png\" style=\"width:100%\"></td>";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-sv_count-240.png\" style=\"width:100%\"></td>";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-histogram-neo_antigens-240.png\" style=\"width:100%\"></td>";
	print $html "</tr>\n";
	print $html "</table>\n";
	
	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-snv_context-400x300.png\"></td>";
	print $html "<td style=\"text-align:center\"><img src=\"$data->{plot_dir}/$data->{tumour}-indel_sv_bins-400x300.png\"></td>";
	print $html "</tr>\n";
	print $html "</table>\n";

	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-snv_vaf-no_sample.png\">\n";
	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}-indel_vaf-no_sample.png\">\n";

}


sub printDriverGenes
{
	my $html = shift;
	my $genes = shift;
	my $driverGeneOrder = shift;
	my $nonsilent = shift;

	print $html "<h3>Driver Gene Status</h3>\n";
	
	print $html "<table style=\"width:100%\" class=\"sortable\">\n";
	print $html "<tr>\n";
	print $html "<th>Gene</th><th>Variant</th><th>Copy Number</th><th>AB Count</th><th>Additional Information</th>";
	print $html "</tr>\n";
	
	my $printed;
	
	for my $gene (@{ $driverGeneOrder })
	{
		$printed = 0;
		for my $v (sort keys %{ $genes->{$gene}{variants} })
		{
			unless ($genes->{$gene}{variants}{$v}{mutation_class} =~ /germline/)
			{
				unless ($genes->{$gene}{variants}{$v}{mutation_class} =~ /NA/)
				{
					printSimpleGeneRow($gene, $v, $genes, $html, "");
					$printed = 1;
				}
			}
		}
	
		for my $v (sort keys %{ $genes->{$gene}{variants} })
		{
			if ($printed == 0)
			{
				printSimpleGeneRow($gene, $v, $genes, $html, "gene only");
				$printed = 1;
			}
		}
	
	}
	print $html "</table>\n";
}

sub printWholeGenome
{
	my $html = shift;
	my $data = shift;

	print $html "<h3>Whole Genome Summary</h3>\n";
	print $html "<img src=\"$data->{plot_dir}/$data->{tumour}.whole_genome.png\" style=\"width:100%\">\n";
}

sub printPloidyAndChromo
{
	my $html = shift;
	my $data = shift;

	print $html "<h3>Ploidy and Chromothripsis</h3>\n";
	print $html "<img src=\"$data->{plot_dir}/$data{tumour}-celluloid_contour.png\" style=\"width:100%\">\n";
}

sub printFullSomaticList
{
	my $html = shift;
	my $data = shift;
	my $genes = shift;
	my $nonsilent = shift;

	print $html "<h3>Full Deleterious Somatic Variant List</h3>\n";
	
	print $html "<table style=\"width:100%\">\n";
	print $html "<tr>\n";
	print $html "<td><b>Nonsynonymous-SNVs:</b> $data->{nonsyn_count}</td>\n";
	print $html "<td><b>Frameshift-deletions:</b> $data->{del_frameshift_count}</td>\n";
	print $html "<td><b>Deletion-Gene-BPs:</b> $data->{sv_del_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Stopgain-SNVs:</b> $data->{stopgain_count}</td>\n";
	print $html "<td><b>Nonframeshift-deletions:</b> $data->{del_nonframeshift_count}</td>\n";
	print $html "<td><b>Duplication-Gene-BPs:</b> $data->{sv_dup_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Stoploss-SNVs:</b> $data->{stoploss_count}</td>\n";
	print $html "<td><b>Frameshift-Insertions:</b> $data->{ins_frameshift_count}</td>\n";
	print $html "<td><b>Inversion-Gene-BPs</b> $data->{sv_inv_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Splicing-SNVs:</b> $data->{splice_count}</td>\n";
	print $html "<td><b>Nonframeshift-Insertions:</b> $data->{ins_nonframeshift_count}</td>\n";
	print $html "<td><b>Translocation-Gene-BPs:</b> $data->{sv_tra_bp_gene_count}</td>\n";
	print $html "</tr>\n";
	print $html "<tr>\n";
	print $html "<td><b>Total Deleterious SNVs:</b> " . ($data->{nonsyn_count} + $data->{stopgain_count} + $data->{stoploss_count} + $data->{splice_count}) . "</td>\n";
	print $html "<td><b>Total Deleterious Indels:</b> " . ($data->{del_frameshift_count} + $data->{del_nonframeshift_count} + $data->{ins_frameshift_count} + $data->{ins_nonframeshift_count}) . "</td>\n";
	print $html "<td><b>Total Deleterious SV Breakpoints:</b> " . ($data->{sv_del_bp_gene_count} + $data->{sv_dup_bp_gene_count} + $data->{sv_inv_bp_gene_count} + $data->{sv_tra_bp_gene_count}) . "</td>\n";
	print $html "</tr>\n";
	print $html "</table>\n";
	
	print $html "<br>\n";
	
	print $html "<table style=\"width:100%\" class=\"sortable\">\n";
	print $html "<tr>\n";
	print $html "<th>Gene</th><th>Variant</th><th>Copy Number</th><th>AB Count</th><th>Additional Information</th>";
	print $html "</tr>\n";
	
	for my $gene (sort keys %{ $genes })
	{
		for my $v (sort keys %{ $genes->{$gene}{variants} })
		{
			if (exists $nonsilent->{ $genes->{$gene}{variants}{$v}{mutation_type} })
			{
				unless ($genes->{$gene}{variants}{$v}{mutation_class} =~ /germline/)
				{
					printSimpleGeneRow($gene, $v, $genes, $html, "");
				}
			}
		}
	}
	print $html "</table>\n";

}


sub printWholeChroms
{
	my $html = shift;
	my $data = shift;

	print $html "<h3>Chromosome Summary</h3>\n";


	for my $chr (qw/chr1 chr2 chr3 chr4 chr5 chr6 chr7-8 chr9-10 chr11-12 chr13-14 chr15-17 chr18-22 chrX-Y/)
	{
		print $html "<img src=\"$data->{plot_dir}/$data->{tumour}.whole_$chr.png\"><hr>\n";
	}
}





sub printFooter
{
	my $html = shift;
	my $data = shift;

	print $html "<div id=\"footer\">\n";
	print $html "<td style=\"width=33%;text-align:center;vertical-align:middle\"><img src=\"$data->{common_dir}/OICR_logo.jpg\"style=\"height:125px;\"></td>";
	print $html "<td style=\"width=33%;text-align:center;vertical-align:middle\"><img src=\"$data->{common_dir}/UHN_logo.png\" style=\"height:50px;\"></td>";
	print $html "<td style=\"width=33%;text-align:center;vertical-align:middle\"><img src=\"$data->{common_dir}/COMPASS_logo.png\" style=\"height:150px;\"></td>";
	print $html "</tr></table>\n";
	print $html "</div>\n";

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


#printSimpleGeneRow($gene, $variant, \%vars, $html);
sub printSimpleGeneRow
{
	my $gene = shift;
	my $v = shift;
	my $vars = shift;
	my $html = shift;
	my $geneOnlyFlag = shift;

	my $geneOnly = 0;
	if ($geneOnlyFlag eq "gene only")
	{
		$geneOnly = 1;
	}

	my ($nuc,$aa,$cons);

	my $img;

	print $html "<tr>";
	print $html "<td>$gene</td>";

	if ($geneOnly == 0)
	{
		if ($vars->{$gene}{variants}{$v}{mutation_class} eq "somatic sv")
		{
			print $html "<td><img src=\"./common/breakpoint.png\" style=\"height:18px;vertical-align:middle\"> $vars->{$gene}{variants}{$v}{mutation_type}</td>";
		}
		elsif (($vars->{$gene}{variants}{$v}{mutation_class} eq "somatic snv") or ($vars->{$gene}{variants}{$v}{mutation_class} eq "somatic indel") or ($vars->{$gene}{variants}{$v}{mutation_class} eq "germline snp") or ($vars->{$gene}{variants}{$v}{mutation_class} eq "germline indel"))
		{
			$nuc = $vars->{$gene}{variants}{$v}{nuc_context};
			$nuc =~ s/\|.*$//;		# first isoform only
			$nuc =~ s/NM.*?://;		# drop the isoform name
	
			$aa = $vars->{$gene}{variants}{$v}{aa_context};
			$aa =~ s/\|.*$//;		# first isoform only
			$aa =~ s/NM.*?://;		# drop the isoform name
	
			if (($vars->{$gene}{variants}{$v}{mutation_type} eq "stopgain") or ($vars->{$gene}{variants}{$v}{mutation_type} eq "frameshift"))
			{
				$img = "<img src=\"./common/nonsense.png\" style=\"height:18px;vertical-align:middle\">";
			}
			elsif (($vars->{$gene}{variants}{$v}{mutation_type} eq "nonsynonymous") or ($vars->{$gene}{variants}{$v}{mutation_type} eq "nonframeshift") or ($vars->{$gene}{variants}{$v}{mutation_type} eq "stoploss"))
			{
				$img = "<img src=\"./common/missense.png\" style=\"height:18px;vertical-align:middle\">";
			}
			elsif ($vars->{$gene}{variants}{$v}{mutation_type} eq "splicing")
			{
				$img = "<img src=\"./common/splicing.png\" style=\"height:18px;vertical-align:middle\">";
				$aa = $nuc;		# use nucleotide context for splicing variants
			}
			elsif ($vars->{$gene}{variants}{$v}{mutation_type} eq "altered promoter")
			{
				$img = "<img src=\"./common/noncoding.png\" style=\"height:18px;vertical-align:middle\">";
			}
			else
			{
				$img = "";
			}


			print $html "<td>$img $vars->{$gene}{variants}{$v}{mutation_type} ($aa)</td>";
		}
		elsif ($vars->{$gene}{variants}{$v}{mutation_type} eq "strong amplification")
		{
			print $html "<td><img src=\"./common/strongamp.png\" style=\"height:18px;vertical-align:middle\"> $vars->{$gene}{variants}{$v}{mutation_type}</td>";
		}
		elsif ($vars->{$gene}{variants}{$v}{mutation_type} eq "homozygous deletion")
		{
			print $html "<td><img src=\"./common/copyloss.png\" style=\"height:18px;vertical-align:middle\"> $vars->{$gene}{variants}{$v}{mutation_type}</td>";
		}
		else
		{
			print $html "<td>$vars->{$gene}{variants}{$v}{mutation_type}</td>";
		}
	}
	else
	{
		print $html "<td></td>";
	}

	my $cn = 99999;
	my $cnCount = 0;
	my $cnFlag = "";
	my $ab;

	my @abs = split(/\|/, $vars->{$gene}{ab_counts});

	for my $state (split(/\|/, $vars->{$gene}{copy_number}))
	{
		unless ($state eq "NA")
		{
			if ($state < $cn)
			{
				$cn = sprintf("%0.3f", $state);
				$ab = $abs[$cnCount];
			}
		}
		$cnCount++;
	}

	if ($cnCount > 1)
	{
		$cnFlag = "*";
	}

	if ($cn == 99999)
	{
		$cn = "NA";
		$ab = "NA";
	}

	if ($cn eq "NA")
	{
		print $html "<td style=\"text-align:center\">$cn$cnFlag ($ab)</td>";
	}
	elsif ($cn < 0.5)
	{
		print $html "<td class=\"cn_loss\" style=\"text-align:center\">$cn$cnFlag ($ab)</td>";
	}
	elsif (($cn < 1.5) or ($ab =~ /^0\./))
	{
		print $html "<td class=\"cn_loh\" style=\"text-align:center\">$cn$cnFlag ($ab)</td>";
	}
	elsif ($cn >= ($vars->{$gene}{ploidy}*4))
	{
		print $html "<td class=\"cn_very_gain\" style=\"text-align:center\">$cn$cnFlag ($ab)</td>";
	}
	elsif ($cn > ($vars->{$gene}{ploidy} + 1))
	{
		print $html "<td class=\"cn_gain\" style=\"text-align:center\">$cn$cnFlag ($ab)</td>";
	}
	else
	{
		print $html "<td style=\"text-align:center\">$cn$cnFlag ($ab)</td>";
	}


	print $html "</tr>";
	


}




