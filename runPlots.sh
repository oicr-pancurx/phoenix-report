
sum=$1
vars=$2
outpath=$3

prodpath=/.mounts/labs/PCSI/production/phoenix-report

${prodpath}/plotIndelSVbars.pl $sum $outpath
#${prodpath}/plotChromothripsis.pl
${prodpath}/plotHistograms.pl $sum $outpath
${prodpath}/plotHallmarkBars.pl $sum $outpath
${prodpath}/plotCosmicSignnlsBars.pl $sum $outpath
${prodpath}/plotCelluloid.pl $sum $outpath
${prodpath}/plotOncoSlice.pl $sum $vars $outpath
${prodpath}/plotSNVbars.pl $sum $outpath
${prodpath}/plotSSMfreqs.pl $sum $outpath

# last 'cause slow
${prodpath}/plotWholeGenome.pl $sum $outpath

