
sum=$1
vars=$2
outpath=$3

/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotIndelSVbars.pl $sum $outpath
#/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotChromothripsis.pl
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotHistograms.pl $sum $outpath
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotHallmarkBars.pl $sum $outpath
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotCosmicSignnlsBars.pl $sum $outpath
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotCelluloid.pl $sum $outpath
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotOncoSlice.pl $sum $vars $outpath
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotSNVbars.pl $sum $outpath
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotSSMfreqs.pl $sum $outpath

# last 'cause slow
/.mounts/labs/PCSI/users/rdenroche/phoenixReporter/plotWholeGenome.pl $sum $outpath

