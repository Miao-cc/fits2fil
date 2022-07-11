rootpath=/data_zfs/fast/miaocc/ZD2021_5/VT1137-0337/20220702
fitspath=${rootpath}/OffCal
gaindiff=1.3971754587448721
baselindiff=-0.36975335723765035
outfile=${rootpath}/PRESTO/VT1137-0337_20220702
ra=19:59:36.76
dec=20:48:15.12
cd ${fitspath}
ls *.fits > filelist.txt
#python psrfits2fil_combine.py -o ${outfile} -g ${gaindiff} -b ${baselindiff} --ra ${ra} --dec ${dec} --noweights --noscales --nooffsets filelist.txt
python psrfits2fil_combine.py -o ${outfile} -g ${gaindiff} -b ${baselindiff} --noweights --noscales --nooffsets filelist.txt
