
ls /data_zfs/low/miaocc/VT1137-0337/20220702/*_000*.fits | xargs -i python cutPol.py {}

#python fold_cal.py 8.192 0.201326592
python fold_cal.py 49.152 0.201326592


