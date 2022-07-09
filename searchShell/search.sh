#filename=ZTFJ1406+1222_20220525.fits

#rfifind -ignorechan 0:200,3900:4096 -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -mask ${filename}*.mask -ncpus 32 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan 0:200,3900:4096 -ncpus 32 -o ${filename} ${filename}



#for i in {1..19}
#do
    #prepfold -accelcand ${i} -accelfile temp_ACCEL_0.cand -noxwin temp.dat
#done


#filename=ZTFJ1406+1222_20220525_0.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}


#filename=ZTFJ1406+1222_20220525_1.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_2.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_3.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_4.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_5.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_6.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_7.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_8.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}
#
#
#filename=ZTFJ1406+1222_20220525_9.fits
##python wide_rfi.py ${filename}
#rfilist=`python getrfitxt.py ${filename}`
#echo ${rfilist}
##rfifind -ignorechan ${rfilist} -time 0.5 -o ${filename} ${filename}
#prepsubband -nobary -lodm 20 -dmstep 1 -numdms 10 -ignorechan ${rfilist} -ncpus 32 -mask ${filename}*.mask -o ${filename} ${filename}



ls *.dat | xargs -i realfft {}

ls *.fft | xargs -i accelsearch -zmax 600 -wmax 600 -inmem -ncpus 32 {}
