name=ZTFJ1406+1222_20220525_0.fits_DM23.00
for i in {1..11}
do
   prepfold -accelcand ${i} -accelfile ${name}_ACCEL_600.cand -noxwin ${name}.dat
done
