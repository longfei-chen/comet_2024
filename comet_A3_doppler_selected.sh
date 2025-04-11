UWB=("UWB2 1300 1350 1350 1450 1630 1664 1664 1670 1705 1750" \
"UWB3 1750 1800 1975 2060 2060 2100 2100 2150 2235 2325" \
"UWB4 2600 2660 2700 2750 2750 2845 2870 2940 2960 3000 3050 3120 3120 3150 3150 3300 3320 3450")


for band in "${UWB[@]}"
do
    band_i=(${band[@]})
    N=${#band_i[@]}
    band_name=${band_i[0]}
    for idx in $(seq 1 2 $((N-1)))
    do
        band_begin=${band_i[$idx]}
        band_end=${band_i[$idx+1]}

        echo "Processing" ${band_name} ${band_begin} ${band_end} "..."
        ~/Applications/miniconda3/bin/python3 comet_A3_doppler.py ${band_name} ${band_begin} ${band_end}
        echo "Done"
    done
done
