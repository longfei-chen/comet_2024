UWB=("UWB1 500 650 800 950" \
"UWB2 950 1150 1350 1550 1750" \
"UWB3 1750 1950 2150 2350 2550" \
"UWB4 2550 2700 2850 3000 3150 3300 3450")


for band in "${UWB[@]}"
do
    band_i=(${band[@]})
    N=${#band_i[@]}
    band_name=${band_i[0]}
    for idx in $(seq 1 $((N-2)))
    do
        band_begin=${band_i[$idx]}
        band_end=${band_i[$idx+1]}

        echo "Processing" ${band_name} ${band_begin} ${band_end} "..."
        ~/Applications/miniconda3/bin/python3 comet_A3_doppler.py ${band_name} ${band_begin} ${band_end}
        echo "Done"
    done
done
