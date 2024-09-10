# mol=("CH3OCHO_1" "CH3OCHO_2")
# mol=("18OH_1" "18OH_2" "18OH_3")
# mol=("OH1665" "OH1667")
# mol=("13CH3OH_1" "13CH3OH_2")
# mol=("HC5N_1" "HC5N_2" "HC9N")
# mol=("CH3264" "CH3335" "CH3349")
# mol=("CH3CHOHCH2OH")
# mol=("c-C3H_1" "c-C3H_2" "c-C3H_3" "c-C3H_4")
# mol=("H2SO4")
mol=("HI")


dt=("20240417" "20240424" "20240510" "20240511" "20240513")

for m in ${mol[@]}
do
    for d in ${dt[@]}
    do
        python comet_12P_OH_CH.py ${m} ${d}
    done
done
