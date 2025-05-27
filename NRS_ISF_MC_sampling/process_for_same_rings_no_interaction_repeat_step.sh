#!/bin/bash
#exec 2>/dev/null

# Check if the correct number of arguments is provided

#./process_for_same_rings_no_interaction_repeat_step.sh T6000_3000-10_time_50000 10 50 0.01 55 12



 
if [ ! -d "data/$1" ]; then
    mkdir "data/$1"
fi

if [ ! -d "data/$1/img_pairdf_fix_max_angle" ]; then
    mkdir "data/$1/img_pairdf_fix_max_angle"
fi

if [ ! -d "data/$1/img_folding_fix_max_angle" ]; then
    mkdir "data/$1/img_folding_fix_max_angle"
fi

current_dir=$(pwd)
echo "Current working directory: $current_dir"

# Print parent directory
p_parent_dir=$(dirname "$current_dir")
parent_dir_origin="$p_parent_dir/amorphous_ring_stat_analysis/data/$1"
echo "Parent directory: $parent_dir_origin"


# Assign parent directory and total_num from command line arguments

total_num=$2
sample_repeat_low=$3
max_angle=$4
step=$5
n=$6



if [ ! -d "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}" ]; then
    mkdir "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}"
fi


if [ ! -d "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/img_pairdf_fix_max_angle/" ]; then
    mkdir "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/img_pairdf_fix_max_angle/"
fi

if [ ! -d "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/img_folding_fix_max_angle/" ]; then
    mkdir "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/img_folding_fix_max_angle/"
fi



# Loop from start to end index
for ((i=0; i<=total_num; i++))
    
do  
    sample_repeat=$(printf "%.5f" $(echo "$sample_repeat_low + $i * $step" | bc -l))    
    

    echo "Copy files for sample_repeat '$sample_repeat' "
 
    if [ ! -d "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}" ]; then
        mkdir "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}"
    fi

    cp -r "$parent_dir_origin/n${n}/out/coord_real.dat" "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/coord_real.dat"
    cp  -r "$parent_dir_origin/n${n}/out/round_rough_real.dat" "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/round_rough_real.dat"

    file="input_file.dat"

    ct_val="${n}"
    awk -v key="n" -v ct_val="$n" '$1 == key {$2 = ct_val} 1' "$file" > temp && mv temp "$file"
    ct_val="${sample_repeat}"
    awk -v key="sample_repeat" -v ct_val="$sample_repeat" '$1 == key {$2 = ct_val} 1' "$file" > temp && mv temp "$file"
    ct_val="${max_angle}"
    awk -v key="max_angle" -v ct_val="$max_angle" '$1 == key {$2 = ct_val} 1' "$file" > temp && mv temp "$file"
    
    python3 gen_second_max_random_structures_step1.py "$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/"

    python3 pair_distribution_nearby_element.py "$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/"

    #python3 plot_coord.py "$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/"

    # use this for figures
    python3 plot_coord_with_real_val.py  "$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/"


    cp -r "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/pair_distribution_hist_.png" "data/$1/img_pairdf_fix_max_angle/n${n}_max_angle_${max_angle}_repeat_${sample_repeat}.png"
    cp -r "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/distance_distribution_.png" "data/$1/img_folding_fix_max_angle/n${n}_max_angle_${max_angle}_repeat_${sample_repeat}.png"

    cp -r "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/pair_distribution_hist_.png" "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/img_pairdf_fix_max_angle/n${n}_max_angle_${max_angle}_repeat_${sample_repeat}.png"
    cp -r "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/distance_distribution_.png" "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/img_folding_fix_max_angle/n${n}_max_angle_${max_angle}_repeat_${sample_repeat}.png"

    cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/skew_out.dat" >> data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/skew_out.dat
     cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/skew_out.dat" >> data/$1/skew_out.dat

    cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/direct_out.dat" >> data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/direct_out.dat
     cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/direct_out.dat" >> data/$1/direct_out.dat

    cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/percentile_out.dat" >> data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/percentile_out.dat
     cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/percentile_out.dat" >> data/$1/percentile_out.dat

    cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/bare_out.dat" >> data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/bare_out.dat
     cat "data/$1/n_${n}_sample_repeat_low_${sample_repeat_low}_step_${step}_max_angle_${max_angle}/sample_repeat${i}/bare_out.dat" >> data/$1/bare_out.dat


done