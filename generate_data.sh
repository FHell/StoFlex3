mkdir -p /p/projects/coen/projects/stochasticflexibility2/exp_results/

module load julia/1.11.0

for F_index in {1..5}; do
    for f_int_index in {1..3}; do
        for n_sample_index in {1..5}; do
            (
            echo "Index of F capacity: $F_index, index of flexibility interval: $f_int_index, index of n samples $n_sample_index"
            julia ./scripts/run_experiment.jl $F_index $f_int_index $n_sample_index
            ) &

            if [[ $(jobs -r -p | wc -l) -ge 7 ]]; then
                wait -n
            fi
        done
    done
done


wait