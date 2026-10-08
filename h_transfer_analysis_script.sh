nan_if_empty() {
    [ -z "$1" ] && echo "NaN" || echo "$1"
}

selected=(
"1 7373"
"2 1121"
"2 6196"
"2 9348"
"3 680"
"3 3270"
"4 6755"
"4 8982"
"5 1508"
"5 2181"
"5 4271"
"5 9713"
"6 3250"
"6 6047"
"7 6431"
"9 3620"
"9 4761"
"9 5017"
"9 6533"
"10 602"
"10 4293"
)
for pair in "${selected[@]}"; do
    read -r i j <<< "$pair"

    dir="set$i/$j"
    cd $dir
        # REL.TRANS.ENERGY
    grep "REL.TRANS.ENERGY= " prod1.out > o
    trans_energy=$(cut -c 34-47 o)
    trans_energy=$(nan_if_empty "$trans_energy")

    # Fragment A
    grep "RESULTS FOR FRAGMENT A" -A 1 prod1.out | grep "TOTAL_ENERGY=" > frag_a
    total_energy_a=$(cut -c 30-44 frag_a)
    vib_energy_a=$(cut -c 58-72 frag_a)
    rot_energy_a=$(cut -c 86-98 frag_a)

    total_energy_a=$(nan_if_empty "$total_energy_a")
    vib_energy_a=$(nan_if_empty "$vib_energy_a")
    rot_energy_a=$(nan_if_empty "$rot_energy_a")

    # Fragment B
    grep "RESULTS FOR FRAGMENT B" -A 1 prod1.out | grep "TOTAL_ENERGY=" > frag_b
    total_energy_b=$(cut -c 30-44 frag_b)
    vib_energy_b=$(cut -c 58-72 frag_b)
    rot_energy_b=$(cut -c 86-98 frag_b)

    total_energy_b=$(nan_if_empty "$total_energy_b")
    vib_energy_b=$(nan_if_empty "$vib_energy_b")
    rot_energy_b=$(nan_if_empty "$rot_energy_b")
    grep "REACTION OCCURRED FOR PATH" prod1.out > p
    path=$(cut -c 30-32 p)
    path=$(nan_if_empty "$path")
    # Output
    echo -e "$i\t$j\t$path\t$trans_energy\t$total_energy_a\t$vib_energy_a\t$rot_energy_a\t$total_energy_b\t$vib_energy_b\t$rot_energy_b" \
        >> ../../prod_energies1.txt
    cd ../../
done
