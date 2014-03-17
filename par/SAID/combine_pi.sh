energy="290"
file1="sigma_T_$energy"
file2="T_$energy"
file3="E_$energy"
file4="F_$energy"
for i in `seq 0.0 180.0`
do
j=$(awk '$1=='$i' {print $2}' "$file1")
k=$(awk '$1=='$i' {print $2}' "$file2")
l=$(awk '$1=='$i' {print $2}' "$file3")
m=$(awk '$1=='$i' {print $2}' "$file4")
echo -e "$i\t$j\t$k\t$l\t$m" >> "Pi0_cm_$energy"
done