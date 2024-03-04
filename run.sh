c=1
for x in `ls testcases`
do  
    echo $c
    start=$(date +%s%N)
    p=`./st-exact -s 1234 "testcases/$x"`
    end=$(date +%s%N)
    input_file="testcases/$x"
    output_file="output.csv"
    echo -n "$c," >> "$output_file"
    awk 'NR==2 {printf "%s,", $2} NR==3 {printf "%s,", $2} /Terminals / {printf "%s,", $2}' "$input_file" >> "$output_file"
    echo "$(($(($end-$start))/1000000)) ms" >> "$output_file"
    c=$((c + 1))
done