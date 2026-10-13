file="NEB_OCCUP1.txt"
test -f PART1.txt && rm PART1.txt
test -f COMBINEPART751.txt && rm COMBINEPART751.txt
echo "Content of the file $file:"
cat "$file"

line_count=0

exec 5< "$file"  # Open file descriptor 3 to read from the file

while read -r line <&5; do
    echo "Processing line: $line"

    echo -e "$line" > occupation_in_1NNN1.txt
    echo -e "$line" 
##### Generate unique filename for each output file
   # output_file="output_files/occupation_in_1NNN1_${line_count}.txt"
    
    ## call F90 code to generate two occupation files
    #./print_occupation_files.x
    ./H751.x
    #./print_occupation_files.x

    ## perform NEB and get the barrier
    initial_file="array${line_count}_initial.xyz"
    final_file="array${line_count}_final.xyz"
    ## echo ${initial_file}
     ./HH751.x ${initial_file} ${final_file} > PART1.txt
    #./nNN.x ${initial_file} ${final_file} > barrier_data7.txt
    ## exit_code=$?
    tail -n 1 PART1.txt >> COMBINEPART751.txt
    #if [ $exit_code -ne 0 ]; then
        #echo "Barrier value: $barrier"
        #echo $barrier 
        ##echo $barrier > barrier_data1.txt
    #    tail -n 2 barrier_data1.txt >> all_barrier.txt
    #else
    #    echo "Error executing nNN.x program. Exit code: $exit_code"
    #fi

    line_count=$((line_count + 1))
    if [ "$line_count" -eq 5 ]; then
        break
    fi
done

exec 2048<&-  # Close file descriptor 2048

# Zip all the output files
#zip -r output_files.zip "$output_dir"



