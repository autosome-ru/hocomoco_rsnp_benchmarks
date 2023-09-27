ls -1 ./motifs/pwm/ > ./scripts/tf_list.txt
ls -1R ./motifs/pwm | grep "\.pwm$" > ./scripts/motifs_list.txt
cd ./scripts
python3 ./ape_list.py
python3 ./selex_benchmark_list.py
sed 's/\/selex/\/adastra/g' ./selex_benchmark_list.txt > ./adastra_benchmark_list.txt
cat ./ape_list.txt | parallel -j ../procfile &> ../logs/ape.log
cat ./selex_benchmark_list.txt | parallel -j ../procfile > ../logs/selex_messages.log 2> ../logs/selex_errors.log
cat ./adastra_benchmark_list.txt | parallel -j ../procfile > ../logs/adastra_messages.log 2> ../logs/adastra_errors.log
cd ..
mv ./scripts/*tsv ./results
