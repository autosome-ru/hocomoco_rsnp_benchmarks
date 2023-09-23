ls -1 ./motifs/pwm/ > ./scripts/tf_list.txt
ls -1R ./motifs/pwm | grep "\.pwm$" > motifs_list.txt
python3 ./scripts/ape_list.py
python3 ./scripts/selex_benchmark_list.py
sed 's/\/selex/\/adastra/g' ./scripts/selex_benchmark_list.txt ./scripts/adastra_benchmark_list.txt
cat ./scripts/ape_list.txt | parallel -j ./scripts/procfile &> ./logs/ape.log
cat ./scripts/selex_benchmark_list.txt | parallel -j ./scripts/procfile > ./logs/selex_messages.log 2> ./logs/selex_errors.log
cat ./scripts/adastra_benchmark_list.txt | parallel -j ./scripts/procfile > ./logs/adastra_messages.log 2> ./logs/adastra_errors.log
mv ./scripts/*tsv ./results
