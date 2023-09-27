# hocomoco_rsnp_benchmarks
A pipeline for ChIP-Seq and HT-SELEX motif benchmarking for HOCOMOCO v12.
It requires additional software to be placed in the specified directories.

## Python requirements
Your Python version must be 3.8 or higher.
### Python libraries
You must have these Python packages installed:
* NumPy
* pandas
* SciPy
* scikit-learn

## Additional software
1. [*MACRO-PERFECTOS-APE*](https://github.com/autosome-ru/macro-perfectos-ape) —
   place the file ```ape.jar``` into the directory ```./external_programs```.
2. [*Bedtools*](https://bedtools.readthedocs.io/en/latest/content/installation.html) —
   place the file ```bedtools.static``` into the directory ```./external_programs```.
3. [*SPRY-SARUS __v2.0.2__*](https://github.com/autosome-ru/sarus) —
   place the file ```sarus-2.0.2.jar``` into the directory ```./external_programs```.

## Adjusting the threads number
To increase the number of threads for computing write the exact number into ```./procfile```
without any other symbols. The default value is __1__ which means single-threaded computing.

## Usage
Execute ```./autorun.sh``` __in this very directory!__

## Input data
Ten motifs of the transcription factor FOXA2 were chosen as demonstration motifs.
These matrices are placed in the ```./pwm``` directory.

### Custom motifs
* You can benchmark your own models placing them into the ```./pwm``` directory.
* The file name must start with transcription factor name separated from the rest of the PWM name with ```@``` symbol. The extension of the file must be ```.pwm```.
* ADASTRA and HT-SELEX data for a custom transcription factor must be placed in the directories ```./adastra/TF``` and ```./selex/batchX``` respectively where ```batchX``` refers to batch1 or batch2 depending on the set of HT-SELEX experiments. The names of these files must include transcription factor name only. The extension of the files must be ```.tsv```.
* Genome files must be placed in the ```./assembly``` directory. These must be [GRCh37](https://www.ncbi.nlm.nih.gov/assembly/2758/) (hg19) and [GRCh38](https://www.ncbi.nlm.nih.gov/assembly/88331) (hg38) human genome assemblies. The file names must be ```hg19.fa``` and ```hg38.fa``` respectively.


## Output data
The output data is stored in the ```./results``` directory.
The results for ChIP-Seq both batches of HT-SELEX are written in the files ```adastra_motifs.tsv``` and ```selex_motifs.tsv``` respectively. 

## Authors
The benchmark was written by Mikhail Nikonov.
