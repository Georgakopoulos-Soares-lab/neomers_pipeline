# Neomer Extraction Pipeline

[second-neomers-extraction-drawio-2.png](https://postimg.cc/SXTvhDS9)

## Summary
This pipeline orchestrates the extraction of neomers from mutation data through the following high-level steps:

- **KMC Database Preparation:**  
  The script first checks for the presence of precomputed KMC suffix and prefix files. If these files are missing, it uses `seqtk` to generate a reverse complement of the reference genome FASTA file, creates an input list of the genome files, and then runs the KMC tool to build a k-mer database.

- **Environment Setup:**  
  After ensuring the KMC database is ready, the script sources the necessary KMC path settings and prepares directories for logs and outputs.

- **Neomer Extraction:**  
  A Python CLI tool is executed with parameters such as k-mer size, data directory, MAF file, FASTA file, and output directory. This tool:
  - Reads mutation data from a MAF file.
  - Extracts reference sequences based on mutation coordinates.
  - Applies mutation logic to generate altered (mutated) sequences.
  - Generates k-mers from both original and mutated sequences.
  - Identifies potential neomers by comparing these k-mer sets.
  - Filters out k-mers that are already present in the precomputed KMC database.
  - Aggregates and logs the results from multiple parallel processing workers.



## Preparation
1. **Build seqtk:**  
```sh
cd seqtk && make
```
2. **Set up Python Environment:**
Create a virtual environment called neomer_venv, activate it, and install the required dependencies from the repository’s top-level requirements.txt:  
```sh
python3 -m venv neomer_venv
source neomer_venv/bin/activate
pip install -r requirements.txt
```

3. **Prepare Genome FASTA:**  
Move your reference FASTA file (e.g., hg19.fa) into the data_working_dir directory:

```sh
mv /path/to/hg19.fa data_working_dir/
```

## Run Extraction
```sh
chmod +x orchestration_script.sh
./orchestration_script.sh
```

