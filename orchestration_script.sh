#!/bin/sh

# ---------------------------- CONFIGURATION ----------------------------

# Set the K-mer size
K=11

# Define the data directory
DATA_WORKING_DIRECTORY="./data_working_dir"

# Define the genome that the mutations map to
FASTA_FILE="hg19.fa"

# Define the MAF file path
MAF_FILEPATH="./test_files/sample_1000_first_final_consensus_passonly.snv_mnv_indel.icgc.public.maf"

# Define the neomers output directory
NEOMERS_OUTPUT_DIR="./neomer_results"

# Define the logs directory
LOGS_DIR="./neomer_logs"

# Define the Python script name (assumed to be in the current directory)
PYTHON_SCRIPT="neomer_extraction.py"

# Define the hg19 FASTA file and its reverse complement
HGENOME_FA="${DATA_WORKING_DIRECTORY}/${FASTA_FILE}"
HGENOME_REV_FA="${DATA_WORKING_DIRECTORY}/rev_${FASTA_FILE}"

# Define the input files list for KMC
KMC_INPUT_LIST="${DATA_WORKING_DIRECTORY}/input_files.txt"

# Define output log filenames based on K
OUT_LOG="${LOGS_DIR}/neomer_algo_${K}.out"
ERR_LOG="${LOGS_DIR}/neomer_algo_${K}.err"

# ---------------------------- FUNCTION DEFINITIONS ----------------------------

# Function to run KMC command
run_kmc() {
    echo "KMC suffix or prefix files not found. Preparing to run KMC..."

    # Generate reverse complement FASTA using seqtk
    echo "Generating reverse complement FASTA: ${HGENOME_REV_FA}..."
    ./seqtk/seqtk seq -r "${HGENOME_FA}" > "${HGENOME_REV_FA}"

    # Check if seqtk command was successful
    if [ $? -ne 0 ]; then
        echo "ERROR: seqtk command failed while generating ${HGENOME_REV_FA}."
        exit 1
    else
        echo "Reverse complement FASTA generated successfully."
    fi

    # Create the input files list with absolute paths
    echo "Creating KMC input files list: ${KMC_INPUT_LIST}..."
    echo "$(realpath "${HGENOME_FA}")" > "${KMC_INPUT_LIST}"
    echo "$(realpath "${HGENOME_REV_FA}")" >> "${KMC_INPUT_LIST}"

    # Verify that the input list was created successfully
    if [ ! -f "${KMC_INPUT_LIST}" ]; then
        echo "ERROR: Failed to create KMC input files list at ${KMC_INPUT_LIST}."
        exit 1
    fi

    echo "KMC input files list created successfully."

    # Run KMC on the list of input files
    echo "Running KMC with input list ${KMC_INPUT_LIST}..."
    kmc -b -k${K} -ci1 -v -fm "@${KMC_INPUT_LIST}" "${DATA_WORKING_DIRECTORY}/${K}mers.res" "${DATA_WORKING_DIRECTORY}/"

    # Check if KMC command was successful
    if [ $? -ne 0 ]; then
        echo "ERROR: KMC command failed. Check logs for details."
        exit 1
    else
        echo "KMC command completed successfully."
    fi
}

# ---------------------------- MAIN EXECUTION ----------------------------

# Check for the existence of KMC suffix and prefix files
SUF_FILE="${DATA_WORKING_DIRECTORY}/${K}mers.res.kmc_suf"
PRE_FILE="${DATA_WORKING_DIRECTORY}/${K}mers.res.kmc_pre"

if [ ! -f "${SUF_FILE}" ] || [ ! -f "${PRE_FILE}" ]; then
    run_kmc
else
    echo "KMC suffix and prefix files already exist. Skipping KMC command."
fi

# Source the KMC path settings
echo "Sourcing KMC path settings..."
source KMC/py_kmc_api/set_path.sh

# Ensure that the logs directory exists
mkdir -p "${LOGS_DIR}"

# Run the Python script with the parameterized arguments
echo "Running Python script: ${PYTHON_SCRIPT} with parameters:"
echo "  --K ${K}"
echo "  --data_dir ${DATA_WORKING_DIRECTORY}"
echo "  --maf_filepath ${MAF_FILEPATH}"
echo "  --fasta_file ${FASTA_FILE}"
echo "  --neomers_output_dir ${NEOMERS_OUTPUT_DIR}"

python3 "${PYTHON_SCRIPT}" \
    --K "${K}" \
    --data_dir "${DATA_WORKING_DIRECTORY}" \
    --maf_filepath "${MAF_FILEPATH}" \
    --fasta_file "${FASTA_FILE}" \
    --neomers_output_dir "${NEOMERS_OUTPUT_DIR}" \
    1> "${OUT_LOG}" 2> "${ERR_LOG}"

# Check if Python script executed successfully
if [ $? -ne 0 ]; then
    echo "ERROR: Python script '${PYTHON_SCRIPT}' failed. Check '${ERR_LOG}' for details."
    exit 1
else
    echo "Python script '${PYTHON_SCRIPT}' executed successfully. Output logged to '${OUT_LOG}' and '${ERR_LOG}'."
fi

# Final message
echo "Script execution completed successfully."
