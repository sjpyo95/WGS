# Name of the conda environment
ENV_NAME=$1

# Create a new conda environment named 'CHIP_analysis' with Python 3.9
conda create -n $ENV_NAME python=3.9 -y

# Activate the environment
conda activate $ENV_NAME

# Install the required packages
conda install -c conda-forge pandas tqdm -y
conda install -c bioconda fastqc trim-galore bwa picard gatk4 annovar snpeff -y

echo "Conda environment '$ENV_NAME' setup complete with all required packages."
echo "To activate the environment, use: conda activate $ENV_NAME"

conda deactivate

