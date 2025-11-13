#!/bin/bash
set -e

echo "Starting installation setup..."

# Check Nextflow installation
if ! command -v nextflow &> /dev/null; then
  echo "Nextflow is not installed. Please install Nextflow (>=24.10.5) (https://www.nextflow.io/docs/latest/install.html#install-nextflow) and rerun the script."
  exit 1
fi

# Check Docker installation
if ! command -v docker &> /dev/null; then
  echo "Docker is not installed. Please install Docker (>=28.0.1) (https://docs.docker.com/desktop/) and rerun the script."
  exit 1
fi

# Install required packages kallisto and gffread
echo "Installing required packages: kallisto, gffread..."
sudo apt update
sudo apt install -y kallisto gffread

# Pull required Docker images
echo "Pulling Docker images..."
docker pull johnnyfan0875/de_analysis:0.1.0
docker pull johnnyfan0875/enrichment_analysis:0.1.0

# Run kallisto index creation script
echo "Creating kallisto index..."
bash scripts/kallisto_index.sh
# alternative: kb ref (https://www.kallistobus.tools/kb_usage/kb_ref/)

echo ""
echo "Optional: To get the latest tx2gene.tsv file, run:"
echo "  Rscript scripts/generate_tx2gene.R"
echo ""
echo "Installation and setup complete!"
