# Dockerized Nextflow Pipeline for ATAC-seq and CUT&Tag Data Analysis

This repository contains a dockerized Nextflow pipeline designed to analyze either ATAC-seq or CUT&Tag Nanopore sequencing data.

## Prerequisites

Before running the pipeline, ensure you have the following files in the same directory on your machine:
- Input file: either FASTQ or BAM file
- `sequencing_summary.txt` file 
- Reference genome FASTA file
- Conda `env.yml` files needed for the workflow
- `Barplot_Chr.R` and `PieChart.R` scripts
- `params.json` file (must be adjusted based on input files, experiment type, and data path)

## Adjusting 'params.json':

Ensure that the params.json file is correctly adjusted based on your input files, experiment type, and data path. This file controls various parameters needed for the pipeline execution.

## Running the Pipeline

### 1. Install Docker

Ensure Docker is installed on your system. You can download and install Docker from [here](https://docs.docker.com/get-docker/).

### 2. Build the Docker Image

Navigate to the directory containing the `Dockerfile` and build the Docker image with the following command:

````sh
docker build -t nf-pipeline .
`````
### 3. Start the Docker Container

Once the image is built, start the Docker container using the command:
````sh
docker run -it nf-pipeline
`````
This command will enter the working directory in your Docker container.

### 4. Run the Nextflow Pipeline

Within the Docker container, you can run the Nextflow pipeline with the parameters you configured in the params.json file using the following command:
````sh
nextflow run ATAC_CUT_pipeline.nf -params-file params.json
`````
## Conclusion

By following these steps, you should be able to run the Dockerized Nextflow pipeline to analyze your ATAC-seq or CUT&Tag data efficiently. For any issues or questions, please refer to the [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) or raise an issue in this repository.

