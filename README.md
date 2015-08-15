RNA-Seq toy pipeline 
======================

A proof of concept of a RNA-Seq pipeline intended to show Nextflow
scripting and reproducibility capabilities.


This version of the pipeline was used to assess the impact of Docker containers on the 
performance of genomic tools. See https://peerj.com/preprints/1171/ and subsequent manuscript 
peer review https://peerj.com/manuscripts/5515/

To replicate the pipeline execution see the following section. 

Quick start
-----------

Make sure you have the required dependencies listed in the last section. 

1. Pull the required Docker image: 

    `$ docker pull nextflow/rnatoy:peerj5515`

2. Download the benchmark dataset: 

    `$ aws s3 sync s3://cbcrg-eu/rnatoy-mouse/ data`

3. Launch the pipeline execution with the benchmark dataset
	
    `$ nextflow run rnatoy -revision peerj5515 -with-docker -bg > log.txt`

The pipeline will run as a background process and the output redirected to the `log.txt` file. 

It will also create a file named `trace.csv` containing the run time information for each executed task.

Note: the execution with Docker will produce some intermediate results with *root* ownership. You will need sudo/root permissions to delete them.


Dependencies
------------
 * Java VM 7 or later
 * [Docker 1.0 or later](http://www.docker.io)
 * [Nextflow 0.15.0 or later](http://nextflow.io)
 
For the execution of the pipeline without using the Docker engine follow the installation steps in the included Dockerfile 
    
    

   
