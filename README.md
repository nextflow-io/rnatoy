RNA-Seq demo 
=============

A toy implementation of a RNA-Seq pipeline intended to show Nextflow
scripting and reproducibility capabilities.


How execute it
----------------

1) Install Docker on your computer. Read more here https://docs.docker.com/

2) Install Nextflow (version 0.10.0 or higher)

3) Pull the required Docker image as shown below 

    docker pull nextflow/examples


4) Launch the pipeline execution 

    nextflow run nextflow-io/rnaseq-demo -with-docker 
    
    
    

   