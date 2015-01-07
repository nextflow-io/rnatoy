Toy RNA pipeline 
=================

A proof of concept of a RNA-Seq pipeline intended to show Nextflow
scripting and reproducibility capabilities.


How execute it
----------------

1) Install Docker on your computer. Read more here https://docs.docker.com/

2) Install Nextflow (version 0.10.0 or higher)

    curl -fsSL get.nextflow.io | bash

3) Pull the required Docker image as shown below 

    docker pull nextflow/rnatoy:1.1


4) Launch the pipeline execution 

    nextflow run nextflow-io/rnatoy -with-docker 
    
    
    

   
