RNA-Seq toy pipeline 
======================

A proof of concept of a RNA-Seq pipeline intended to show Nextflow
scripting and reproducibility capabilities.

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.24.0-brightgreen.svg)](http://nextflow.io)
[![Build Status](https://travis-ci.org/nextflow-io/rnatoy.svg?branch=master)](https://travis-ci.org/nextflow-io/rnatoy)

How execute it
----------------

1) Install Docker on your computer. Read more here https://docs.docker.com/

2) Install Nextflow (version 20.01.0 or later)

    `curl -fsSL get.nextflow.io | bash`

3) Pull the required Docker image as shown below: 

    `docker pull nextflow/rnatoy:1.3`


4) Launch the pipeline execution: 

    `nextflow run nextflow-io/rnatoy -with-docker` 
    
    
