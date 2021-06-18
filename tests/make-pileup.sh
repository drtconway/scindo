#!/bin/bash

samtools mpileup -Q 0 $1/data/tiny-bam.bam > $1/data/tiny-bam.pileup
