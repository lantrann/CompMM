#!/bin/bash
#PBS -N prog_name
#PBS -j oe
#PBS -m abe
#PBS -M tnlan@hcmip.vast.vn
#PBS -l select=1:ncpus=1:mem=2G 
#PBS -q para_cpu

/datausers/hcmip/tnlan/MMC/ExchangeBias/excbias > /datausers/hcmip/tnlan/MMC/ExchangeBias/test.out3

