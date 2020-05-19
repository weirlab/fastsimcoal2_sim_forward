# What this pipeline does

This pipeline takes output from a fastSIMCOAL2 model (run separate from this pipeline) and simulates SNPs under that model to a variety of time points both past and future. The data is analyzed using Hudsons Fst and ADMIXTURE at each time point allowing the user to see how genome-wide Fst and admixture values are expected to change through time between the populations.

---

## Limitations

The pipeline is currently coded only for a fastSIMCOAL2 model with 5 populations. Recoding would be necessary in order to allow other numbers of populations in the model.

## Setup
The following R packages must be installed:
dplyr

---
# PIPELINE
---
## First, source the R code
```bash
source("/home/0_PROGRAMS/0_Weir_Lab_Custom_Code/fastSIMCOAL2_sim_forward_Pipeline/R_functions")
```

## USER SETS OPTIONS

```bash
   WD="/home/0_PROGRAMS/fastsimcoal2/fsc26_linux64/SIMULATE_FORWARD_5POPS/"
   
   BATCHFILE="MODEL_BATCH_5POP.par"
   
   TIMES=c(0,1,2, 3, 4, 5, 10, 100,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,
      70000,80000,90000,100000,105701,110000,120000, 130000, 140000, 150000, 160000, 170000, 180000, 190000, 200000, 
      250000, 300000, 400000, 450000, 500000, 600000, 700000, 800000, 900000, 1000000, 1500000, 2000000, 2500000, 
      3000000, 4000000, 5000000)
      
   ADD_TIME=-(105701 - TIMES) #amount of time to add to TIME1 and TIME2 (for projecting into future)
   
   MODEL_NAME=""
   
   RESULTS="SIMULATION_RESULTS"
   
   NCORES=24'
```


### USER SETS fastSIMCOAL parameters

```bash
   #Ne of current populations
   
      Ne_POP0=800000
      
      Ne_POP1=847665
      
      Ne_POP2=850487
      
      Ne_POP3=769159
      
      Ne_POP4=800000

   #note, the code below is only set up for samples sizes of 100. Changing to another value will require reprogramming.
   
      SAMP_POP0=100
      SAMP_POP1=100
      SAMP_POP2=100
      SAMP_POP3=100
      SAMP_POP4=100

   #migration rates (backwards in time: M01 in forward time is migration from 1 to 0)
   
      FACTOR=1 #leave at 1 to keep rates at their actual values
      MG01=0
      MG02=0
      MG03=0
      MG04=0
      MG10=0
      MG12=1.71978E-05*FACTOR
      MG13=0
      MG14=0
      MG20=0
      MG21=8.86588E-06*FACTOR
      MG23=1.96296E-05*FACTOR
      MG24=0
      MG30=0
      MG31=0
      MG32=5.94791E-05*FACTOR
      MG34=0
      MG40=0
      MG41=0
      MG42=0
      MG43=0

   #Ancestral Ne for ancestor of R. melanogaster and R. carbo as a ratio of the the Ne for R. melanogaster
      NewDemeSize_ANC01432=0.0094129
      NewDemeSize_ANC32=1
      NewDemeSize_ANC43=0.1569759
      NewDemeSize_ANC01=1
     
   #time of split in generations for Amazonian and Tarapoto R. carbo populations
      TIME2=105701

   #time of split in generations for R. melanogaster and R. carbo
      TIME1=508799

   #number of loci (each a single base pair) to sequence
      NLOCI=300000

   #Mutation rate
      MUTATION_RATE=3.79e-9
```

### Run Pipeline

```bash
fastSIMCOAL2_sim_forward(WD, BATCHFILE, TIMES, ADD_TIME, MODEL_NAME, RESULTS, NCORES, Ne_POP0, Ne_POP2,  Ne_POP3, Ne_POP4, SAMP_POP0, SAMP_POP1, SAMP_POP2, SAMP_POP3, SAMP_POP4, MG01, MG02, MG03, MG04, MG10, MG12, MG13, MG14, MG20, MG21, MG22, MG23, MG24, MG30, MG31, MG32, MG34, MG40, MG41, MG42, MG43, NewDemeSize_ANC01432, NewDemeSize_ANC32, NewDemeSize_ANC43, NewDemeSize_ANC01, TIME1, TIME2, NLOCI, MUTATION_RATE)
```

