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
  #Setup RESULTS matrix
      RESULTS <- matrix(NA, nrow=length(ADD_TIME), ncol=36)
      colnames(RESULTS) <- c("ADD_TIME", "BEST_K", "ADMIX_MEAN_POP0", "ADMIX_MEAN_POP1","ADMIX_MEAN_POP2","ADMIX_MEAN_POP3","ADMIX_MEAN_POP4",
         "ADMIX_VAR_POP0", "ADMIX_VAR_POP1","ADMIX_VAR_POP2","ADMIX_VAR_POP3","ADMIX_VAR_POP4",
         "Fst_POP0_POP1", "Fst_POP0_POP2", "Fst_POP0_POP3", "Fst_POP0_POP4", 
         "Fst_POP1_POP2", "Fst_POP1_POP3", "Fst_POP1_POP4", 
         "Fst_POP2_POP3", "Fst_POP2_POP4", 
         "Fst_POP3_POP4", "BREAK", 
         "ADMIX_MEAN_POP1__POP1&2", "ADMIX_MEAN_POP2__POP1&2","BREAK","ADMIX_MEAN_POP1__POP1&2&3","ADMIX_MEAN_POP2__POP1&2&3",
         "ADMIX_MEAN_POP3__POP1&2&3", "BREAK", 
         "ADMIX_VAR_POP1__POP1&2", "ADMIX_VAR_POP2__POP1&2","BREAK","ADMIX_VAR_POP1__POP1&2&3","ADMIX_VAR_POP2__POP1&2&3",
         "ADMIX_VAR_POP3__POP1&2&3" 
   )
      RESULTS[,"ADD_TIME"] <- ADD_TIME

   #loop through TIMES
      for(z in 1:length(ADD_TIME)){
        ADD_TIME_z = ADD_TIME[z]

#loop through TIMES
for(z in 1:length(ADD_TIME)){
   ADD_TIME_z = ADD_TIME[z]

###STEP0 fastSIMCOAL2 within R
   ###Setup the .par file for input into fastSIMCOAL2
   setwd(WD)
   BATCH1 <- readLines(BATCHFILE)
   BATCH2  <- gsub(pattern = "TIME2", replace = format(TIME2+ADD_TIME_z, scientific=FALSE), x = BATCH1)
   BATCH3  <- gsub(pattern = "TIME1", replace = format(TIME1+ADD_TIME_z, scientific=FALSE), x = BATCH2)
   BATCH4  <- gsub(pattern = "Ne_POP0", replace = format(Ne_POP0, scientific=FALSE), x = BATCH3)
   BATCH5  <- gsub(pattern = "Ne_POP1", replace = format(Ne_POP1, scientific=FALSE), x = BATCH4)
   BATCH6  <- gsub(pattern = "Ne_POP2", replace = format(Ne_POP2, scientific=FALSE), x = BATCH5)
   BATCH7  <- gsub(pattern = "Ne_POP3", replace = format(Ne_POP3, scientific=FALSE), x = BATCH6)
   BATCH8  <- gsub(pattern = "Ne_POP4", replace = format(Ne_POP4, scientific=FALSE), x = BATCH7)
   BATCH9  <- gsub(pattern = "SAMP_POP0", replace = format(SAMP_POP0, scientific=FALSE), x = BATCH8)
   BATCH10  <- gsub(pattern = "SAMP_POP1", replace = format(SAMP_POP1, scientific=FALSE), x = BATCH9)
   BATCH11  <- gsub(pattern = "SAMP_POP2", replace = format(SAMP_POP2, scientific=FALSE), x = BATCH10)
   BATCH12  <- gsub(pattern = "SAMP_POP3", replace = format(SAMP_POP3, scientific=FALSE), x = BATCH11)
   BATCH13  <- gsub(pattern = "SAMP_POP4", replace = format(SAMP_POP4, scientific=FALSE), x = BATCH12)
   BATCH14  <- gsub(pattern = "MG01", replace = format(MG01, scientific=FALSE), x = BATCH13)
   BATCH15  <- gsub(pattern = "MG02", replace = format(MG02, scientific=FALSE), x = BATCH14)
   BATCH16  <- gsub(pattern = "MG03", replace = format(MG03, scientific=FALSE), x = BATCH15)
   BATCH17  <- gsub(pattern = "MG04", replace = format(MG04, scientific=FALSE), x = BATCH16)
   BATCH18  <- gsub(pattern = "MG10", replace = format(MG10, scientific=FALSE), x = BATCH17)
   BATCH19  <- gsub(pattern = "MG12", replace = format(MG12, scientific=FALSE), x = BATCH18)
   BATCH20  <- gsub(pattern = "MG13", replace = format(MG13, scientific=FALSE), x = BATCH19)
   BATCH21  <- gsub(pattern = "MG14", replace = format(MG14, scientific=FALSE), x = BATCH20)
   BATCH22  <- gsub(pattern = "MG20", replace = format(MG20, scientific=FALSE), x = BATCH21)
   BATCH23  <- gsub(pattern = "MG21", replace = format(MG21, scientific=FALSE), x = BATCH22)
   BATCH24  <- gsub(pattern = "MG23", replace = format(MG23, scientific=FALSE), x = BATCH23)
   BATCH25  <- gsub(pattern = "MG24", replace = format(MG24, scientific=FALSE), x = BATCH24)
   BATCH26  <- gsub(pattern = "MG30", replace = format(MG30, scientific=FALSE), x = BATCH25)
   BATCH27  <- gsub(pattern = "MG31", replace = format(MG31, scientific=FALSE), x = BATCH26)
   BATCH28  <- gsub(pattern = "MG32", replace = format(MG32, scientific=FALSE), x = BATCH27)
   BATCH29  <- gsub(pattern = "MG34", replace = format(MG34, scientific=FALSE), x = BATCH28)
   BATCH30  <- gsub(pattern = "MG40", replace = format(MG40, scientific=FALSE), x = BATCH29)
   BATCH31  <- gsub(pattern = "MG41", replace = format(MG41, scientific=FALSE), x = BATCH30)
   BATCH32  <- gsub(pattern = "MG42", replace = format(MG42, scientific=FALSE), x = BATCH31)
   BATCH33  <- gsub(pattern = "MG43", replace = format(MG43, scientific=FALSE), x = BATCH32)
   BATCH34  <- gsub(pattern = "NewDemeSize_ANC01432", replace = format(NewDemeSize_ANC01432, scientific=FALSE), x = BATCH33)
   BATCH35  <- gsub(pattern = "NewDemeSize_ANC32", replace = format(NewDemeSize_ANC32, scientific=FALSE), x = BATCH34)
   BATCH36  <- gsub(pattern = "NewDemeSize_ANC01", replace = format(NewDemeSize_ANC01, scientific=FALSE), x = BATCH35)
   BATCH37  <- gsub(pattern = "NewDemeSize_ANC43", replace = format(NewDemeSize_ANC43, scientific=FALSE), x = BATCH36)
   BATCH38  <- gsub(pattern = "NLOCI", replace = format(NLOCI, scientific=FALSE), x = BATCH37)
   BATCH39  <- gsub(pattern = "MUTATION_RATE", replace = format(MUTATION_RATE, scientific=FALSE), x = BATCH38)

   writeLines(BATCH39, con=paste("MODEL__", format(ADD_TIME_z, scientific=FALSE), MODEL_NAME,".par", sep=""))

   ###run fastSIMCOAL2
   PAR_WD = WD
   PAR_FILE = paste("MODEL__", format(ADD_TIME_z, scientific=FALSE), MODEL_NAME, ".par", sep="")
   cmd1 <- paste("cd ", PAR_WD, "; ./fsc26 -i ", PAR_FILE, " -q -c 20 -n 1 -s 10000", sep="")
   system(cmd1)

###STEP1 read fastSIMCOAL2 simulated SNPs into R and transform
   INPUT=paste(WD, "MODEL__", format(ADD_TIME_z, scientific=FALSE), MODEL_NAME, "/MODEL__", 
      format(ADD_TIME_z, scientific=FALSE), MODEL_NAME, "_1_1.arp", sep="")
   ARL <- arlequinRead(INPUT)
   GENOTYPES2 <- ARL$data.info[[1]][4][1]
   ggg <- as.matrix(GENOTYPES2)
   fff <- strsplit(ggg, "")
   eee <- do.call(rbind, fff)
   ddd <- matrix(as.numeric(eee), nrow=nrow(eee), ncol=ncol(eee))

   #This gives 1 and 0 for two alleles of a bialleleic snp and deletes non biallelic snps
   for(i in ncol(eee):1){
      jjj <- unique(eee[,i])
      ddd[eee[,i] == jjj[1],i] <- 1
      ddd[eee[,i] == jjj[2],i] <- 0
      if(length(jjj) != 2){ddd <- ddd[,-i]} # delete if not biallelic
   }

   #puts SNPs into 0,1,2  for homo, hetero, homozygote states for use latter in Hudson's Fst etc
      sss <- matrix(NA, nrow=nrow(ddd)/2, ncol=ncol(ddd))
      for(i in 1:50){
          sss[i,]<- as.vector(ddd[i,])+ as.vector(ddd[i+50,])
          sss[i+50,]<- as.vector(ddd[i+100,])+ as.vector(ddd[i+150,])
          sss[i+100,]<- as.vector(ddd[i+200,])+ as.vector(ddd[i+250,])
          sss[i+150,]<- as.vector(ddd[i+300,])+ as.vector(ddd[i+350,])
          sss[i+200,]<- as.vector(ddd[i+400,])+ as.vector(ddd[i+450,])
      }

###STEP2 run ADMIXTURE from within R
   WD2 <- paste(WD, "MODEL__", format(ADD_TIME_z, scientific=FALSE), MODEL_NAME, sep="")
   setwd(WD2)

   ###STEP2A Make Diploid individuals out of haploids in ped format which is 2 2 or 1 1 and haploid which is 1 2 or 2 1. Missing is 0 0
      ttt <- matrix(NA, nrow=nrow(ddd)/2, ncol=ncol(ddd)*2)
      for(i in 1:50){
          ttt[i,]<- c(rbind(as.vector(ddd[i,]), as.vector(ddd[i+50,])))
          ttt[i+50,]<- c(rbind(as.vector(ddd[i+100,]), as.vector(ddd[i+150,])))
          ttt[i+100,]<- c(rbind(as.vector(ddd[i+200,]), as.vector(ddd[i+250,])))
          ttt[i+150,]<- c(rbind(as.vector(ddd[i+300,]), as.vector(ddd[i+350,])))
          ttt[i+200,]<- c(rbind(as.vector(ddd[i+400,]), as.vector(ddd[i+450,])))
      }
      for(i in 1:nrow(ttt)){
         ttt[i,ttt[i,]==0] <- 2
      }

   ###STEP2B write output as PLINK/PED files
      ccc <- matrix(NA, nrow=nrow(ttt), ncol=6)
      VEC1 <- c(1:nrow(ccc))
      ccc[,1] <- "0"
      ccc[,2] <- paste("INDIVIDUAL_", VEC1, sep="")
      ccc[,3] <- "0"
      ccc[,4] <- "0"
      ccc[,5] <- "0"
      ccc[,6] <- "1"
      bbb <- cbind(ccc,ttt)
 
     write.table(bbb, "SIMULATED_SNPS_PLINK.ped", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
     write.table(bbb[51:150,], "SIMULATED_SNPS_PLINK_pop1and2.ped", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
     write.table(bbb[51:200,], "SIMULATED_SNPS_PLINK_pop1and2and3.ped", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)

     #STEP2C write the PLINK.map file
      aaa <- matrix("0", nrow=ncol(ttt)/2, ncol=4)
      VEC2 <- c(1:(ncol(ttt)/2))
      aaa[,1] <- "1" #Chromosome
      aaa[,2] <- paste(VEC2, "__", 1, sep="") #paste("SNP_", VEC2, sep="")
      aaa[,3] <- "0" #Genetic_distance
      aaa[,4] <- "0" #Physical_position
      write.table(aaa, "SIMULATED_SNPS_PLINK.map", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
      write.table(aaa, "SIMULATED_SNPS_PLINK_pop1and2.map", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
      write.table(aaa, "SIMULATED_SNPS_PLINK_pop1and2and3.map", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)

   #STEP2D.1 run ADMIXTURE with K=1to5 to test for best K
      PLINKPED=paste(WD2, "/SIMULATED_SNPS_PLINK.ped", sep="")
      cmd=(paste("cd ", WD2, " ; for K in 1 2 3 4 5; \ do /opt/POPGEN/admixture_linux-1.3.0/admixture -j20 --cv ", 
          PLINKPED, " $K | tee ", WD2, "/log${K}.out; done", sep=""))
      system(cmd)

      x <- readLines("log1.out")
      K1 <- as.numeric(strsplit(x[length(x)-1], ": ")[[1]][2])

      x <- readLines("log2.out")
      K2 <- as.numeric(strsplit(x[length(x)-1], ": ")[[1]][2])

      x <- readLines("log3.out")
      K3 <- as.numeric(strsplit(x[length(x)-1], ": ")[[1]][2])

      x <- readLines("log4.out")
      K4 <- as.numeric(strsplit(x[length(x)-1], ": ")[[1]][2])

      x <- readLines("log5.out")
      K5 <- as.numeric(strsplit(x[length(x)-1], ": ")[[1]][2])

      CV_K <- c(K1, K2, K3, K4, K5)
      K=c(1,2,3,4,5)
      RESULTS[,"BEST_K"] <- K[CV_K==min(CV_K)] #best K


   #STEP2D run ADMIXTURE with K=2
      #cmd=paste("cd ", WD2, "; /opt/POPGEN/admixture_linux-1.3.0/admixture --cv ", PLINKPED, " 2", sep="")
      #system(cmd)

   ###STEP2F get stats on ADMIXTURE k=2 output
      ADMIX1 <- read.table("SIMULATED_SNPS_PLINK.2.Q", header=FALSE)
      #par(mfrow = c(2, 2)) 
      #plot(ADMIX1[,1], main="Populations 0 to 4")

      #mean POP0 melanogaster "Allopatric"
      RESULTS[z,"ADMIX_MEAN_POP0"] <- mean(ADMIX1[1:50,1])

      #mean POP1 melanogaster Mayo Valley
      RESULTS[z,"ADMIX_MEAN_POP1"] <- mean(ADMIX1[51:100,1])

      #mean POP2 Tarapoto carbo
      RESULTS[z,"ADMIX_MEAN_POP2"] <- mean(ADMIX1[101:150,1])

      #mean POP3 Amazonian carbo
      RESULTS[z,"ADMIX_MEAN_POP3"] <- mean(ADMIX1[151:200,1])

      #mean POP4 "Allopatric" carbo
      RESULTS[z,"ADMIX_MEAN_POP4"] <- mean(ADMIX1[201:250,1])

      #var POP0 melanogaster "Allopatric"
      RESULTS[z,"ADMIX_VAR_POP0"] <- var(ADMIX1[1:50,1])

      #var POP1 melanogaster Mayo Valley
      RESULTS[z,"ADMIX_VAR_POP1"] <- var(ADMIX1[51:100,1])

      #var POP2 Tarapoto carbo
      RESULTS[z,"ADMIX_VAR_POP2"] <- var(ADMIX1[101:150,1])

      #var POP3 Amazonian carbo
      RESULTS[z,"ADMIX_VAR_POP3"] <- var(ADMIX1[151:200,1])

      #var POP4 "Allopatric" carbo
      RESULTS[z,"ADMIX_VAR_POP4"] <- var(ADMIX1[201:250,1])

############OPTIONAL: rerun ADMIXTURE without population 0 and 3/4 (i.e. the allopatric pops)
   #STEP2D run ADMIXTURE with K=2
      PLINKPED=paste(WD2, "/SIMULATED_SNPS_PLINK_pop1and2.ped", sep="")
      cmd=paste("cd ", WD2, "; /opt/POPGEN/admixture_linux-1.3.0/admixture --cv ", PLINKPED, " 2", sep="")
      system(cmd)

   ###STEP2E get stats on ADMIXTURE k=2 output
      ADMIX2 <- read.table("SIMULATED_SNPS_PLINK_pop1and2.2.Q", header=FALSE)
      #plot(ADMIX2[,1], "Populations 1 and 2 excluding 0 and 3")

      RESULTS[z,"ADMIX_MEAN_POP1__POP1&2"] <- mean(ADMIX2[1:50,1])
      RESULTS[z,"ADMIX_MEAN_POP2__POP1&2"] <- mean(ADMIX2[51:100,1])

      RESULTS[z,"ADMIX_VAR_POP1__POP1&2"] <- var(ADMIX2[1:50,1])
      RESULTS[z,"ADMIX_VAR_POP2__POP1&2"] <- var(ADMIX2[51:100,1])

############OPTIONAL: rerun ADMIXTURE with pops 1 2 and 3
   #STEP2D run ADMIXTURE with K=2
      PLINKPED=paste(WD2, "/SIMULATED_SNPS_PLINK_pop1and2and3.ped", sep="")
      cmd=paste("cd ", WD2, "; /opt/POPGEN/admixture_linux-1.3.0/admixture --cv ", PLINKPED, " 2", sep="")
      system(cmd)

   ###STEP2E get stats on ADMIXTURE k=2 output
      ADMIX3 <- read.table("SIMULATED_SNPS_PLINK_pop1and2and3.2.Q", header=FALSE)
      #plot(ADMIX3[,1], "Allopatric Populations 0 and 4 excluding pops 1 and 2")


      RESULTS[z,"ADMIX_MEAN_POP1__POP1&2&3"] <- mean(ADMIX3[1:50,1])
      RESULTS[z,"ADMIX_MEAN_POP2__POP1&2&3"] <- mean(ADMIX3[51:100,1])
      RESULTS[z,"ADMIX_MEAN_POP3__POP1&2&3"] <- mean(ADMIX3[101:150,1])

      RESULTS[z,"ADMIX_VAR_POP1__POP1&2&3"] <- var(ADMIX3[1:50,1])
      RESULTS[z,"ADMIX_VAR_POP2__POP1&2&3"] <- var(ADMIX3[51:100,1])
      RESULTS[z,"ADMIX_VAR_POP3__POP1&2&3"] <- var(ADMIX3[101:150,1])


############
###Calculate Hudson's Fst
DATA <- sss

POPS <- vector(mode = "numeric", length = nrow(sss))
POPS[1:50] <- 1
POPS[51:100] <- 2
POPS[101:150] <- 3
POPS[151:200] <- 4
POPS[201:250] <- 5

NAMES <- vector(mode = "numeric", length = nrow(sss))
for(i in 1:nrow(sss)){NAMES[i] <- i}

DATA2 <- cbind(NAMES, POPS, DATA)

#NBOOT = 200
#RESULT <- Pairwise_Hudsons_Fst(DATA2, BOOT = TRUE, NBOOT, NCORES)
FST <- Pairwise_Hudsons_Fst(DATA2, BOOT = FALSE, NBOOT, NCORES)
FST

      RESULTS[z,"Fst_POP0_POP1"] <- FST[2,1]
      RESULTS[z,"Fst_POP0_POP2"] <- FST[3,1]
      RESULTS[z,"Fst_POP0_POP3"] <- FST[4,1]
      RESULTS[z,"Fst_POP0_POP4"] <- FST[5,1]

      RESULTS[z,"Fst_POP1_POP2"] <- FST[3,2]
      RESULTS[z,"Fst_POP1_POP3"] <- FST[4,2]
      RESULTS[z,"Fst_POP1_POP4"] <- FST[5,2]

      RESULTS[z,"Fst_POP2_POP3"] <- FST[4,3]
      RESULTS[z,"Fst_POP2_POP4"] <- FST[5,3]

      RESULTS[z,"Fst_POP3_POP4"] <- FST[5,4]

write.table(RESULTS, "SIMULATION_RESULTS")


print("###############################")
print(z)
print("###############################")
}
```

