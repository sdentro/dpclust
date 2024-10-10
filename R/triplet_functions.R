
#########################################################################################################
# Pipeline functions
#########################################################################################################

#' Step 1 of the triplet pipeline : run DPC core algorithm
#'
#' @param samplename Identifier to be used to name output files
#' @param no.iters The total number of iterations the DP should be run for. Burin is determined automatically
#' @param output_folder Directory where the output is to be stored
#' @param max.considered.clusters The maximum to be considered clusters (Default: 30)
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
run_dpc_core_triplet = function(samplename, no.iters, output_folder, max.considered.clusters=30) {

dat = read_triplet_pipeline_input()
attach(dat)

if(!file.exists(output_folder)){ dir.create(output_folder) }

DPClust:::DirichletProcessClustering(mutCount = mutCount,
                                     WTCount = WTCount,
                                     totalCopyNumber = totalCopyNumber,
                                     copyNumberAdjustment = copyNumberAdjustment,
                                     mutation.copy.number = mutation.copy.number,
                                     cellularity = cellularity,
                                     output_folder = output_folder,
                                     no.iters = no.iters,
                                     no.iters.burn.in = ceiling(no.iters/5),
                                     subsamplesrun = subsamplenames,
                                     samplename=samplename,
                                     conc_param = 1,
                                     cluster_conc = 5,
                                     mut.assignment.type = 1,
                                     most.similar.mut = NA,
                                     mutationTypes = NA,
                                     max.considered.clusters=max.considered.clusters)
}

#' Step 2 of the triplet pipeline : estimate density by samples grouped in triplets
#'
#' @param donorname Identifier to be used to name output files, should be the same as "samplename" in run_dpc_core_triplet
#' @param run Not needed parameter, should really always be set to 1. This code could be run across multiple donors, but has not been adapted to continue to work as such
#' @param subsample.indices The indices related to the samples to be included in this particular triplet
#' @param noiters The number of total iterations, should be the same as "no.iters" in run_dpc_core_triplet
#' @param burn.in Number of iterations to be used for burnin. This should be the same value as in run_dpc_core_triplet, in which it is defined as ceiling(no.iters/5)
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
run_parallel_density_est = function(donorname, run, subsample.indices, noiters, burn.in) {

gsdata_file = file.path("DP_output", paste0(donorname, "_gsdata.RData"))

conc_param=1
cluster_conc = 5
density.smooth = 0.01

dat = read_triplet_pipeline_density_input(donorname, run, subsample.indices)
attach(dat)

no.muts = nrow(mutation.copy.number)

load(gsdata_file)
attach(GS.data)
#S.i=GS.data$S.i
#V.h=GS.data$V.h
#pi.h=GS.data$pi.h
no.subsamples = length(subsample.indices)
pi.h=pi.h[,,subsample.indices]
GS.data = list(S.i=S.i,V.h=V.h,pi.h=pi.h)

opts = list(samplename=samplename,
            subsamplenames=subsamples$V1[subsample.indices],
            no.iters=noiters,
            no.iters.burn.in=burn.in,
            outdir= "md_out")
res = DPClust:::multiDimensionalClustering(mutation.copy.number=mutation.copy.number,
                           copyNumberAdjustment=copyNumberAdjustment,
                           GS.data=GS.data,
                           density.smooth=density.smooth,
                           opts=opts,
                           outfilesuffix=paste0("_",paste(subsample.indices,collapse="_")))
}

#########################################################################################################
# Helper functions
#########################################################################################################

#' Determines the total number of triplets from the given number of samples
#'
#' @param num_samples Integer total number of samples
#' @return Integer total of triplets
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
determine_num_triplets = function(num_samples) {
	return(choose(num_samples,3))
}

#' Determines the indices of the samples to be used for a certain triplet
#'
#' @param choose.index 
#' @param choose.number The number within the total triplets to be picked
#' @param choose.from The total number of triplets
#' @return Vector of integers
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
get.subsample.indices<-function(choose.index,choose.number,choose.from){
        subsample.indices = 1:choose.number
        temp.choose.index = choose.index
        for(i in 1:choose.number){
                last.subtotal = 0
                subtotal = choose(choose.from-subsample.indices[i],choose.number-i)
                while(temp.choose.index>subtotal){
                        subsample.indices[i] = subsample.indices[i] + 1
                        last.subtotal = subtotal
                        subtotal = subtotal + choose(choose.from-subsample.indices[i],choose.number-i)
                }
                if(i<choose.number){
                        subsample.indices[i+1] = subsample.indices[i]+1
                }
                temp.choose.index = temp.choose.index - last.subtotal
        }
        return(subsample.indices)
}

#########################################################################################################
# Data parsers
#########################################################################################################
read_triplet_pipeline_combine_input = function(donorname, run) {

sample=donorname
samplenames=donorname
subsamples=read.table("samplenames.lst")
subsamples[[run]]=as.vector(subsamples$V1)

#get cellularities
ce=read.table("purity.lst")
cellularities = list()
cellularities[[run]]=as.vector(ce$V1)

samplename = samplenames[run]
cellularity = cellularities[[run]]

no.subsamples = length(subsamples[[run]])
choose.number = 3
choose.from = no.subsamples
no.perms = choose(no.subsamples,3)

new_output_folder = "md_out"

print(paste(no.subsamples[[run]]))
#required for mutation.copy.number and copyNumberAdjustment
data=list()
for(s in 1:no.subsamples){
  data[[s]]=read.table(paste(subsamples[[run]][s],"_allDirichletProcessInfo.txt",sep=""),sep="\t",header=T)
}

#required for plotting the PDF: mutation.copy.number/copyNumberAdjustment
WTCount = data[[1]]$WT.count
mutCount = data[[1]]$mut.count
totalCopyNumber = data[[1]]$subclonal.CN
copyNumberAdjustment = data[[1]]$no.chrs.bearing.mut
normalCopyNumber = data[[1]]$normalCopyNumber
mutation.copy.number = data[[1]]$mutation.copy.number
for(s in 2:length(subsamples[[run]])){
  WTCount = cbind(WTCount,data[[s]]$WT.count)
  mutCount = cbind(mutCount,data[[s]]$mut.count)
  totalCopyNumber = cbind(totalCopyNumber,data[[s]]$subclonal.CN)
  copyNumberAdjustment = cbind(copyNumberAdjustment,data[[s]]$no.chrs.bearing.mut)
  normalCopyNumber = cbind(normalCopyNumber,data[[1]]$normalCopyNumber)
  mutation.copy.number = cbind(mutation.copy.number, data[[1]]$mutation.copy.number)
}

non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment==0)==0)
mutCount = mutCount[non.zero,]
WTCount = WTCount[non.zero,]
totalCopyNumber = totalCopyNumber[non.zero,]
copyNumberAdjustment = copyNumberAdjustment[non.zero,]
normalCopyNumber = normalCopyNumber[non.zero,]
mutation.copy.number = mutation.copy.number[non.zero,]
for(s in 1:length(subsamples[[run]])){
  data[[s]]=data[[s]][non.zero,]
}

return(list(data=data, mutCount=mutCount, WTCount=WTCount, totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, mutation.copy.number=mutation.copy.number, cellularity=cellularity, samplename=samplename, subsamples=subsamples))
}


read_triplet_pipeline_density_input = function(donorname, run, subsample.indices) {
sample=donorname
samplenames=donorname
subsamples=read.table("samplenames.lst")
subsamples[[run]]=as.vector(subsamples$V1)

#get cellularities
ce=read.table("purity.lst")
cellularities = list()
cellularities[[run]]=as.vector(ce$V1)


samplename = samplenames[run]
cellularity = cellularities[[run]]

no.subsamples = length(subsamples[[run]])

data=list()
for(s in 1:no.subsamples){
  data[[s]]=read.table(paste0(subsamples[[run]][s],"_allDirichletProcessInfo.txt"),sep="\t",header=T)
}

new_output_folder = "md_out"
if(!file.exists(new_output_folder)){
  dir.create(new_output_folder)
}

WTCount = data[[1]]$WT.count
mutCount = data[[1]]$mut.count
totalCopyNumber = data[[1]]$subclonal.CN
copyNumberAdjustment = data[[1]]$no.chrs.bearing.mut
normalCopyNumber = data[[1]]$normalCopyNumber
mutation.copy.number = data[[1]]$mutation.copy.number
for(s in 2:length(subsamples[[run]])){
  WTCount = cbind(WTCount,data[[s]]$WT.count)
  mutCount = cbind(mutCount,data[[s]]$mut.count)
  totalCopyNumber = cbind(totalCopyNumber,data[[s]]$subclonal.CN)
  copyNumberAdjustment = cbind(copyNumberAdjustment,data[[s]]$no.chrs.bearing.mut)
  normalCopyNumber = cbind(normalCopyNumber,data[[1]]$normalCopyNumber)
  mutation.copy.number = cbind(mutation.copy.number, data[[1]]$mutation.copy.number)
}

non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment==0)==0)
mutCount = mutCount[non.zero,]
WTCount = WTCount[non.zero,]
totalCopyNumber = totalCopyNumber[non.zero,]
copyNumberAdjustment = copyNumberAdjustment[non.zero,]
normalCopyNumber = normalCopyNumber[non.zero,]
mutation.copy.number = mutation.copy.number[non.zero,]
for(s in 1:length(subsamples[[run]])){
  data[[s]]=data[[s]][non.zero,]
}

mutation.copy.number = mutation.copy.number[,subsample.indices]
copyNumberAdjustment = copyNumberAdjustment[,subsample.indices]
temp.data = list()
for(i in 1:length(subsample.indices)){
  temp.data[[i]] = data[[subsample.indices[i]]]
}
data = temp.data

return(list(data=data, mutCount=mutCount, WTCount=WTCount, totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, mutation.copy.number=mutation.copy.number, cellularity=cellularity, samplename=samplename, subsamples=subsamples))
}

read_triplet_pipeline_input = function() {
  subsamplenames=read.table("samplenames.lst")
  subsamplenames=as.vector(subsamplenames$V1)
  
  #get cellularities
  ce=read.table("purity.lst")
  cellularity=as.vector(ce$V1)
  
  
  #get DP info
  DP.files = paste(subsamplenames,"_allDirichletProcessInfo.txt",sep="")
  print(paste(DP.files))
  
  #naser: exclude unwanted variants using the "non.zero" filter (as opposed to the original "-exclude.indices") as applied in the naser_multidimclustparallel_updated_final.R
  data=list()
  for (s in 1:length(DP.files)){
    data[[s]]=read.table(DP.files[s],sep="\t",header=T)
  }
  #check original dimension
  print(paste(dim(data[[1]])))
  #
  WTCount = data[[1]]$WT.count
  mutCount = data[[1]]$mut.count
  totalCopyNumber = data[[1]]$subclonal.CN
  copyNumberAdjustment = data[[1]]$no.chrs.bearing.mut
  normalCopyNumber = data[[1]]$normalCopyNumber
  mutation.copy.number = data[[1]]$mutation.copy.number
  
  if (length(subsamplenames) > 1) {
    for(s in 2:length(subsamplenames)){
      WTCount = cbind(WTCount,data[[s]]$WT.count)
      mutCount = cbind(mutCount,data[[s]]$mut.count)
      totalCopyNumber = cbind(totalCopyNumber,data[[s]]$subclonal.CN)
      copyNumberAdjustment = cbind(copyNumberAdjustment,data[[s]]$no.chrs.bearing.mut)
      mutation.copy.number = cbind(mutation.copy.number,data[[s]]$mutation.copy.number)
    }
  } else {
    WTCount = as.matrix(as.vector(WTCount))
    mutCount = as.matrix(as.vector(mutCount))
    totalCopyNumber = as.matrix(as.vector(totalCopyNumber))
    copyNumberAdjustment = as.matrix(as.vector(copyNumberAdjustment))
    mutation.copy.number = as.matrix(as.vector(mutation.copy.number))
    
  }
  
  non.zero = which(rowSums(mutCount)>0 & !is.na(rowSums(totalCopyNumber)) & rowSums(copyNumberAdjustment==0)==0)
  mutCount = mutCount[non.zero,,drop=F]
  WTCount = WTCount[non.zero,,drop=F]
  totalCopyNumber = totalCopyNumber[non.zero,,drop=F]
  copyNumberAdjustment = copyNumberAdjustment[non.zero,,drop=F]
  #normalCopyNumber = normalCopyNumber[non.zero,,drop=F]
  mutation.copy.number = mutation.copy.number[non.zero,,drop=F]
  for(s in 1:length(subsamplenames)){
    data[[s]]=data[[s]][non.zero,,drop=F]
  }
  
  
  return(list(data=data, mutCount=mutCount, WTCount=WTCount, totalCopyNumber=totalCopyNumber, copyNumberAdjustment=copyNumberAdjustment, mutation.copy.number=mutation.copy.number, cellularity=cellularity))
}


#' Step 3 of the triplet pipeline : Combine output of the density across triplets into a single consensus
#'
#' @param donorname Identifier to be used to name output files, should be the same as "samplename" in run_dpc_core_triplet
#' @param run Not needed parameter, should really always be set to 1. This code could be run across multiple donors, but has not been adapted to continue to work as such
#' @param noiters The number of total iterations, should be the same as "no.iters" in run_dpc_core_triplet
#' @param burn.in Number of iterations to be used for burnin. This should be the same value as in run_dpc_core_triplet, in which it is defined as ceiling(no.iters/5)
#' @author Naser Ansari-Pour (BDI, Oxford)
#' @export
combine_triplet_runs = function(donorname, run, noiters, burn.in) {

  conc_param=1
  cluster_conc = 5
  density.smooth = 0.01

  dat = read_triplet_pipeline_combine_input(donorname, run)
  attach(dat)

  no.subsamples = length(subsamples[[run]])
  choose.number = 3
  choose.from = no.subsamples
  no.perms = choose(no.subsamples,3)
  new_output_folder = "md_out"

  no.muts = nrow(data[[1]]) #naser: nrow(info) until runDP filter has changed from -exclude.indices to non.zero ??
  print(paste(no.muts))
  node.assignments=NULL
  for(choose.index in 1:no.perms){
    print(choose.index)
    subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
    print(subsample.indices)
    perm.assignments = read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and_cluster_info_0.01.txt",sep=""),header=T,sep="\t")

    if(is.null(node.assignments)){
      node.assignments = array(data.matrix(perm.assignments[ncol(perm.assignments)]),c(no.muts,1))
    }else{
      node.assignments = cbind(node.assignments,data.matrix(perm.assignments[ncol(perm.assignments)]))
    }
  }

  #get different combinations of assignments from the different permutations
  unique.assignments = unique(node.assignments)
  write.table(unique.assignments,paste(new_output_folder,"/",samplename,"_matchedClustersInParallelRuns.txt",sep=""),sep="\t",row.names=F,quote=F)

  no.consensus.nodes = nrow(unique.assignments)
  print(paste("#consensus nodes = ",no.consensus.nodes,sep=""))
  consensus.assignments = vector(mode="numeric",length = no.consensus.nodes)
  for(n in 1:no.consensus.nodes){
    print(n)
    consensus.assignments[sapply(1:no.muts,function(a,u,i){all(u==a[i,])},a=node.assignments,u=unique.assignments[n,])]=n
  }

  #get probabilities of assignment to each set of assignments
  prob.consensus.assignments = array(1,c(no.muts,no.consensus.nodes))
  #prob.consensus.assignments = array(0,c(no.muts,no.consensus.nodes)) #get average probability, rather than product
  for(choose.index in 1:no.perms){
    print(choose.index)
    subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
    perm.assignments = read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_DP_and_cluster_info_0.01.txt",sep=""),header=T,sep="\t")
    for(n in 1:no.consensus.nodes){
      #prob.consensus.assignments[,n] = prob.consensus.assignments[,n] * perm.assignments[,unique.assignments[n,choose.index]+2]
      #prob.consensus.assignments[,n] = prob.consensus.assignments[,n] + perm.assignments[,unique.assignments[n,choose.index]+2]

      prob.consensus.assignments[,n] = prob.consensus.assignments[,n] * perm.assignments[,unique.assignments[n,choose.index]]
    }
  }
  #prob.consensus.assignments = prob.consensus.assignments / no.perms


  out = data[[1]][,1:2]
  for(s in 1:no.subsamples){
    out = cbind(out,data[[s]]$subclonal.fraction)
  }
  out = cbind(out,prob.consensus.assignments)
  names(out) = c("chr","pos",paste(samplename,subsamples[[run]],"_subclonalFraction",sep=""),paste("prob.cluster",1:no.consensus.nodes,sep=""))
  write.table(out,paste(new_output_folder,"/",samplename,"_allClusterProbabilitiesFromParallelRuns_16Oct2014.txt",sep=""),sep="\t",row.names=F,quote=F)

  ##############################################################################################################################################
  all.CIs = array(NA,c(no.consensus.nodes,no.subsamples,no.perms,2))
  for(choose.index in 1:no.perms){
    subsample.indices = get.subsample.indices(choose.index,choose.number,choose.from)
    print(subsample.indices)
    CIs = data.matrix(read.table(paste(new_output_folder,"/",samplename,"_",paste(subsample.indices,collapse="_"),"_confInts_",density.smooth,".txt",sep=""),header=T,sep="\t",stringsAsFactors=F))
    for(n in 1:no.consensus.nodes){
      all.CIs[n,subsample.indices,choose.index,1] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)-1]
      all.CIs[n,subsample.indices,choose.index,2] = CIs[unique.assignments[n,choose.index],2*(1:choose.number)]
    }
  }
  median.CIs = array(NA,c(no.consensus.nodes,no.subsamples,2))
  for(n in 1:no.consensus.nodes){
    for(s in 1:no.subsamples){
      for(c in 1:2){
        median.CIs[n,s,c] = median(all.CIs[n,s,,c],na.rm=T)
      }
    }
  }
  median.CIs.2D = t(sapply(1:no.consensus.nodes,function(m,i){as.vector(t(m[i,,]))},m=median.CIs))
  write.table(cbind(1:no.consensus.nodes,median.CIs.2D,table(consensus.assignments)),sep="\t",row.names=F,quote=F,col.names=c("cluster.no",paste(rep(paste(samplename,subsamples[[run]],sep=""),each=2),rep(c("lowerCI","upperCI"),times=no.subsamples),sep="_"),"no.muts"),paste(new_output_folder,"/",samplename,"_consensusClustersByParallelNodeAssignment_16Oct2014.txt",sep=""))
  out = data[[1]][,1:2]
  for(s in 1:no.subsamples){
    out = cbind(out,data[[s]]$subclonal.fraction)
  }
  out = cbind(out,consensus.assignments)
  names(out) = c("chr","pos",paste(samplename,subsamples[[run]],"_subclonalFraction",sep=""),"cluster.no")
  write.table(out,paste(new_output_folder,"/",samplename,"_allClusterassignmentsFromParallelRuns_23Nov2018.txt",sep=""),sep="\t",row.names=F,quote=F)

  pdf(paste(new_output_folder,"/","consensus_cluster_assignment_",samplename,"_combined_",density.smooth,"_23Nov2018.pdf",sep=""),height=4,width=4)
  #its hard to distinguish more than 8 different colours
  max.cols = 8
  cols = rainbow(min(max.cols,no.consensus.nodes))
  plot.data = mutation.copy.number/copyNumberAdjustment
  plot.data[is.na(plot.data)]=0

  for(i in 1:(no.subsamples-1)){
    for(j in (i+1):no.subsamples){
      plot(plot.data[,i],plot.data[,j],type = "n",xlab = paste(samplename,(subsamples[[run]][subsample.indices])[i]," subclonal fraction",sep=""), ylab = paste(samplename,subsamples[[run]][j]," subclonal fraction",sep=""),xlim = c(0,max(plot.data[,i])*1.25))
      for(n in 1:no.consensus.nodes){
        pch=20 + floor((n-1)/max.cols)
        #pch is not implmeneted above 25
        if(pch>25){
          pch=pch-20
        }
        points(plot.data[,i][consensus.assignments==n],plot.data[,j][consensus.assignments==n],col=cols[(n-1) %% max.cols + 1],pch=pch,cex=0.5)
      }
      pch=20 + floor((0:(no.consensus.nodes-1))/max.cols)
      pch[pch>25] = pch[pch>25]-20
      legend(max(plot.data[,i])*1.05,max(plot.data[,j]),legend = 1:no.consensus.nodes,col=cols[(0:(no.consensus.nodes-1)) %% max.cols + 1],pch=pch,cex=1)
    }
  }
  dev.off()
}
