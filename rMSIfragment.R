##################################################
## LIPID FRAGMENTATION IDENTIFICATION WORKFLOW  ##
## Gerard Baquer Gómez                          ##
## 29/09/20                                     ##
##################################################

# 0. LITERATURE DATA
lipids<-data.frame(name=c("DG","TG","PA","PE","PC","PG","PI","PS","CL","Cer","CerP","SM","HexCer","SFT","GM3","FRAG"),stringsAsFactors = F)
lipids$standard <- c("DG 18:1","TG 18:1","PA 18:1", "PE 16:0/18:1 & PE P-18:1","PC 16:0/18:1 & PC P-18:1","PG 18:0","PI 16:0/18:1","PS 16:0","CL 16:1","Cer d18:/18:0","CerP d18:1/16:0","SM d18:1/12:0","GlcCer d18:1/18:0","SFT d18:1/12:0","Ganglioside GM3")
lipids$synonyms <-c("Dag|Diacylglycerol","Tag|Triacylglycerol","Phosphatidic acid|Phosphatidate",
                    "Phophatidylethanolamine|GPEtn","Phosphatidylcholine|GPCho","Phosphatidylglycerol|GPG",
                    "Phosphatidylinositol|Pino","Phosphatidylserine|Pser","Cardiolipin","DHCer|Ceramide",
                    "Ceramide phosphate","Sphingomyelin","","","","")
lipid_dict<-as.list(lipids$name)
names(lipid_dict)<-lipids$name
lipid_dict$PEtOH <- "PE"
lipid_dict$DGTS <- "DG"
lipid_dict$SHexCer <- "HexCer"
lipid_dict$exCer<-"HexCer"

adducts_pos<-data.frame(name=c("M+Na-COOH","M-H2O+H","M+H","M-H2O+Na","M-H20+K","M+Na","M+K","M-H+2Na","M-H+Na+K","M-2H+3Na","M-H+2K","M-2H+2Na+K","M-2H+Na+2K","M-2H+3K","M+"),stringsAsFactors = F)
adducts_pos$mass<-c(-21.0006,-14.9876,1.0073,6.9943,22.9682,22.9892,38.9632,44.9712,60.9451,66.9531,76.9190,82.9271,98.9010,114.8749,0)
adducts_pos$mode<-rep("pos",nrow(adducts_pos))
adducts_pos$DG<-c(0,1,0,1,1,3,1,0,0,0,0,0,0,0,1)
adducts_pos$TG<-c(0,0,0,0,0,3,1,0,0,0,0,0,0,0,1)
adducts_pos$PA<-c(0,0,1,0,0,1,1,3,1,1,1,1,1,1,1)
adducts_pos$PE<-c(0,0,1,0,0,2,1,3,1,0,1,0,0,0,1)
adducts_pos$PC<-c(0,0,2,0,0,3,1,0,0,0,0,0,0,0,1)
adducts_pos$PG<-c(0,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$PI<-c(0,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$PS<-c(0,0,0,0,0,1,1,1,1,3,1,1,1,1,1)
adducts_pos$CL<-c(0,0,0,0,0,1,1,2,1,3,1,1,1,1,1)
adducts_pos$Cer<-c(0,2,1,1,1,3,1,0,0,0,0,0,0,0,1)
adducts_pos$CerP<-c(0,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$SM<-c(0,0,2,0,0,3,1,0,0,0,0,0,0,0,1)
adducts_pos$HexCer<-c(0,0,0,0,0,3,2,0,0,0,0,0,0,0,1)
adducts_pos$SFT<-c(0,0,0,0,0,0,0,3,0,0,1,0,0,0,1)
adducts_pos$GM3<-c(1,0,0,0,0,1,1,3,1,0,1,0,0,0,1)
adducts_pos$FRAG<-c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)

adducts_neg<-data.frame(name=c("M-2(H2O)","M-oh","M-Ch3","M-H","M+Na-2h","M+K-2h","M+Dan-H","M-"),stringsAsFactors = F)
adducts_neg$mass<-c(-36.0212,-17.0027,-15.0229,-1.0073,20.9747,36.9486,157.0771,0)
adducts_neg$mode<-rep("neg",nrow(adducts_neg))
adducts_neg$DG<-c(0,0,0,0,0,0,0,1)
adducts_neg$TG<-c(0,0,0,0,0,0,0,1)
adducts_neg$PA<-c(0,0,0,3,0,0,1,1)
adducts_neg$PE<-c(0,0,0,3,0,0,0,1)
adducts_neg$PC<-c(0,0,3,0,0,0,0,1)
adducts_neg$PG<-c(0,0,0,3,0,0,0,1)
adducts_neg$PI<-c(0,0,0,3,0,0,0,1)
adducts_neg$PS<-c(0,0,0,3,1,1,0,1)
adducts_neg$CL<-c(0,0,0,1,3,1,0,1)
adducts_neg$Cer<-c(0,0,0,3,0,0,0,1)
adducts_neg$CerP<-c(0,0,0,3,0,0,1,1)
adducts_neg$SM<-c(0,0,3,0,0,0,0,1)
adducts_neg$HexCer<-c(0,0,0,3,0,0,0,1)
adducts_neg$SFT<-c(1,1,0,3,0,0,0,1)
adducts_neg$GM3<-c(0,0,0,3,0,0,0,1)
adducts_neg$FRAG<-c(0,0,0,0,0,0,0,1)
  
adducts<-rbind(adducts_pos,adducts_neg)
rownames(adducts)<-adducts$name

fragmentation<-data.frame(matrix("",ncol=nrow(lipids),nrow=nrow(lipids)),stringsAsFactors = F)
colnames(fragmentation)<-lipids$name
rownames(fragmentation)<-lipids$name

fragmentation_pos<-fragmentation
fragmentation_neg<-fragmentation
rm(fragmentation)

fragmentation_pos['PC','PA']<-"-N(CH3)3-Cho"
fragmentation_pos['PE','PA']<-"-EtnA"
fragmentation_pos['PI','PA']<-"-Inositol"
fragmentation_pos['PS','PA']<-"-Serine"

fragmentation_pos['PC','DG']<-"-P_Cho"
fragmentation_pos['PE','DG']<-"-P_EtnA"
fragmentation_pos['PI','DG']<-"-P_Inositol"

fragmentation_pos['SM','CerP']<-"-N(CH3)3-Cho"
fragmentation_pos['SM','Cer']<-"-P_Cho"

fragmentation_pos['GM3','FRAG']<-"-Sialic_Acid"

fragmentation_pos['CL','FRAG']<-"-DG+H2O | -DG | -RCOOH"

fragmentation_pos[,'FRAG']<-paste(fragmentation_pos[,'FRAG'],"| -NH2 | -NH3")

fragmentation_neg['PC','PA']<-"-CN(CH3)3 | -N(CH3)3 | -Cho"
fragmentation_neg['PS','PA']<-"-Serine"

fragmentation_neg['PI','FRAG']<-"-C4H10O5"

fragmentation_neg['SM','CerP']<-"-CN(CH3)3' | -N(CH3)3 | -Cho"

fragmentation_neg['GM3','Cer']<-"-GM3HG"
fragmentation_neg['GM3','HexCer']<-"-Sialic_Acid-Hex"
fragmentation_neg['GM3','FRAG']<-"-Sialic_Acid"
fragmentation_neg['GM3','Cer']<-"-Sialic_Acid-Lac"

fragmentation_neg['HexCer','Cer']<-"-Hex"

fragmentation_neg['SFT','FRAG']<-"-OH | 2H2O"

fragmentation_neg['CL','FRAG']<-"-RCOOH-RCOOH | -RCOOH-RCOH | -RCOOH | -RCOH | -DG+H2O | -DG"

fragmentation_neg[,'FRAG']<-paste(fragmentation_neg[,'FRAG'],"| -NH2 | -NH3")

losses<-NULL

losses['P_Cho']<-183.0660444
losses['P_EtnA']<-141.0190942
losses['P_Inositol']<-259.0218935
losses['N(CH3)3']<-59.0734993
losses['Cho']<-85.08914936
losses['EtnA']<-43.04219917
losses['Inositol']<-162.0528234
losses['Serine']<-87.03257699
losses['Sialic_Acid']<-291.0954165
losses['CN(CH3)3']<-71.0734993
losses['Hex']<-162.0528234
losses['Lac']<-324.1056468
losses['OH']<-17.00273965
losses['2H2O']<-36.02112937
losses['C4H10O5']<-138.0528234
losses['GM3HG']<-615.2010634
losses['NH2']<- 16.01872407 
losses['NH3']<- 17.02654911

# 1. HELPER FUNCTIONS
parse_fragmentation<- function(s)
{
  s<- gsub(" ", "", s, fixed = T)
  n<-unlist(strsplit(s,"|",fixed=T))
  v<-sapply(strsplit(n,"-",fixed = T),function(x)(-1)*sum(losses[x],na.rm = T))
  names(v)<-n
  return(v)
}

find_c_delta <- function(n)
{
  pattern<-"([0-9]+:[a-z]?[0-9]+)"
  matches<-gsub("[a-z]","",unlist(regmatches(n, gregexpr(pattern, n))))
  if(length(matches)==0)
    result=c(NA,NA)
  else
    result<-apply(matrix(as.numeric(unlist(strsplit(matches,":"))),ncol=length(matches)),1,sum)
  return(result)
}
find_lipid <- function(n,s)
{
  pattern<-paste("(",s,")+",sep="") #pattern<-paste("[",s,"]\\w+",sep="")
  return(unlist(regmatches(n, gregexpr(pattern, n)))[1])
}
find_lipid_c_delta <- function(n,s=paste(lipids$name,collapse ="|")){
  return(paste(find_lipid(n,s),paste(find_c_delta(n),collapse = ":")))
}
find_adduct <- function(n)
{
    pattern<-"([+|-][^0-9][^]]+)"
    return(gsub(" ","",paste("M",unlist(regmatches(n, gregexpr(pattern, n)))[1],sep="")))
}
get_tol<-function(m1,m2){
  return(abs((m1-m2)/m2)*10^6)
}
withintol<-function(m1,m2,tol){
  return(get_tol(m1,m2)<tol)
}
load_gt<-function(path="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_neg.csv"){
  table<-read.csv(path,sep=";")
  #lipid
  table$lipid<-sapply(table$assignement,function(s)find_lipid(s,paste(names(adducts)[-(1:3)],collapse="|")))
  #C_delta
  c_delta<-sapply(table$assignement,find_c_delta)
  table$c<-c_delta[1,]
  table$delta<-c_delta[2,]
  #adduct
  table$adduct<-sapply(table$assignement,find_adduct)
  
  table<-subset(table,!is.na(lipid)&!grepl("?",adduct,fixed = T)&!grepl(":",adduct,fixed = T))
  return(table)
}
#This should be parsing, load lipidmaps should only load the already stored db
load_lipidMAPS <- function(path="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/db/LMSD.sdf"){
  #Load sdf file
  lipidMAPS<-ChemmineR::read.SDFset(path)
  metadata<-ChemmineR::datablock(lipidMAPS)
  table<-data.frame(t(sapply(metadata,function(x)x[c("LM_ID","NAME","CATEGORY","EXACT_MASS","FORMULA","ABBREVIATION","INCHI_KEY","INCHI" )])),stringsAsFactors = F)
  #Remove overhead
  rm(lipidMAPS)
  rm(metadata)
  
  #Refactor names to match the package's standards
  names(table)<-tolower(names(table))
  names(table)[names(table)=="exact_mass"]<-"exactmass"
  names(table)[names(table)=="lm_id"]<-"id"
  table$exactmass<-as.numeric(table$exactmass)
  
  #parse lipid, c and delta
  table[c("c","delta")]<-t(sapply(table$abbreviation,find_c_delta))
  table$lipid<-sapply(table$abbreviation,function(x)find_lipid(x,paste(lipids$name,collapse ="|")))
  
  return(table)
}
load_RefMet<-function(path="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/refmet.csv"){
  RefMet<-read.csv(path)
  #Trim RefMet
  lipids_pattern<-paste(lipids$name,collapse ="|")
  RefMet_lipids<-subset(RefMet,grepl(lipids_pattern,refmet_name))
  RefMet_lipids[c("c","delta")]<-t(sapply(RefMet_lipids$refmet_name,find_c_delta))
  RefMet_lipids<-subset(RefMet_lipids,!is.na(c))
  RefMet_lipids$lipid<-sapply(RefMet_lipids$refmet_name,function(x)find_lipid(x,lipids_pattern))
  RefMet_lipids$change<-"None"
  
  #Temporary fix
  RefMet_lipids_min_oxigen<-RefMet_lipids
  RefMet_lipids_min_oxigen$exactmass<-RefMet_lipids$exactmass-15.99491462
  RefMet_lipids_min_oxigen$change<-"O"
  RefMet_lipids<-rbind(RefMet_lipids,RefMet_lipids_min_oxigen)
  
  RefMet_lipids$exactmass<-round(RefMet_lipids$exactmass,4)
  RefMet_lipids<-RefMet_lipids[rownames(unique(RefMet_lipids[c("lipid","c","delta","change","exactmass")])),]
  RefMet_lipids<-RefMet_lipids[order(RefMet_lipids$exactmass),]
  
  
  #Do not filter (or maybe increase the subclasses available)
  #sc<-lipids$name
  #sc[which(!lipids$name%in%unique(RefMet$sub_class))]<-c("DAG","TAG","Cer-1-P", "", "Gangliosides", "")
  #RefMet_lipids<-subset(RefMet_lipids,sub_class%in%sc)
  return(RefMet_lipids)
}
get_tree<-function(parent,relation){
  
  offspring=unlist(subset(relation,TYPE=="is_a"&INIT_ID%in%parent[[1]],FINAL_ID))
  if((length(offspring))!=0){
    names(offspring)<-paste("LVL",length(parent),"_ID",1:length(offspring),sep="")
    parent=c(list(offspring),parent)
    return(get_tree(parent,relation))
  }
  else
    return(unlist(parent,relation))
  
}
load_decoy_non_animal<-function(r=NA,path="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/chebi/"){
  ## Load ChEBI
  relation<-read.table(file = paste(path,"relation.tsv",sep=""), sep = '\t', header = TRUE)
  chemical_data<-read.table(file = paste(path,"chemical_data_3star.tsv",sep=""), sep = '\t', header = TRUE)
  #n<-read.table(file = paste(path,"names_3star.tsv",sep=""), sep = '\t', header = TRUE,fill=T)
  
  ids=c(76924,84735,76969,76946,35703,78298)
  names(ids)<-c("plant metabolite","algal metabolite","bacterial metabolite","fungal metabolite","xenobiotic","environmental contaminant")
  id_trees<-lapply(ids,function(i)get_tree(i,relation))
  non_animal_compounds<-unlist(subset(relation,TYPE=="has_role"&INIT_ID%in%unlist(id_trees),FINAL_ID))
  decoy_db<-subset(chemical_data,TYPE=="MONOISOTOPIC MASS"&COMPOUND_ID%in%non_animal_compounds,c(COMPOUND_ID,CHEMICAL_DATA))
  formulas<-subset(chemical_data,TYPE=="FORMULA",c(COMPOUND_ID,CHEMICAL_DATA))
  decoy_db$formula<-formulas$CHEMICAL_DATA[match(decoy_db$COMPOUND_ID,formulas$COMPOUND_ID)]
  #decoy_db$abbreviation<-n$NAME[match(decoy_db$COMPOUND_ID,n$COMPOUND_ID)]
  decoy_db$abbreviation<-decoy_db$COMPOUND_ID
  colnames(decoy_db)<-c("id","exactmass","formula","abbreviation")
  decoy_db$exactmass<-as.numeric(as.character(decoy_db$exactmass))
  
  
  if(!is.na(r))
    decoy_db<-subset(decoy_db,exactmass>r[1]&exactmass<r[2])
  
  #Expand with synthetic_processes
  synthetic_processes<-c(0,15.99491462,-15.99491462, 30.01056468, 162.0528234, -18.01056468, 15.0234751,1.007825032)
  names(synthetic_processes)<-c("","O","-O","CH2O","C6H10O5","-H2O","CH3","H")#Added methylation and protonation
  
  big_decoy_db<-do.call(rbind,replicate(length(synthetic_processes),decoy_db,simplify = FALSE))
  big_decoy_db$exactmass<-big_decoy_db$exactmass+rep(synthetic_processes,each=nrow(decoy_db))
  big_decoy_db$synth_proc<-rep(names(synthetic_processes),each=nrow(decoy_db))
  
  big_decoy_db$lipid<-with(big_decoy_db,paste(id,synth_proc,sep="_"))
  big_decoy_db$c<-0
  big_decoy_db$delta<-0
  
  return(big_decoy_db)
}
get_elements_decoy<-function(elems_target = c("H","D","C","O","N","Na","K","Cl")){
  data("isotopes",package = "enviPat",envir = environment())
  isotopes_o<-isotopes[order(isotopes$abundance,decreasing = T),]
  elems<-with(isotopes_o,sapply(unique(element),function(x)mass[which(element==x)[1]]))
  elems<-sort(elems)
  #remove isotopes
  elems<-elems[which(!(grepl("]",names(elems))))]
  #remove targets
  elems<-elems[which(!(names(elems)%in%elems_target))]
}
get_decoy_adduct<-function(elems,mode="pos"){
  return(sample(elems,1))
}
get_decoy_fragment<-function(elems,mode="pos"){
  tmp<-sample(elems,sample(1:2,1))
  frag<--sum(tmp)
  names(frag)<-paste(c("",names(tmp)),sep="",collapse="-")
  return(frag)
}
get_lipids <- function(s){
  return(unique(sapply(strsplit(as.character(s),",")[[1]],find_lipid_c_delta)))
}
load_mat_pks <- function(path="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/mat/150514_N29_neg25um/No_Norm__Masked__Calibrated_Analysis"){
  
  p<-list()
  p$mass<-c(R.matlab::readMat(paste(path,"/mz.mat",sep=""))$mz)
  p$intensity<-t(R.matlab::readMat(paste(path,"/Matriz.mat",sep=""))$Matriz)
  attr(p$intensity,"Csingle")<-NULL
  p$pos<-cbind(c(R.matlab::readMat(paste(path,"/xy.mat",sep=""))$x),c(R.matlab::readMat(paste(path,"/xy.mat",sep=""))$y))
  colnames(p$pos)<-c("x","y")
  p$posMotors<-p$pos
  p$numPixels<-nrow(p$intensity)
  p$SNR<-array(5,dim(p$intensity))
  p$area<-array(1,dim(p$intensity))
  fs<-strsplit(path,"/")[[1]]
  p$names<-fs[length(fs)-1]
  class(p)<-"rMSIprocPeakMatrix"
  return(p)
}

merge_mat_pks <- function(path="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/mat",mode="pos"){
  paths<-paste(dir(path,full.names = T,pattern = mode),"/No_Norm__Masked__Calibrated_Analysis",sep="")
  pks_list<-lapply(paths,load_mat_pks)
  return(rMSIproc::MergePeakMatrices(pks_list))
}
store_mat_pks<- function(){
  rMSIproc::StorePeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/pks_neg.zip",merge_mat_pks(mode="neg"))
  rMSIproc::StorePeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/pks_pos.zip",merge_mat_pks(mode="pos"))
}

load_db_from_scratch <- function(){
  db_t<-load_lipidMAPS()
  db_t<-db_t[!is.na(db_t$lipid),]
  db_t<-db_t[with(db_t,match(unique(abbreviation),abbreviation)),]
  db_d<-load_decoy_non_animal(c(min(db_t$exactmass),max(db_t$exactmass)))
}
adjust_db_densities <- function(db_1,db_2){
  #sample
  t_density<-density(db_1$exactmass)
  while(nrow(db_2)>=2*nrow(db_1)){
    d_density<-density(db_2$exactmass)
    db_2<-db_2[sample(1:nrow(db_2),nrow(db_2)/2,prob=approx(t_density$x,t_density$y,db_2$exactmass)$y/approx(d_density$x,d_density$y,db_2$exactmass)$y),]
  }
  return(db_2[sample(1:nrow(db_2),nrow(db_1)),])
}


store_db <- function(){
  save(db_t,db_d,file="/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/target_decoy_db.RData")
}
load_db <-function(){
  load("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/target_decoy_db.RData")
}


#Custom class functions
class(res_r_d1)<-c("rMSILipidAnnot","data.frame")
print.rMSILipidAnnot<-function(x){
  title=as.character(substitute(x))
  #class(x)<-"data.frame"
  if(is.null(x$correlation))
    View((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","adduct_score")]),title)
  else
    View((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","lipid_occurences","correlation","adduct_score")]),title)
  
}

simply<-function(x){
  if(is.null(x$correlation))
    return((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","adduct_score")]))
  else
    return((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","lipid_occurences","correlation","adduct_score")]))
  
}
#2. Main functions

search_db<-function(mz,db,tol=5,m="pos",rand_af=F,db_ref=db){
  #Generate a and f
  db_mz<-db$exactmass
  a<-subset(adducts,mode==m)$mass
  names(a)<-subset(adducts,mode==m)$name
  if(m=="neg")
    frag<-fragmentation_neg
  else
    frag<-fragmentation_pos
  f<-unlist(lapply(unique(unlist(frag)),parse_fragmentation))
  f<-f[!(f==0)]
  f<-f[unique(names(f))]
  f<-c(f,0)
  
  if(rand_af){
    elems<-get_elements_decoy()
    offset<-matrix(sapply(1:(length(a)*length(f)),function(x)get_decoy_adduct(elems)+get_decoy_fragment(elems)),nrow = length(a),ncol=length(f))
  }
  else{
    offset<-outer(a,f,function(x,y)x+y)
  }
  
  #Generate target
  target<-outer(mz,offset,function(x,y)x-y)
  
  #Generate matches
  target_min<-target-target*tol*10^-6
  target_max<-target+target*tol*10^-6
  
  matches<-lapply(db_mz,function(x)which((target_min<x)&(target_max>x),arr.ind = T))
  
  #Generate final data.table
  results<-do.call(rbind,lapply(seq_along(matches),function(i,x=matches[[i]])if(nrow(x)==0) NULL else data.frame(db[i,c("lipid","c","delta","exactmass","formula","id","abbreviation")],mz_i=x[,1],adduct_i=x[,2],fragment_i=x[,3],lipid_ref=db_ref[i,"lipid"])))
  results$experimental_mz<-mz[results$mz_i]
  results$adduct<-names(a)[results$adduct_i]
  results$adduct_mode<-subset(adducts,mode==m)$mode[results$adduct_i]
  results$adduct_score<-sapply(1:nrow(results),function(i)subset(adducts,mode==m)[results$adduct_i[i],as.character(results$lipid_ref[i])])
  results$adduct_score<-sapply(results$adduct_score,function(x)if(is.null(x))0 else x)
  results$fragmentation<-names(f)[results$fragment_i]
  results$fragmentation_possible<-sapply(1:nrow(results),function(i)frag_possible(results$adduct_mode[i],as.character(results$lipid_ref[i]),results$fragmentation[i]))
  results$ppm_error=get_tol(results$exactmass,results$experimental_mz-a[results$adduct_i]-f[results$fragment_i])
  r=regexpr("(?<=C)[0-9]*", results$formula, perl=T)
  results$carbon_number=rep(0,nrow(results))
  results$carbon_number[r!=-1]=as.numeric(regmatches(results$formula,r))
  results$carbon_parity=((results$carbon_number%%2)==0)
  class(results)<-c("rMSILipidAnnot","data.frame")
  rownames(results)<-NULL
  return(unique(results)) #Temporary fix. Find where the entries are being doubled.
}
frag_possible<-function(m,l,frag){
  return(frag %in% possible_frags(m,l))
}

possible_frags<-function(m,l){
  if(m=="neg")
    f=fragmentation_neg
  else
    f=fragmentation_pos
  s<- paste(f[l,],collapse="|")
  s<- gsub(" ", "", s, fixed = T)
  n<-unique(unlist(strsplit(s,"|",fixed=T)))
  return(n)
}
cor_score<-function(m){
  if(sum(m!=1)==0)
    return(0)
  else
    return(mean(m[m!=1]))
}
ranking<-function(res,pks,m){
  res<-unique(res)
  #res<-res[!is.na(res$lipid),]
  res$lipid_string<-paste(res$lipid," ",res$c,":",res$delta, " ",round(res$exactmass,2),sep="")
  res$lipid_id<-as.numeric(as.factor(res$lipid_string))
  res$lipid_occurences<-sapply(res$lipid_id,function(x)sum(res$lipid_id==x))
  res$lipid_parents<-sapply(res$lipid_id,function(x)sum(res$lipid_id==x&res$fragmentation==""))
  res$lipid_children<-sapply(res$lipid_id,function(x)sum(res$lipid_id==x&res$fragmentation!=""))
  res$ambiguity<-sapply(res$mz_i,function(x)sum(res$mz_i==x))
  res$dominance<-sapply(1:nrow(res),function(i)res$lipid_occurences[i]/sum(res$lipid_occurences[res$mz_i==res$mz_i[i]]))
  res$parental_percentage<-rep(0,nrow(res))
  res$fragment_percentage<-rep(0,nrow(res))
  for(i in unique(res$lipid_id)){
    p_adduct<-unique(subset(res,lipid_id==i&fragmentation=="")$adduct)
    possible_adducts<-rownames(subset(adducts,mode==m))[which(subset(adducts,mode==m)[,res$lipid_ref[match(i,res$lipid_id)]]!=0)]
    res$parental_percentage[res$lipid_id==i]<-sum(possible_adducts%in%p_adduct)/length(possible_adducts)
    
    f_frag<-unique(subset(res,lipid_id==i&fragmentation!="")$fragmentation)
    possible_fragmentations<-possible_frags(m,res$lipid_ref[match(i,res$lipid_id)])
    possible_fragmentations<-possible_fragmentations[possible_fragmentations!=""]
    res$fragment_percentage[res$lipid_id==i]<-sum(possible_fragmentations%in%f_frag)/length(possible_fragmentations)
  }
  
  c<-cor(pks$intensity)
  
  res$correlation<-0
  res$km_correlation<-0
  res$km_correlation_old<-0
  res$km_correlation_class<-"Not computed"
  
  for(i in unique(res$lipid_id))
  {
    is<-which(res$lipid_id==i)
    initial_mz_i<-unique(res$mz_i[is])
    mz_i=initial_mz_i
    max_cor=cor_score(c[mz_i,mz_i])
    min_cor=max_cor
    res$correlation[is]<-max_cor
    class=rep(1,length(mz_i))
    
    #Current approach
    if(length(mz_i) > 2){
      km<-kmeans(c[mz_i,mz_i],2)
      a=mz_i[km$cluster==1]
      b=mz_i[km$cluster==2]
      res$km_correlation[is[initial_mz_i%in%a]]=cor_score(c[a,a])
      res$km_correlation[is[initial_mz_i%in%b]]=cor_score(c[b,b])
    }
    else{
      res$km_correlation[is]=cor_score(c[mz_i,mz_i])
    }
    #Original approach
    while(max_cor<0.7 && length(mz_i) > 2)
    {
      km<-kmeans(c[mz_i,mz_i],2)
      a=mz_i[km$cluster==1]
      b=mz_i[km$cluster==2]
      a_cor=cor_score(c[a,a])
      b_cor=cor_score(c[b,b])
      max_cor=max(a_cor,b_cor)
      min_cor=min(a_cor,b_cor)
      if(a_cor>b_cor)
        mz_i=a
      else
        mz_i=b
      
    }
    
    res$km_correlation_old[is]=max_cor
    res$km_correlation_class[is]="Excluded"
    if(max_cor>=0.7){
      res$km_correlation_old[is[!initial_mz_i%in%mz_i]]=0
      res$km_correlation_class[is[!initial_mz_i%in%mz_i]]="Excluded"
      res$km_correlation_class[is[initial_mz_i%in%mz_i]]="Selected"
    }
    
  }
  
  res$lipid_occurences_selected<-sapply(res$lipid_id,function(x)sum(res$lipid_id==x&res$km_correlation_class=="Selected"))
  res$S<-with(res,lipid_occurences*(1+correlation))
  return(res)
}
annotate<-function(pks,db,tol=5,m="pos",rand_af=F,db_ref=db){
  res<-search_db(pks$mass,db,tol,m,rand_af,db_ref)
  res_rank<-ranking(subset(res,adduct_mode==m&fragmentation_possible),pks,m)
  return(res_rank)
}

#A. Manual identification validation

#B. Ranking evaluation
target_decoy_validation<-function(pks,db_t,db_d,mode="pos",tol=5,d_type=T){
  
  #Target
  mz<-pks$mass
  t<-annotate(pks,db_t,tol,mode)
  
  if(d_type){
    #Decoy 1: Non-Animal compunds
    d<-annotate(pks,db_d,tol,mode,db_ref = db_t)
  }
  else{
    #Decoy 2: Random adducts and fragments
    d<-annotate(pks,db_t,tol,mode,rand_af = T)
  }
  
  #Final results
  td<-rbind(t,d)
  td$db<-c(rep("T",nrow(t)),rep("D",nrow(d)))
  
  return(td)
}

#ROC curve
td_roc<-function(scores,weights,mz){
  r<-list()
  r$scores<-sort(unique(scores))
  
  r$curve<-data.frame(fpr=rep(0,length(r$scores)))
  r$curve$fpr<-sapply(r$scores,function(x)sum(!weights[scores>=x])/sum(!weights))
  r$curve$tpr<-sapply(r$scores,function(x)sum(weights[scores>=x])/sum(weights))
  
  r$curve2<-data.frame(hits=rep(0,length(r$scores)))
  r$curve2$hits<-sapply(r$scores,function(x)sum(scores>=x)/length(unique(mz)))
  r$curve2$fdr<-sapply(r$scores,function(x)sum(!weights[scores>=x])/sum(weights[scores>=x]))
  r$curve2<-r$curve2[r$curve2$hits>0.05,]
  i<-order(r$curve$fpr)
  r$auc.integral<-sum(diff(r$curve$fpr[i])*(head(r$curve$tpr[i],-1)+tail(r$curve$tpr[i],-1)))/2
  return(r)
}
plot_roc <- function(res,metrics,names=metrics,mode="ROC"){
  if( mode=="ROC")
    roc<-lapply(metrics,function(x)with(res,td_roc(eval(parse(text=x)),db=="T",mz_i)))
  else
    roc<-lapply(metrics,function(x)with(res,PRROC::pr.curve(scores.class0 = eval(parse(text=x)),weights.class0 = db=="T",curve=T)))
  
  if(mode=="ROC"){
    plot(0,type="l",col=1,xlab="Hits per mz",ylab="False Discovery Rate (Decoy/Target)",xlim=c(0,18),ylim=c(0,1))
    for (i in 1:length(roc)){
      lines(roc[[i]]$curve2,col=i,lw=2)
    }
    legend(8.5,1.05,paste(names,round(sapply(roc,function(x)x$auc.integral),2)),col=1:(length(roc)),lw=2)
  }
  
  plot(0:1,0:1,type="l",col=1,asp=1,ylab="True Positive Rate (Target Rate)",xlab="False Positive Rate (Decoy Rate)")
  for (i in 1:length(roc)){
      lines(roc[[i]]$curve,col=i,lw=2)
  }
  legend(0.5,0.4,paste(names,round(sapply(roc,function(x)x$auc.integral),2)),col=1:(length(roc)),lw=2)
  
  
}

fig_dir<-"/home/gbaquer/msidata/1. In-source Fragmentation/1.6. Paper/Figures/"
plot_figure1<-function(res,metrics,names=metrics,label=""){
  #Target Decoy Validation
  roc<-lapply(metrics,function(x)with(res,td_roc(eval(parse(text=x)),db=="T",mz_i)))
  pr<-lapply(metrics,function(x)with(res,PRROC::pr.curve(scores.class0 = eval(parse(text=x)),weights.class0 = db=="T",curve=T)))

  df<-data.frame(x=unlist(sapply(roc,function(x)x$curve$fpr)),y=unlist(sapply(roc,function(x)x$curve$tpr)),var=rep(names,times=sapply(roc,function(x)nrow(x$curve))))
  
  pa<-ggplot(df,aes(x,y,col=var))+geom_line(size=2)+
    xlab("False Positive Rate (Decoy Rate)")+ylab("True Positive Rate (Target Rate)")+
    theme_bw()+theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("A")+
    scale_color_discrete(name = "Score", breaks = names, labels = paste(names," (",round(sapply(roc,function(x)x$auc.integral),2)," AUC)",sep = ""))
  
  df<-data.frame(x=unlist(sapply(pr,function(x)x$curve[,1])),y=unlist(sapply(pr,function(x)x$curve[,2])),var=rep(names,times=sapply(pr,function(x)nrow(x$curve))))
  
  pb<-ggplot(df,aes(x,y,col=var))+geom_line(size=2)+
    xlab("Recall")+ylab("Precision")+
    theme_bw()+theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("B")+
    scale_color_discrete(name = "Score", breaks = names,labels = paste(names," (",round(sapply(pr,function(x)x$auc.integral),2)," AUC)",sep = ""))
  
  
  
  df<-data.frame(x=unlist(sapply(roc,function(x)x$curve2$hits)),y=unlist(sapply(roc,function(x)x$curve2$fdr)),var=rep(names,times=sapply(roc,function(x)nrow(x$curve2))))
  
  pc<-ggplot(df,aes(x,y,col=var))+geom_line(size=2)+
    xlab("Annotations per MS signal")+ylab("False Discovery Rate (Decoy/Target)")+
    theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("C")+
    scale_color_discrete(name = "Score", breaks = names,labels = paste(names," (",round(sapply(roc,function(x)x$auc.integral),2)," AUC)",sep = ""))
  
  p=grid.arrange(pa,pc,nrow=1)
  ggsave(paste("Fig1_Target_Decoy_Validation",".tiff",sep=label),p,path=fig_dir,width=20,height=10)
}

plot_figure2<-function(r,t,s="S",steps=1000,label=""){
  #Manual validation
  df<-data.frame(o=rep(0,steps),o_hplc=rep(0,steps),th=seq(min(r[,s]),max(r[,s]),length.out=steps))
  df$n<-sapply(df$th,function(x)sum(r[,s]>=x))
  for(i in seq_along(df$th)){
    m<-manual_validation(r[r[,s]>=df$th[i],],t)
    df$o[i]<-sum(m$table$match)/nrow(m$table)
    df$o_hplc[i]<-sum(subset(m$table,HPLC_adduct!="_")$match)/nrow(subset(m$table,HPLC_adduct!="_"))
  }
  c<-max(df$n)
  dfp<-data.frame(threshold=rep(df$th,3),y=with(df,c(o,o_hplc,n/c)),var=rep(c("Matches (%)","Matches(%) (HPLC validated)","Annotations per MS signal"),each=steps))
  pa<-ggplot(dfp,aes(threshold,y,col=var))+geom_line(size=2)+
    scale_y_continuous(
    name = "Matches to Garate et al. 2020 (%)",
    sec.axis = sec_axis(~.*c/length(pks$mass), name="Annotations per MS signal")
  ) +xlab("Ranking Score (S) threshold")+theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("A")
  
  m<-manual_validation(r,t)
  dfp<-data.frame(score=m$res[,s],type=as.factor(match(m$res$match_HPLC_adduct,c("","_"),nomatch=3)))
  pb<-ggplot(dfp,aes(type,score,col=type))+geom_violin()+geom_point()+
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.5)+
    xlab("")+ylab("Ranking Score (S)")+theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                          panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("B")
  
  p=grid.arrange(pa,pb)
  ggsave(paste("Fig2_Manual_Validation",".tiff",sep=label),p,path=fig_dir,width=10,height=10)
}

require(data.table)
plot_figure2_topN<-function(r,t,s="S",steps=5,label=""){
  #Manual validation
  df<-data.frame(o=rep(0,steps),o_hplc=rep(0,steps),th=steps:1)
  r<-manual_validation(r,t)$res
  for(i in seq_along(df$th)){
    m<-manual_validation(topN(r,s,df$th[i]),t)
    df$o[i]<-sum(m$table$match)/nrow(m$table)
    df$o_hplc[i]<-sum(subset(m$table,HPLC_adduct!="_")$match)/nrow(subset(m$table,HPLC_adduct!="_"))
  }
  dfp<-data.frame(threshold=rep(df$th,2),y=with(df,c(o,o_hplc)),var=rep(c("Matches (%)","Matches(%) (HPLC validated)"),each=steps))
  pa<-ggplot(dfp,aes(threshold,y,col=var))+geom_line(size=2)+geom_point(size=10)+xlab("Top N for each mz")+ylab("Matches to Garate et al. 2020 (%)")+theme_bw()+theme(legend.justification = c(1, 0), legend.position = c(1, 0),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                            panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("A")+ylim(0,0.75)
  
  m<-manual_validation(r,t)
  dfp<-data.frame(score=m$res[,s],type=as.factor(match(m$res$match_HPLC_adduct,c("","_"),nomatch=3)))
  pb<-ggplot(dfp,aes(type,score,col=type))+geom_violin()+geom_point()+
    stat_summary(fun = "median",
                 geom = "crossbar", 
                 width = 0.5)+
    xlab("")+ylab("Ranking Score (S)")+theme_bw()+theme(legend.justification = c(1, 1), legend.position = c(1, 1),legend.background=element_rect(size=0.2,colour=1),panel.border = element_rect(colour = "black", fill=NA), panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0),text=element_text(size=20)) +
    ggtitle("B")
  
  p=pa
  ggsave(paste("Fig2_Manual_Validation_TopN",".tiff",sep=label),p,path=fig_dir,width=10,height=10)
}

plot_figure3<-function(pks,res,mz,topN=NA,score="lipid_occurences*(1+correlation)",label=""){
  res<-res[order(with(res,eval(parse(text=score))),decreasing=T),]
  i<-which.min(abs(mz-pks$mass))
  x<-subset(res,mz_i==i)
  if(!is.na(topN)&nrow(x)>0){
    x<-x[1:min(nrow(x),topN),]
  }
  k<-list()
  for(i in 1:nrow(x)){
    y<-subset(res,lipid_id==x$lipid_id[i])
    z<-matrix(0,nrow(y)+1,nrow(y)+1)
    n=with(rbind(x[i,],y),paste(round(experimental_mz,2),"\n[",adduct," ",fragmentation,"]",sep=""))
    n[1]<- paste(x$abbreviation[i],x$formula[i],n[1],paste("LO:",x$lipid_occurences[i],"C:",round(x$correlation[i],2)),sep="\n")
    z[1,-1]=1
    k[[i]]<-list(z=z,n=n,t=c(0,1+as.numeric(y$fragmentation=="")))
  }
  net = network(do.call(magic::adiag,sapply(k,function(x)x$z)), directed = T)
  pa<-ggnet2(net,mode="kamadakawai",label=unlist(sapply(k,function(x)x$n)),shape=15,label.size=3.5,size=20,color=1+unlist(sapply(k,function(x)x$t)))+
    ggtitle("A")+
    scale_color_brewer(palette="Set3",labels = c("Annotation for current mz","Fragment","Parental ion"))
  
  
  
  y<-rbind(x[1,],subset(res,lipid_id==x$lipid_id[1]))
  pb<-grid.arrange(top=grid::textGrob("B", x = 0, hjust = 0),grobs=lapply(1:nrow(y),function(i)ggplot_peak_image(one_pks(pks,4),y$mz_i[i],k[[1]]$n[i],k[[1]]$t[i]+1)))
  
  p<-grid.arrange(pa,pb,ncol=1)
  ggsave(paste("Fig3_Example_Annotations",".tiff",sep=label),p,path=fig_dir,width=10,height=20)
}

load_metaspace<-results
update_table_formulas<-function(t,db)
{
  t$formula<-db_t$formula[match(with(t,paste(lipid," ",c,":",delta,sep="")),db$abbreviation)]
  return(t)
}
metaspace_validation<-function(res,metaspace_res_path,t){
  res_meta<-read.csv(metaspace_res_path)
  t<-update_table_formulas(t,db_t)
  res<-manual_validation(res,t)$res
  t<-manual_validation(res,t)$table
  res_meta$table_mz<-sapply(res_meta$mz,function(x)with(t,mz[which.min(abs(mz-x))]))
  a<-with(t,paste(formula,mz))
  b<-with(res_meta,paste(formula,table_mz))
  c<-with(res,paste(formula,table_mz))
  
  res$match_f_t<-c%in%a
  res$match_f_m<-c%in%b
  
  res_meta$match_f_t<-b%in%a
  res_meta$match_f_e<-b%in%c
  
  t$match_f_e<-a%in%c
  t$match_f_m<-a%in%b
  
  return(list(res=res,res_meta=res_meta,table=table))
}
one_pks<-function(pks,i){
  j<-cbind(1+c(0,cumsum(pks$numPixels)[-length(pks$numPixels)]),cumsum(pks$numPixels))
  return(pks[j[i,1]:j[i,2],])
}
ggplot_mean_spectra<- function(pks,i1=NA,i2=NA,i3=NA){
  col=rep("black",length(pks$mass))
  
  p=ggplot(data.frame(mz=pks$mass,i=apply(pks$intensity,2,mean)),
           aes(x=mz, ymax=i, ymin=0,col=col)) +
    geom_linerange()
  
  if(!is.na(i1))
    p<-p+ geom_vline(xintercept = pks$mass[as.numeric(i1)], linetype="dotted", 
                      color = "red", size=1.5)
  if(!is.na(i2))
    p<-p+ geom_vline(xintercept = pks$mass[as.numeric(i2)], linetype="dotted", 
                      color = "green", size=1.5)
  if(!is.na(i3))
    p<-p+ geom_vline(xintercept = pks$mass[as.numeric(i3)], linetype="dotted", 
                      color = "blue", size=1.5)
  return(p)
}
ggplot_peak_image <- function(pks,i,title="",col=1)
{
  palette=c("#8CD3C7","#FDFFB2","#BFBAD9")
  df=data.frame(x=pks$pos[,2],y=pks$pos[,1],z=pks$intensity[,i])
  background=element_rect(fill = palette[col])
  p=ggplot(df, aes(x, y, fill = z)) + geom_raster() +
    coord_fixed(1,expand = F) +
    ggtitle(title) +
    scale_fill_gradientn(colours=viridis(1000))+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
           axis.text.y=element_blank(),axis.ticks=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),legend.position="none",
           panel.background=element_rect(fill = palette[col]),panel.border=element_blank(),panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),plot.background=background,plot.title = element_text(hjust = 0.5,size=9))
  return(p)
}
upsample <- function(res){
  x<-res
  new_rows<-nrow(x)+(1:(with(res,sum(db=="T")-sum(db=="D"))))
  x[new_rows,]<-NA
  x[new_rows,"db"]<-"D"
  x[new_rows,c("lipid_occurences","correlation","km_correlation","km_correlation_old","fragment_percentage","parental_percentage")]<-0
  x[new_rows,"lipid_occurences"]<-1
  return(x)
}
downsample <- function(res){
  return(res[c(sample(which(res$db=="T"),sum(res$db=="D")),which(res$db=="D")),])
}
#Plot curve
tmp<-outer(pks$mass,pks$mass,function(x,y)x-y)
tmp2<-table(round(tmp,2))
plot(tmp2,xlim=c(0.01,500))

##Figures
#Figure 2:
#Normal DB
library(ggnet)
library(network)
library(sna)
library(ggplot2)
library(ggplot2)
library(gridExtra)
pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/pks_neg.zip")
metrics<-c("lipid_occurences","correlation","lipid_occurences*(correlation+1)","parental_percentage","fragment_percentage","km_correlation","km_correlation_old")

res1<-target_decoy_validation(pks,db_t,db_d,mode="neg",tol=5)
res2<-target_decoy_validation(pks,db_t,db_d,mode="neg",tol=5,d_type = F)

plot_figure1(res1,metrics,label="_NEG_D1_N")
plot_figure1(upsample(res1),metrics,label="_NEG_D1_U")
plot_figure1(downsample(res1),metrics,label="_NEG_D1_D")

plot_figure1(res2,metrics,label="_NEG_D2_N")
plot_figure1(upsample(res2),metrics,label="_NEG_D2_U")
plot_figure1(downsample(res2),metrics,label="_NEG_D2_D")

res<-annotate(pks,db_t,20,"neg")
res5ppm<-annotate(pks,db_t,5,"neg")
plot_figure3(pks,res5ppm,744.558,3,label="_5ppm_top3")

plot_figure2(res5ppm,load_gt(),"adduct_score",label="_adduct_score")
plot_figure2_topN(res,load_gt())
#Downsample T
plot_roc(res1[c(sample(which(res$db=="T"),sum(res$db=="D")),which(res$db=="D")),],c("lipid_occurences","correlation","lipid_occurences*correlation","lipid_occurences*(correlation+1)","parental_percentage","fragment_percentage","km_correlation","km_correlation_old"))
plot_roc(res2[c(sample(which(res$db=="T"),sum(res$db=="D")),which(res$db=="D")),],c("lipid_occurences","correlation","lipid_occurences*correlation","lipid_occurences*(correlation+1)","parental_percentage","fragment_percentage","km_correlation","km_correlation_old"))

#Upsample D
plot_roc(tmp,c("lipid_occurences","correlation","lipid_occurences*correlation","lipid_occurences*(correlation+1)*(fragment_percentage+1)","parental_percentage","fragment_percentage","km_correlation","km_correlation_old"))

tmp2<-res2
new_rows2<-nrow(tmp2)+(1:(with(res2,sum(db=="T")-sum(db=="D"))))
tmp[new_rows2,]<-NA
tmp2[new_rows2,"db"]<-"D"
tmp2[new_rows2,c("lipid_occurences","correlation","km_correlation","km_correlation_old","fragment_percentage","parental_percentage")]<-0
tmp2[new_rows2,"lipid_occurences"]<-1
plot_roc(tmp2,c("lipid_occurences","correlation","lipid_occurences*correlation","lipid_occurences*(correlation+1)*(fragment_percentage+1)","parental_percentage","fragment_percentage","km_correlation","km_correlation_old"))

#Smaller DB
table_neg<-load_gt()
is<-with(db_t,paste(lipid,c,delta))%in%unique(with(table_neg,paste(lipid,c,delta)))
res<-target_decoy_validation(pks,db_t[is,],db_d[is,],mode="neg",tol=5)
plot_roc(res,c("lipid_occurences","correlation","lipid_occurences*correlation"))


man_val<-manual_validation(subset(res1,db=="T"),load_gt())
manual_validation<-function(res,table){
  
  res$table_mz<-sapply(res$experimental_mz,function(x)with(table,mz[which.min(abs(mz-x))]))
  a<-with(table,paste(lipid,c,delta,mz))
  b<-with(res,paste(lipid,c,delta,table_mz))
  res$match_HPLC_adduct<-sapply(match(b,a),function(i)if(is.na(i))""else as.character(table$HPLC_adduct[i]))
  res$match<-b%in%a
  table$match<-a%in%b
  
  return(list(res=res,table=table))
}
plot_manual_validation<-function(x){
  res<-x$res
  table<-x$table
  plot(density(subset(res,adduct_mode=="neg"&fragmentation_possible&match)$lipid_occurences))
  lines(density(subset(res,adduct_mode=="neg"&fragmentation_possible&!match)$lipid_occurences),col=2)
  plot(density(subset(res,adduct_mode=="neg"&fragmentation_possible&match)$correlation))
  lines(density(subset(res,adduct_mode=="neg"&fragmentation_possible&!match)$correlation),col=2)
  plot(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&match),ambiguity*lipid_occurences)))
  lines(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&!match),ambiguity*lipid_occurences)),col=2)
  plot(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&match),correlation*lipid_occurences)))
  lines(density(with(subset(res,adduct_mode=="neg"&fragmentation_possible&!match),correlation*lipid_occurences)),col=2)
}
plot_manual_percentage <- function(res,scores,steps=100)
{
  table=load_gt()
  
  for(s in scores)
  {
    o<-rep(0,steps)
    th<-seq(0,max(res[,s]),length.out=steps)
    for(i in seq_along(th)){
      m<-manual_validation(res[res[,s]>th[i],],table)
      o[i]<-sum(m$table$match)/nrow(m$table)
    }
    plot(th,o)
  }
  return(o)
}

plot_manual_percentage_2 <- function(res,scores,steps=100)
{
  table=load_gt()
  
  for(s in scores)
  {
    o<-rep(0,steps)
    th<-seq(0,nrow(res),length.out=steps)
    a<-order(res[,s])
    for(i in seq_along(th)){
      m<-manual_validation(res[a[1:th[i]],],table)
      o[i]<-sum(m$table$match)/nrow(m$table)
    }
    plot(th,o)
  }
  return(o)
}

topN<-function(d,s="S",n=3)
{
  return(do.call("rbind",lapply(unique(d[["mz_i"]]),function(i)head(subset(d[order(d[[s]],decreasing = T),],mz_i==i),n))))
}
#Comparison against selected METASPACE datasets

#####################
## SHINY DASHBOARD ##
#####################
#tagList
library(shiny)
library(shinydashboard)
library(DT)
is<-(1:length(pks$mass))
names(is)<-as.character(round(pks$mass,2))
ui <- dashboardPage(
  dashboardHeader(title = "rMSIfragment explorer"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("MSI Dataset", tabName = "dataset", icon = icon("shapes")),
      menuItem("Annotations", tabName = "annotations", icon = icon("th-list")),
      menuItem("Explorer", tabName = "explorer", icon = icon("project-diagram"))
    )
  ),
  dashboardBody(
    tabItems(
      # Sample Content
      tabItem(tabName = "dataset",
              fluidRow(plotOutput("meanSpectra")),
              fluidRow(
                column(4,selectInput("i1","Ion 1",is),plotOutput("ion1")),
                column(4,selectInput("i2","Ion 2",is),plotOutput("ion2")),
                column(4,selectInput("i3","Ion 3",is),plotOutput("ion3"))
              )
      ),
      #Annotation Content
      tabItem(tabName = "annotations",
              fluidRow(
                column(12,div(DTOutput("res"), style = "overflow-y: auto;"))
              ),
              downloadButton("downloadData", "Download .csv")
      ),
      
      # Explorer Content
      tabItem(tabName = "explorer",
              fluidRow(column(2,selectInput("mz","Choose mz",is)),column(2,selectInput("n","Page",1:5)),column(3,radioButtons("nID","Choose network",choices=1:5,inline=T))),
              fluidRow(div(style = 'overflow-y: scroll',plotOutput("networks"))),
                fluidRow(plotOutput("plots"))
      )
    )
  )
)

server <- function(input, output) {
  
  output$res <- renderDT({datatable(simply(res),filter = "top")%>% 
      formatRound(columns = c(1,2,8), digits = 2)})
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("annotation-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(res, file)
    }
  )
  output$ion1 <- renderPlot({
    ggplot_peak_image(one_pks(pks,4),as.numeric(input$i1))
    
  })
  output$ion2 <- renderPlot({
    ggplot_peak_image(one_pks(pks,4),as.numeric(input$i2))
  })
  output$ion3 <- renderPlot({
    ggplot_peak_image(one_pks(pks,4),as.numeric(input$i3))
  })
  
  output$meanSpectra <-renderPlot({
    ggplot_mean_spectra(one_pks(pks,4),input$i1,input$i2,input$i3)
  })
  
  k<-reactive({plot_networks(one_pks(pks,4),res5ppm,pks$mass[as.numeric(input$mz)])})
  
  output$networks <- renderPlot({
    display_networks(k(),as.numeric(input$n),5)
  })
  
  output$plots <- renderPlot({
    plot_ions(one_pks(pks,4),res5ppm,pks$mass[as.numeric(input$mz)],as.numeric(input$nID)+(as.numeric(input$n)-1)*5)
  })
}
plot_networks<-function(pks,res,mz,topN=NA,score="lipid_occurences*(1+correlation)"){
  res<-res[order(with(res,eval(parse(text=score))),decreasing=T),]
  i<-which.min(abs(mz-pks$mass))
  x<-subset(res,mz_i==i)
  if(!is.na(topN)&nrow(x)>0){
    x<-x[1:min(nrow(x),topN),]
  }
  k<-list()
  for(i in 1:nrow(x)){
    y<-subset(res,lipid_id==x$lipid_id[i])
    z<-matrix(0,nrow(y)+1,nrow(y)+1)
    n=with(rbind(x[i,],y),paste(round(experimental_mz,2),"\n[",adduct," ",fragmentation,"]",sep=""))
    n[1]<- paste(x$abbreviation[i],x$formula[i],n[1],paste("LO:",x$lipid_occurences[i],"C:",round(x$correlation[i],2)),sep="\n")
    z[1,]=1
    net = network(z, directed = T)
    k[[i]]<-ggnet2(net,mode="kamadakawai",label=n,shape=15,label.size=3.5,size=20,color=1+c(0,as.numeric(y$fragmentation=="")))+
      scale_color_brewer(palette="Set3",labels = c("Annotation for current mz","Fragment","Parental ion"))+theme(aspect.ratio = 1)
  }
  return(k)
}
display_networks<-function(k,i,n=5){
  p<-grid.arrange(grobs=k[1:length(k)%in%(((i-1)*5+1):(i*n))],nrow=1)
  return(p)
}
plot_ions<-function(pks,res,mz,j,topN=NA,score="lipid_occurences*(1+correlation)"){
  
  res<-res[order(with(res,eval(parse(text=score))),decreasing=T),]
  i<-which.min(abs(mz-pks$mass))
  x<-subset(res,mz_i==i)
  if(!is.na(topN)&nrow(x)>0){
    x<-x[1:min(nrow(x),topN),]
  }
  y<-rbind(x[j,],subset(res,lipid_id==x$lipid_id[j]))
  p<-grid.arrange(grobs=lapply(1:nrow(y),function(i)ggplot_peak_image(pks,y$mz_i[i])))
  return(p)
}
shinyApp(ui, server)

##

###############
## SANDBOX 1 ##
###############

#Validate results

gt <- table

results$matched_lipid<-paste(results$lipid,results$c,results$delta,results$experimental_mz)%in%paste(gt$lipid,gt$c,gt$delta,gt$mz)
r<-subset(results,adduct_mode=="neg")
gt$matched_lipid<-paste(gt$lipid,gt$c,gt$delta,gt$mz)%in%paste(r$lipid,r$c,r$delta,r$experimental_mz)
gt$in_RefMet<-paste(gt$lipid,gt$c,gt$delta)%in%paste(RefMet_lipids$lipid,RefMet_lipids$c,RefMet_lipids$delta)
subset(gt,!matched_lipid)
#gt$matched_lipid_MS2ID<-table$matched_lipid


#Cross Validation with METASPACE
results_pos<-read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_pos.csv")
results_neg<-read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_neg.csv")

metaspace_pos_LM <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_pos_LM.csv",sep=";")
metaspace_pos_HMDB <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_pos_HMDB.csv",sep=";")
metaspace_pos_ChEBI<- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_pos_ChEBI.csv",sep=";")
metaspace_neg_LM <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_neg_LM.csv",sep=";")
metaspace_neg_HMDB <- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_neg_HMDB.csv",sep=";")
metaspace_neg_ChEBI<- read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/metaspace_annotations_neg_ChEBI.csv",sep=";")

metaspace_pos_LM$moleculeNames[1:10]

metaspace_pos_LM_lipids<-lapply(metaspace_pos_LM$moleculeNames,get_lipids)
metaspace_neg_LM_lipids<-lapply(metaspace_neg_LM$moleculeNames,get_lipids)

metaspace_pos_LM$paper_mz<-sapply(metaspace_pos_LM$mz,function(x)results_pos$mz[which.min(abs(x-results_pos$mz))])
metaspace_neg_LM$paper_mz<-sapply(metaspace_neg_LM$mz,function(x)results_neg$mz[which.min(abs(x-results_neg$mz))])

results_pos$matched_METASPACE<-paste(results_pos$lipid,paste(results_pos$c,results_pos$delta,sep=":"),results_pos$mz)%in%unlist(lapply(metaspace_pos_LM_lipids,function(x)paste(x,metaspace_pos_LM$paper_mz)))
results_neg$matched_METASPACE<-paste(results_neg$lipid,paste(results_neg$c,results_neg$delta,sep=":"),results_neg$mz)%in%unlist(lapply(metaspace_neg_LM_lipids,function(x)paste(x,metaspace_neg_LM$paper_mz)))


write.csv(results_pos,"/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_pos_metaspace.csv")
write.csv(results_neg,"/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results_neg_metaspace.csv")

# 3. DISPLAY FUNCTIONS
bd<-"/home/gbaquer/msidata/1. In-source Fragmentation/1.3. AuBSi Wells/20200701_AuBSi_wells_std/20200701_AuBSi_wells_std/"
pks<-rMSIproc::LoadPeakMatrix(paste(bd,"20200701_AuBSi_wells_std-peaks.zip",sep=""))
pks_list<-lapply(1:3,function(i)rMSIcleanup:::get_one_peakMatrix(pks,i))


#Edit position

pks_neg[[1]]$pos[1:1800,1]<-rep(1:90,each=20)
pks_neg[[1]]$pos[1:1800,2]<-rep(1:20,each=90)
pks_neg[[1]]$pos[1801,]<-c(91,1)
pks_neg[[1]]$pos[1802,]<-c(91,2)


#Check results (method pere, lluc, maria)

fn<-subset(gt,!matched_lipid)
unique(fn$adduct)[-c(1,2,3,4,14)]

names(offset[1,]) #names of fragmentations
names(offset[,1]) #names of adducts

#M-N(CH3)3(CH2)-H is missing

table(fn$adduct)
table(gt$adduct)
#Most of the missing adducts are M-CH3 and M-H (focus on this ones)
# M-NH2+DAN-H was not found at all [SOLVED]



subset(fn,adduct=="M-H")


Refsubset(gt,adduct=="M-H"&matched_lipid)


 #Results metaspace
pks_METASPACE<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/3. METASPACE/Brain01_animal1_s1235.zip")
results_METASPACE<-search_db(mz=pks_METASPACE$mass)
results_METASPACE_filtered<-subset(results_METASPACE,fragmentation_possible==T&adduct_score!="0"&adduct_mode=="pos"&adduct!="M+")
results_METASPACE_ranked<-ranking(results_METASPACE_filtered,pks_METASPACE)
nrow(subset(results_METASPACE_ranked,correlation_class=="Selected"&lipid_occurences_selected>1))

plot(cumsum(table(table(results_METASPACE_filtered$experimental_mz)))/671*100,ylab = "Cumulative ion percentage",xlab="Number of annotations / ion",ylim=c(0,100),type="l")
lines(cumsum(table(table(subset(results_METASPACE_ranked,correlation_class=="Selected"&lipid_occurences_selected>1)$experimental_mz)))/355*100,ylab = "Cumulative ion percentage",xlab="Number of annotations / ion",ylim=c(0,100),col="red")
lines(cumsum(table(table(subset(results_METASPACE_ranked,correlation_class=="Selected"&lipid_occurences_selected>1&adduct_score==3)$experimental_mz)))/155*100,ylab = "Cumulative ion percentage",xlab="Number of annotations / ion",ylim=c(0,100),col="blue")
a<-nrow(results_METASPACE_filtered)/length(pks_METASPACE$mass)
b<-nrow(subset(results_METASPACE_ranked,correlation_class=="Selected"&lipid_occurences_selected>1))/length(pks_METASPACE$mass)
c<-nrow(subset(results_METASPACE_ranked,correlation_class=="Selected"&lipid_occurences_selected>1&adduct_score==3))/length(pks_METASPACE$mass)
legend(10,30,paste(c("No filter","Cor>0.8 & Lipid>1","Cor>0.8 & Lipid>1 & Max adduct score"),"(",round(c(a,b,c),1),"mean ann/ion)"),c("black","red","blue"))


table_pos<-load_gt("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_pos.csv")
results_Garete_pos<-search_db(mz=unique(table_pos$mz))
results_Garete_pos_filtered<-subset(results_Garete_pos,fragmentation_possible==T&adduct_score!="0"&adduct_mode=="pos"&adduct!="M-")

#pos
pks_pos<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/imzml/161108_38002_N95_pos_25um-peaks.zip")
results_pos<-search_db(mz=pks_pos$mass)
results_pos_filtered<-subset(results_pos,fragmentation_possible==T&adduct_score!="0"&adduct_mode=="pos"&adduct!="M+")
results_pos_ranked<-ranking(results_pos_filtered,pks_pos)
nrow(subset(results_pos_ranked,correlation_class=="Selected"&lipid_occurences_selected>1))

plot(cumsum(table(table(results_pos_filtered$experimental_mz)))/1258*100,ylab = "Cumulative ion percentage",xlab="Number of annotations / ion",ylim=c(0,100),type="l")
lines(cumsum(table(table(subset(results_pos_ranked,correlation_class=="Selected"&lipid_occurences_selected>1)$experimental_mz)))/567*100,ylab = "Cumulative ion percentage",xlab="Number of annotations / ion",ylim=c(0,100),col="red")
lines(cumsum(table(table(subset(results_pos_ranked,correlation_class=="Selected"&lipid_occurences_selected>1&adduct_score==3)$experimental_mz)))/204*100,ylab = "Cumulative ion percentage",xlab="Number of annotations / ion",ylim=c(0,100),col="blue")
a<-nrow(results_pos_filtered)/length(pks_pos$mass)
b<-nrow(subset(results_pos_ranked,correlation_class=="Selected"&lipid_occurences_selected>1))/length(pks_pos$mass)
c<-nrow(subset(results_pos_ranked,correlation_class=="Selected"&lipid_occurences_selected>1&adduct_score==3))/length(pks_pos$mass)
legend(10,30,paste(c("No filter","Cor>0.8 & Lipid>1","Cor>0.8 & Lipid>1 & Max adduct score"),"(",round(c(a,b,c),1),"mean ann/ion)"),c("black","red","blue"))

gt<-load_gt("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_pos.csv")
results<-results_pos_ranked
results$matched_lipid<-paste(results$lipid,results$c,results$delta,results$experimental_mz)%in%paste(gt$lipid,gt$c,gt$delta,gt$mz)
gt$experimental_mz<-unlist(sapply(gt$mz,function(x)if(!is.na(x))pks_pos$mass[which.min(abs(pks_pos$mass-x))] else 0))
r<-results_pos_ranked
r$gt_mz<-unlist(sapply(r$experimental_mz,function(x)if(!is.na(x))unique(gt$mz)[which.min(abs(unique(gt$mz)-x))] else 0))

r<-results_pos_ranked
gt$matched_lipid_filter_0<-paste(gt$lipid,gt$c,gt$delta)%in%paste(r$lipid,r$c,r$delta)
r<-subset(results_pos_ranked,correlation_class=="Selected"&lipid_occurences_selected>1)
gt$matched_lipid_filter_1<-paste(gt$lipid,gt$c,gt$delta)%in%paste(r$lipid,r$c,r$delta)
r<-subset(results_pos_ranked,correlation_class=="Selected"&lipid_occurences_selected>1&adduct_score==3)
gt$matched_lipid_filter_2<-paste(gt$lipid,gt$c,gt$delta)%in%paste(r$lipid,r$c,r$delta)


sum(gt$matched_lipid_filter_0)/nrow(gt)*100
sum(gt$matched_lipid_filter_1)/nrow(gt)*100
sum(gt$matched_lipid_filter_2)/nrow(gt)*100

gt$matched_lipid_2<-sapply(1:nrow(gt),function(i)any(present(gt[i,],r,tol=30)))
gt$in_RefMet<-paste(gt$lipid,gt$c,gt$delta)%in%paste(RefMet_lipids$lipid,RefMet_lipids$c,RefMet_lipids$delta)
subset(gt,!matched_lipid)

gts<-rep(F,nrow(gt))
for (i in 1:nrow(gt)){
  rs<-rep(F,nrow(r))
  for(j in 1:nrow(r)){
    rs[j]<-present(gt[i,],r[j,],15)
  }
  gts[i]<-any(rs)
}

present<-function(x,y,tol=15){
  return(x$lipid==y$lipid & x$c==y$c & x$delta==y$delta & withintol(x$mz,y$experimental_mz,tol))
}
#Ranking workflow
pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/imzml/151030_30794_N95_Maldi1_25um_neg-peaks.zip")
res<-search_db(mz=pks$mass)

nrow(subset(res,correlation_class=="Selected"&lipid_occurences_selected>2))



#A Validation Ranking HPLC
#A1. Run annotation
pks_pos<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/imzml/161108_38002_N95_pos_25um-peaks.zip")
results_pos<-search_db(mz=pks_pos$mass)

#A2. Match against gt
table_pos<-load_gt("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_pos.csv")

results_pos$gt_i<-sapply(results_pos$experimental_mz,function(x)if(min(get_tol(x,table_pos$mz),na.rm = T)<5) which.min(get_tol(x,table_pos$mz)) else 0)

results_pos$matched_id<-paste(round(results_pos$experimental_mz,2),results_pos$lipid,results_pos$c,results_pos$delta)%in%paste(round(table_pos$mz,2),table_pos$lipid,table$c,table_pos$delta)


results_pos_filtered<-subset(results_pos,fragmentation_possible==T&adduct_score!="0"&adduct_mode=="pos"&adduct!="M+")
results_pos_filtered_gtmz<-subset(results_pos_filtered,gt_i!=0)


#A3. Run ranking
results_pos_ranked<-ranking(results_pos_filtered,pks_pos)

#A4. Ranking performance metrics plots
HPLC_i<-which(table_pos$HPLC_adduct!="-"& table_pos$HPLC_adduct!="")


###############
## SANDBOX 2 ##
###############

#Garete tables for calibration

pos<-load_gt("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_pos.csv")
neg<-load_gt("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_neg.csv")

p<-cbind(1:length(unique(pos$mz)),unique(pos$mz))
n<-cbind(1:length(unique(neg$mz)),unique(neg$mz))

write.table(p, "/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_pos.ref", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)
write.table(n, "/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/gt_neg.ref", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE)

################
## DEPRECATED ##
################

RefMet<-read.csv("/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/refmet.csv")
load("/home/gbaquer/msidata/sp_MS2ID_20191002_125914.RData")
pks<-rMSIproc::LoadPeakMatrix("/home/gbaquer/msidata/1. In-source Fragmentation/1.2. Standard Study /Munster_Au_Res_140k_10um.zip")

files_neg<-list.files(path = "/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/imzml",pattern=".zip",full.names = T)
pks_neg<-lapply(files_neg,rMSIproc::LoadPeakMatrix)

results_neg<-lapply(pks_neg,function(x)identify_fragments(x,20,"neg",F,F))
T1<-do.call(rbind,lapply(results_neg,function(x)cbind(x$T1,data.frame(dataset=x$parameters$names[1]))))
T2<-do.call(rbind,lapply(results_neg,function(x)cbind(x$T2,data.frame(dataset=x$parameters$names[1]))))
T1$complete_id<-paste(T1$id,T1$dataset,sep="_")
T2$complete_parental_id<-paste(T2$parental_id,T2$dataset,sep="_")
T2$complete_fragment_id<-paste(T2$fragment_id,T2$dataset,sep="_")
T2_filters<-filter_fragments(T2,mode="neg")

MS2ID_lipids<-do.call(rbind,lapply(1:nrow(lipids),function(i)generate_lipid_list(lipids$name[i],lipids$synonyms[i])))
#T1$dataset<-rep(sapply(pks_neg,function(x)x$names[1]),each=sapply(results_neg,function(x)nrow(x$T1)))
#T2$dataset<-rep(sapply(pks_neg,function(x)x$names[1]),each=sapply(results_neg,function(x)nrow(x$T2)))

searchMS2ID<-function(mz,m_adduct,tol=5,substring=""){
  matches<-subset(subset(df_metametabolits,withintol(monoisotopic_molecular_weight,mz-m_adduct,tol)),grepl(substring,nommetabolit))
  unique_result=NULL
  if(nrow(matches)>0){
    c_delta<-matrix(unlist(lapply(matches$nommetabolit,find_c_delta)),nrow=2)
    
    result<-data.frame(metabolite_id=matches$idmetabolit,mz=matches$monoisotopic_molecular_weight,stringsAsFactors = F)
    result$original_lipid<-unlist(lapply(matches$nommetabolit,function(n)find_lipid(n,substring)))
    result$lipid<-unlist(lapply(result$original_lipid,function(x)if(is.null(unlist(lipid_dict[x])))"FRAG" else unlist(lipid_dict[x])))
    result$c<-c_delta[1,]
    result$delta<-c_delta[2,]
    result$ppm_error<-get_tol(result$mz,(mz-m_adduct))
    
    unique_result<-result[!duplicated(result[,c('lipid','c','delta')]),]
  }
  return(unique_result)
}

# 2. IDENTIFICATION FUNCTIONS
identify_fragments<-function(pks,tol=5, m="pos",filter_monoisotopic=T,perform_T2=T){
  #Get Monoisotopic
  if(filter_monoisotopic){
    monoisotopic_i<-(rMSIproc:::isotopeAnnotation(pks,tolerance = tol))$monoisotopicPeaks
  }else{
    monoisotopic_i<-1:length(pks$mass)
  }
  
  a<-subset(adducts,mode==m)
  if(m=="pos"){
    fragmentation<-fragmentation_pos
  }else{
    fragmentation<-fragmentation_neg
  }
  
  #Perform all DB searches [It can be optimized]
  
  T1_db_matches<- data.frame(NULL,stringsAsFactors = F)
  for(i in monoisotopic_i){
    for(j in 1:nrow(a)){
      #before: (l in colnames(a)[which(a[j,-(1:3)]!=0)+3])
      for(l in colnames(a)[-(1:3)]){
        fs=c(0,parse_fragmentation(fragmentation[l,]))
        for(k in seq_along(fs))
        {
          r<-searchMS2ID(pks$mass[i],a$mass[j]+as.numeric(fs[k]),tol,substring=l)
          if(!is.null(r)){
            r$i<-i
            r$adduct<-a$name[j]
            r$experimental_mz<-pks$mass[i]
            r$fragmentation<-names(fs)[k]
          }
          T1_db_matches<-rbind(T1_db_matches,r)
        }
      }
    }
  }
  
  # T1_db_matches<- data.frame(NULL,stringsAsFactors = F)
  # for(i in monoisotopic_i){
  #   for(j in 1:nrow(a)){
  #     r<-searchMS2ID(pks$mass[i],a$mass[j],tol,substring=paste(colnames(a)[which(a[j,-(1:3)]!=0)+3],collapse = "|"))
  #     if(!is.null(r)){
  #       r$i<-i
  #       r$adduct<-a$name[j]
  #       r$experimental_mz<-pks$mass[i]
  #     }
  #     T1_db_matches<-rbind(T1_db_matches,r)
  #   }
  # }
  
  T1_db_matches<-T1_db_matches[order(-T1_db_matches$mz),]
  T1_db_matches$id<-1:nrow(T1_db_matches)
  
  #Update adduct likelyhood and times found
  T1_db_matches$adduct_likelyhood<-unlist(lapply(1:nrow(T1_db_matches),function(i){a[T1_db_matches[i,]$adduct,T1_db_matches[i,]$lipid]}))
  T1_db_matches$times_found<-0
  times_found<-aggregate(list(numdup=rep(1,nrow(T1_db_matches))), T1_db_matches[c('lipid','c','delta')], length)
  for(i in 1:nrow(times_found)){
    tf<-times_found[i,]
    T1_db_matches[rownames(subset(T1_db_matches,lipid==tf$lipid&c==tf$c&delta==tf$delta)),]$times_found<-tf$numdup
  }
  
  #Perform all correlation checks
  T2_possible_fragments<- data.frame(NULL,stringsAsFactors = F)
  if(perform_T2)
  {
    c<-cor(pks$intensity[monoisotopic_i,monoisotopic_i])
    
    
    for(l in lipids$name){
      parentals<-subset(T1_db_matches,lipid==l)
      fragment_lipids<-names(fragmentation)[which(fragmentation[l,]!="")]
      fragments<-subset(T1_db_matches,lipid%in%fragment_lipids)
      if(nrow(parentals)!=0&nrow(fragments)!=0){
        for(i in 1:nrow(parentals)){
          parental=parentals[i,]
          for(j in 1:nrow(fragments)){
            fragment=fragments[j,]
            if(parental$mz>fragment$mz){
              f<-data.frame(parental_id=parental$id,parental_experimental_mz=parental$experimental_mz,parental_lipid=parental$lipid,parental_c=parental$c,parental_delta=parental$delta,
                            fragment_id=fragment$id,fragment_experimental_mz=fragment$experimental_mz,fragment_lipid=fragment$lipid,fragment_c=fragment$c,fragment_delta=fragment$delta,
                            cor=c[match(parental$i,monoisotopic_i),match(fragment$i,monoisotopic_i)],diff=parental$mz-fragment$mz, combined_adduct_likelyhood=NA,
                            combined_ppm_error=parental$ppm_error+fragment$ppm_error,adduct_match=(parental$adduct==fragment$adduct),stringsAsFactors = F)
              T2_possible_fragments<-rbind(T2_possible_fragments,f)
            }
          }
        }
      }
      
    }
    #Update putative fragments and putative parentals
    
    T1_db_matches$putative_fragments<-unlist(lapply(T1_db_matches$id,function(x)sum(T2_possible_fragments$parental_id==x)))
    T1_db_matches$putative_parents<-unlist(lapply(T1_db_matches$id,function(x)sum(T2_possible_fragments$fragment_id==x)))
    
    
  }
  
  
  
  
  
  return(list(T1=T1_db_matches,T2=T2_possible_fragments, parameters=list(names=pks$names,tol=tol,mode=mode,filter_monoisotopic=filter_monoisotopic,monoisotopic_i=monoisotopic_i)))
}

#It probably should be included in the previous function

filter_fragments<-function(T2,tol=5,mode="pos")
{
  if(mode=="pos")
    fragmentation<-fragmentation_pos
  else
    fragmentation<-fragmentation_neg
  
  f<-NULL
  for(i in 1:nrow(T2)){
    if(any(withintol(T2[i,]$diff,abs(parse_fragmentation(fragmentation[T2[i,]$parental_lipid,T2[i,]$fragment_lipid])),tol)))
      f<-c(f,i)
  }
  
  return(T2[f,])
}

#Generate T3
generate_lipid_list<-function(lipid,synonyms="")
{
  if(synonyms!="")
    matches<-subset(df_metametabolits,grepl(lipid,nommetabolit)|(grepl(tolower(synonyms),tolower(nommetabolit))))
  else
    matches<-subset(df_metametabolits,grepl(lipid,nommetabolit))
  
  l<-data.frame(t(sapply(matches$nommetabolit,find_c_delta)))
  if(ncol(l)==2){
    names(l)<-c("c","delta")
    l$lipid<-lipid
    l$mz<-matches$monoisotopic_molecular_weight
  }
  else{
    l=NULL
  }
  return(unique(l))
  
}
generate_T3<-function(T1,T2)
{
  
  #Match paper mz
  MS2ID_lipids<-do.call(rbind,lapply(lipids$name,generate_lipid_list))
  
  tols<-outer(T1$experimental_mz,table$mz,get_tol)
  closest<-apply(tols,1,which.min)
  T1$paper_mz<-sapply(seq_along(closest),function(i) if(tols[i,closest[i]]>5) NA else table$mz[closest[i]])
  
  tols<-outer(T2_filters$fragment_experimental_mz....fragment.experimental_mz,table$mz,get_tol)
  closest<-apply(tols,1,which.min)
  T2_filters$paper_fragment_mz<-sapply(seq_along(closest),function(i) if(tols[i,closest[i]]>5) NA else table$mz[closest[i]])
  
  table$matched_mz<-table$mz%in%T1$paper_mz
  
  table$matched_id<-paste(table$mz,table$lipid,table$c,table$delta,tolower(table$adduct))%in%paste(T1$paper_mz,T1$lipid,T1$c,T1$delta,tolower(T1$adduct))
  
  table$matched_lipid<-paste(table$mz,table$lipid,table$c,table$delta)%in%paste(T1$paper_mz,T1$lipid,T1$c,T1$delta)
  
  table$lipid_in_MS2ID_old<-table$lipid_in_MS2ID
  
  table$lipid_in_MS2ID<-paste(table$lipid,table$c,table$delta) %in% paste(MS2ID_lipids$lipid,MS2ID_lipids$c,MS2ID_lipids$delta)
  
  subset(table,!matched_id&matched_mz)
  
  table$adduct_considered<-tolower(table$adduct)%in%tolower(subset(adducts,mode=="neg")$name)
  
  #Solving the missed matched that had an adduct considered (This is mostly a tolerance issue)
  s1<-subset(table,!matched_id&adduct_considered&matched_mz)
  subset(T1,paper_mz==s1$mz[2])
  
  table$matched_fragment_mz<-table$mz%in%T2_filters$paper_fragment_mz
  table$matched_fragment_id<-paste(table$mz,table$lipid,table$c,table$delta)%in%paste(T2_filters$paper_fragment_mz,T2_filters$parental_lipid,T2_filters$parental_c,T2_filters$parental_delta)
  
  write.csv(table,"/home/gbaquer/msidata/1. In-source Fragmentation/1.5. Validation Datasets/2. Nevus articulo fragmentación/results.csv")
  
  #Total Cases: 227
  nrow(table)
  
  #Different YELLOW cases 37 (16.3%)
  #YELLOW Case 1: Masses not accounted for 11 (4.8%)
  nrow(subset(table,!matched_mz))
  #YELLOW Case 2: Lipid not present in DB 26 (11.5%)
  nrow(subset(table,!lipid_in_MS2ID_old))
  
  #Case 2.1: No matches 2
  nrow(subset(table,is.na(lipid)))
  #Case 2.2: Found synonym 6
  nrow(subset(table,!lipid_in_MS2ID_old& lipid_in_MS2ID))
  #Case 2.3: Truly not found 18
  nrow(subset(table,!lipid_in_MS2ID&!is.na(lipid)))
  tmp<-subset(table,lipid_in_MS2ID&!is.na(lipid))
  tmp2<-(paste(tmp$lipid," ",tmp$c,":",tmp$delta,sep=""))
  length(tmp2)
  
  
  #Different GREEN cases:  
  
  
  #Different RED cases
  #RED Case 1: Adduct or fragment not accounted for
  #RED Case 2: Tolerance
  
  return(table)
}


generate_T1<-function(mz,tol=5,m="neg")
{
  ix<-1:length(mz)
  
  a<-subset(adducts,mode==m)
  if(m=="pos"){
    fragmentation<-fragmentation_pos
  }else{
    fragmentation<-fragmentation_neg
  }
  
  target_mz<-outer(outer(mz,adducts$mass,function(x,y)x+y),losses,function(x,y)x-y)
  results<-lapply(MS2ID_lipids$mz,function(x)withintol(x,target_mz,tol))
  return(results)
}
#Perform all DB searches [It can be optimized]

display_fragmentation_network<-function(T1,T2){
  library('igraph')
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
  library(igraph)
  
  # Create data
  data <- matrix(sample(0:1, 400, replace=TRUE, prob=c(0.8,0.2)), nrow=20)
  data <- (c>0.9)
  losses<-round(c(59.0734993,85.08914936 ),2)
  T2_visualisation<-subset(T2_possible_fragments,cor>0.7&round(diff,2)%in%losses)
  
  edgelist<-matrix("",nrow(T2_visualisation),2)
  edgelist[,1]<-unlist(lapply(1:nrow(T2_visualisation),function(i)paste(T2_visualisation[i,c('parental_lipid','parental_c','parental_delta')],collapse=":")))
  edgelist[,2]<-unlist(lapply(1:nrow(T2_visualisation),function(i)paste(T2_visualisation[i,c('fragment_lipid','fragment_c','fragment_delta')],collapse=":")))
  
  T1_db_matches[unique(match(T2_visualisation$parental_id,T1_db_matches$id)),c("i","adduct")]
  T1_db_matches[unique(match(T2_visualisation$fragment_id,T1_db_matches$id)),c("i","adduct")]
  
  
  network<-graph_from_edgelist(edgelist)
  plot(network, layout=layout.auto, main="fruchterman.reingold")
  
  # random graph
  net = rgraph(10, mode = "graph", tprob = 0.5)
  net = network(net, directed = FALSE)
  
  # vertex names
  network.vertex.names(net) = letters[1:10]
  ggnet2(net)
}

# DEPRECATED
# db_rm<-load_RefMet()
# plot(density(db_d$exactmass))
# lines(density(db_t$exactmass),col=2)
# lines(density(pks$mass),col=3)
# lines(density(db_rm$exactmass),col=4)
# lines(density(load_decoy_non_animal(c(500,1000))$exactmass),col=5)
# legend(2000,0.002,c("decoy","lipidmaps","pks","refmet","capped decoy"),col=1:5,lw=2)
# 
