# utild.R : Helper functions
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
simply<-function(x){
  if(is.null(x$correlation))
    return((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","adduct_score")]))
  else
    return((x[,c("experimental_mz","ppm_error","abbreviation","formula","adduct","fragmentation","lipid_occurences","correlation","adduct_score")]))

}
cor_score<-function(m){
  if(sum(m!=1)==0)
    return(0)
  else
    return(mean(m[m!=1]))
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
topN<-function(d,s="S",n=3)
{
  return(do.call("rbind",lapply(unique(d[["mz_i"]]),function(i)head(subset(d[order(d[[s]],decreasing = T),],mz_i==i),n))))
}
