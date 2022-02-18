#API - main functions exported
env <- new.env(parent = emptyenv())

updateEnv <-function(x){
  for (i in 1:length(x)){
    env[[names(x)[i]]]<-x[[i]]
  }
}

#' Load fragmentation and adduct data
#'
#' Creates a plot of the crayon colors in \code{\link{brocolors}}
#'
#' @return A list containing lipis, lipiddict, adducts, fragmentation pos, fragmentation neg, and losses
#'
#' @examples
#' loadData()
#'
#' @export
loadData<-function(){
  lipids<-data.frame(name=c("DG","TG","PA","PE","PC","PG","PI","PS","CL","Cer","CerP","SM","HexCer","SFT","GM3","FRAG"),stringsAsFactors = F)
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

  return(list(lipids=lipids,lipid_dict=lipid_dict,adducts=adducts,fragmentation_pos=fragmentation_pos,fragmentation_neg=fragmentation_neg,losses=losses))
}
#' Search adducts and fragments in database
#'
#' @param mz list of masses in the dataset
#' @param db database to search
#' @param tol tolerance in ppm. Default is 5ppm.
#' @param m ionization mode ("pos" or "neg"). Default is positive ion mode ("pos").
#' @param rand_af boolean indicating whether to use random adducts and fragments (for target-decoy validation method 1). During normal operation it should be set to False. Default False.
#' @param db_ref reference db (for target-decoy validation method 2). Default is db (No alternative reference db)
#'
#' @return a data.frame of class "rMSILipidAnnot". Containing information about all hits. Including theoretical mz, molecular formula, compound name, lipid family and the ranking scores.
#'
#' @examples
#' Missing examples
#'
#' @export
search_db<-function(mz,db,tol=5,m="pos",rand_af=F,db_ref=db,d){
  #Generate a and f
  db_mz<-db$exactmass
  a<-subset(env$adducts,mode==m)$mass
  names(a)<-subset(env$adducts,mode==m)$name
  if(m=="neg")
    frag<-env$fragmentation_neg
  else
    frag<-env$fragmentation_pos
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
  results$adduct_mode<-subset(env$adducts,mode==m)$mode[results$adduct_i]
  results$adduct_score<-sapply(1:nrow(results),function(i)subset(env$adducts,mode==m)[results$adduct_i[i],as.character(results$lipid_ref[i])])
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

#' Add extra ranking scores
#'
#' @param res data.frame of class "rMSILipidAnnot". Usually obtained from rMSIfragment::search_db().
#' @param pks peak matrix of the MSI experiment. Usually obtained from rMSIproc
#' @param m ionization mode ("pos" or "neg"). Default is positive ion mode ("pos").
#'
#' @return a data.frame of class "rMSILipidAnnot". With updated ranking scores.
#'
#' @examples
#' Missing examples
#'
#' @export
ranking<-function(res,pks,m,d){
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
    #possible_adducts<-rownames(subset(env$adducts,mode==m))[which(subset(env$adducts,mode==m)[,res$lipid_ref[match(i,res$lipid_id)]]!=0)]
    #res$parental_percentage[res$lipid_id==i]<-sum(possible_adducts%in%p_adduct)/length(possible_adducts)

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
#' Perform annotation (combined call to search_db and ranking)
#'
#' @param pks peak matrix of the MSI experiment. Usually obtained from rMSIproc
#' @param db database to search
#' @param tol tolerance in ppm. Default is 5ppm.
#' @param m ionization mode ("pos" or "neg"). Default is positive ion mode ("pos").
#' @param rand_af boolean indicating whether to use random adducts and fragments (for target-decoy validation method 1). During normal operation it should be set to False. Default False.
#' @param db_ref reference db (for target-decoy validation method 2). Default is db (No alternative reference db)
#'
#' @return a data.frame of class "rMSILipidAnnot". With the results of annotation.
#'
#' @examples
#' Missing examples
#'
#' @export
annotate<-function(pks,db,tol=5,m="pos",rand_af=F,db_ref=db){
  res<-search_db(pks$mass,db,tol,m,rand_af,db_ref)
  res_rank<-ranking(subset(res,adduct_mode==m&fragmentation_possible),pks,m)
  return(res_rank)
}

#' Manual validation
#'
#' @export
manual_validation<-function(res,table){
  res$table_mz<-sapply(res$experimental_mz,function(x)with(table,mz[which.min(abs(mz-x))]))
  a<-with(table,paste(lipid,c,delta,mz))
  b<-with(res,paste(lipid,c,delta,table_mz))
  res$match_HPLC_adduct<-sapply(match(b,a),function(i)if(is.na(i))""else as.character(table$HPLC_adduct[i]))
  res$match<-b%in%a
  table$match<-a%in%b

  return(list(res=res,table=table))
}
