 igmoExon<-function(
     cel.path
     ,SampleNameTable
     ,exontype = c("Human", "Mouse", "Rat")
     ,background=FALSE
     ,gsnorm=c("median", "none", "mean", "meanlog")
     ,savepar=FALSE
     ,eps=1.0e-6
     ,addConstant = 0
     ,condition=c("Yes","No")
     ,cl=NULL
     ,BatchFold=10
)
{     
    if(!is.null(cl))
	     clusterEvalQ(cl, library(puma))
    library(pumadata);

   if(missing(cel.path)){
	stop("'cel.path' is missing! Please give the path to your CEL files!")
     }
   
   if(missing(SampleNameTable)){
	stop("'SampleNameTable' is missing! Please give the path of the your sample table!")
     }


   setwd(cel.path);
   desc <- read.table(SampleNameTable,header=TRUE,sep="\t",stringsAsFactors=FALSE)

   celnames<-desc[,1];

   condition_all<-desc[,2];

   object_condition_names<-c();

  

    expr_g_all<-c()
     se_g_all<-c()
     prc5_g_all<-c()
     prc25_g_all<-c()
     prc50_g_all<-c()
     prc75_g_all<-c()
     prc95_g_all<-c()


     expr_t_all<-c()
     se_t_all<-c()
     prc5_t_all<-c()
     prc25_t_all<-c()
     prc50_t_all<-c()
     prc75_t_all<-c()
     prc95_t_all<-c()
  
    
for(i in c(1:max(condition_all))){   

   

      celname_condition=celnames[which(condition_all==i)]
  
   
      object_condition<-read.celfiles(celname_condition);

    if(condition[1]=="Yes"){
         object_condition_names<-c(object_condition_names,celname_condition); 
    }

    if(condition[1]=="Yes"&&i==max(condition_all)){
         object<-read.celfiles(object_condition_names);
    }
   if(condition[1]!="Yes"){
        object<-object_condition;
    }
   
    chipnum<-length(sampleNames(object_condition));
    chipnum_all<-length(celnames);
        
              Human_transcript_NO <- NULL
            Human_transcript_name <-    NULL  
              Human_probes_transcripts <- NULL
             Human_Location <-NULL
            rm(Human_transcript_NO);
            rm(Human_probes_transcripts);
            rm(Human_transcript_name);
           rm(Human_Location);

              Mouse_transcript_NO <- NULL
            Mouse_transcript_name <-    NULL  
              Mouse_probes_transcripts <- NULL
             Mouse_Location <-NULL
            rm(Mouse_transcript_NO);
            rm(Mouse_probes_transcripts);
            rm(Mouse_transcript_name);
           rm(Mouse_Location);

        Rat_transcript_NO <- NULL
            Rat_transcript_name <-    NULL  
              Rat_probes_transcripts <- NULL
             Rat_Location <-NULL
            rm(Rat_transcript_NO);
            rm(Rat_probes_transcripts);
            rm(Rat_transcript_name);
            rm(Rat_Location);

###load corresponding between gene and transcirpt and probes and alpha num of every gene
    
    if(exontype[1]=="Human")                     
       {   
            cat('Exon type is Human');      

            data(Human_transcript_NO); 
            Gene_T_NO =Human_transcript_NO ;
               
            data(Human_probes_transcripts);
            Probes_A_NUM = Human_probes_transcripts;
       
            data(Human_transcript_name)
            transcript_name <- Human_transcript_name;          
           
            data(Human_Location);
            location=Human_Location;
            
         
       }else if(exontype[1]=="Mouse")
         {  
           cat('Exon type is Mouse');      

            data(Mouse_transcript_NO); 
            Gene_T_NO =Mouse_transcript_NO ;

            data(Mouse_probes_transcripts);
            Probes_A_NUM = Mouse_probes_transcripts;

            data(Mouse_transcript_name)
            transcript_name <- Mouse_transcript_name;
            
            data(Mouse_Location);
            location=Mouse_Location;
              

        } else if(exontype[1]=="Rat")
        {
            cat('Exon type is Rat');      

            data(Rat_transcript_NO); 
            Gene_T_NO =Rat_transcript_NO ;
               
            data(Rat_probes_transcripts);
            Probes_A_NUM = Rat_probes_transcripts;
         
            data(Rat_transcript_name)
            transcript_name <- Rat_transcript_name; 

            data(Rat_Location) 
            location=Rat_Location

       }
   

    cat('\n');
 ###find pm of every gene######




    #All_index <- xy2indices(Gene_T_NO[,1],Gene_T_NO[,2],abatch=object_condition);    ###pos_x and pos_y to indices
    All_index <-Gene_T_NO[,2]*2560+Gene_T_NO[,1]+1
    Gene_T_INdex = cbind( Gene_T_NO[,3],All_index);  ###table of gene transcript index   

   
    unique_index  <- unique(All_index);
    
   if(background == TRUE)
   {
    for (i in c(1:chipnum)){
     
      m<-min(oligo::pm(object_condition,target='probeset')[,i])
      oligo::pm(object_condition,target='probeset')[,i]<-oligo::pm(object_condition,target='probeset')[,i]-m+1
      
    }
  }

    pm_g <- oligo::pm(object_condition,target='probeset');                           ##pm of cel 
   

    Length_unique_index <- length(unique_index);
      
   if(chipnum==1)
   {
       pm <- pm_g[location]         ##nofind pm use 1 to substitute for one chips
       for (j in c(1:Length_unique_index))
       {
           if(is.na(pm[j]))
             pm[j] <- 1;
       }
   }else
   {  
       pm <- pm_g[location,]         ##nofind pm use 1 to substitute for multi chips
       for (j in c(1:Length_unique_index))
       {
           if(is.na(pm[j,2]))
           {            
             pm[j,] <- 1;
            }
       }

   }

  
  
   rm(location)
   rm(Gene_T_NO)

   pm <- cbind(unique_index,pm);  ##Table of PM
   rm(unique_index)

   total_probes_everygene <- Probes_A_NUM[,1];    ##total num of probes for everygene
   total_isoform<-dim(Human_transcript_name)[1];
   rm(pm_g);
  

   
   st <- 0;
   ed <- 0;

   Total_genes <- dim(Probes_A_NUM)[1];  ##total num of genes
   unique_probes_everygene <- seq(1,Total_genes,by=1);  

   for(i in c(1:Total_genes ))
   { 
      ed <- st + total_probes_everygene[i];
      st <- st + 1;
     
      unique_probes_everygene[i] <- length(unique(All_index[st:ed]));
     
      st<-ed;
   }

  alpha_num_everygene <- Probes_A_NUM[,2];  ###alpha num of everygene

  rm(All_index)

  gc();

  ####optimise#####
 prctiles <- 0.01*c(5, 25, 50, 75, 95);
  len_pc<-length(prctiles);
   pm_index <- c();
   gt_index <- c();
   pmst=0;
   gtst=0;
    for(i in c(1:Total_genes ))
   { 
      pm_index[i] <- pmst;
      pmst <- pmst + unique_probes_everygene[i];
      
      gt_index[i] <- gtst;
      
      gtst <- gtst+ total_probes_everygene[i];
      
   }
   pm_index[Total_genes+1]<-pmst
   gt_index[Total_genes+1]<-gtst
   
   expectedCompletionTime <-	system.time(
								res <-
  	      							.Call(
  	  								"gme_c"
         								, pm
         								, Gene_T_INdex
         								, unique_probes_everygene
         								, total_probes_everygene
         								, alpha_num_everygene
         								, 10                        
                        , prctiles
  	                    , len_pc
         								, savepar
         								, eps
  	 								)
							)[3];
   if(is.null(cl))
   {
      expectedCompletionTime <- expectedCompletionTime*Total_genes/10;
   }else
   {
      expectedCompletionTime <- expectedCompletionTime*Total_genes/10/length(cl);
   }
   if (expectedCompletionTime < 120){
		expectedTimeString <- paste(round(expectedCompletionTime, 0), "seconds")
	}else if (expectedCompletionTime < 7200){
		expectedTimeString <- paste(round(expectedCompletionTime/60, 0), "minutes")
	}else if (expectedCompletionTime < 172800){
		expectedTimeString <- paste(round(expectedCompletionTime/3600, 0), "hours")
	}else
	     {	expectedTimeString <- paste(round(expectedCompletionTime/172800, 0), "days")
             }	
    cat(paste("optimise expected completion time is", expectedTimeString,"\n"))


   if(is.null(cl))
   {
        ####optimise#####
        res <-
  	      .Call(
  	  	"gme_c"
         	, pm
         	, Gene_T_INdex
         	, unique_probes_everygene
         	, total_probes_everygene
         	, alpha_num_everygene
         	, Total_genes
                , prctiles
  	       , len_pc
         	, savepar
         	, eps
  	 	)
   }else
   {
           cl_genes_num<-matrix(0,1,ncol=BatchFold*length(cl))
          cl_genes_num_temp = trunc(Total_genes/(BatchFold*length(cl)));
   
  
     for (j in c(1:(BatchFold*length(cl)-1))){
          cl_genes_num[1,j]<-cl_genes_num_temp
      }
     cl_genes_num[1,BatchFold*length(cl)]<-Total_genes-cl_genes_num_temp*(BatchFold*length(cl)-1)


		paramsList<-list();   
		if(length(cl)>1)
		{
			for (i in 1:(BatchFold*length(cl)))
			{
                                
				paramsList[[i]] <- list(
				 pm[(pm_index[sum(cl_genes_num[1,0:(i-1)])+1]+1):pm_index[sum(cl_genes_num[1,0:i])+1],]			###pm
				, Gene_T_INdex[(gt_index[sum(cl_genes_num[1,0:(i-1)])+1]+1):gt_index[sum(cl_genes_num[1,0:i])+1],]	###Gene_T_INdex
				, unique_probes_everygene[(sum(cl_genes_num[1,0:(i-1)])+1):(sum(cl_genes_num[1,0:i]))]    		###unique_probes_everygene
         			, total_probes_everygene[(sum(cl_genes_num[1,0:(i-1)])+1):(sum(cl_genes_num[1,0:i]))]      		###total_probes_everygene
         			, alpha_num_everygene[(sum(cl_genes_num[1,0:(i-1)])+1):(sum(cl_genes_num[1,0:i]))]        		###alpha_num_everygene
			  	, cl_genes_num[1,i]											###total_genes
				);
				
			}
		}
		
		cat('start parallel compute...');	   
            temp_res <- clusterApplyLB(
						 cl
						,paramsList
						,function(
                                                  params        
                                                                , prctiles
  	                                                        , len_pc
								, savepar
         							, eps
								) 
							{
     							 	res <- .Call(
  	  							"gme_c"
         							, params[[1]]                 ###pm
         							, params[[2]]               	###Gene_T_INdex
         							, params[[3]]   			###unique_probes_everygene
         							, params[[4]]     		###total_probes_everygene
         							, params[[5]]       		###alpha_num_everygene
          							, params[[6]] 
                                                                , prctiles
  	                                                        , len_pc              	
          							, savepar                    	###savepar
          							, eps                        	###eps
     	 							)
							}
                                                , prctiles
  	                                        , len_pc
						, savepar
         					, eps
						)
            

            cl_isoform_num<-c()
            for(i in 1:(BatchFold*length(cl))){

               cl_isoform_num[i]<-(length(temp_res[[i]][,1])-cl_genes_num[1,i]*7)/7

            }

            expr_g<-c();
            se_g<-c();
            expr_t<-c();
            se_t<-c();
            prc5_g<-c();
            prc25_g<-c();
            prc50_g<-c();
            prc75_g<-c();
            prc95_g<-c();
            prc5_t<-c();
            prc25_t<-c();
            prc50_t<-c();
            prc75_t<-c();
            prc95_t<-c();
            res_t<-c();
		for(i in 1:(BatchFold*length(cl)))
  		{
			expr_g <- rbind(expr_g,temp_res[[i]][1:cl_genes_num[1,i],])
                        se_g <- rbind(se_g,temp_res[[i]][(cl_genes_num[1,i]+1):(2*cl_genes_num[1,i]),])
                       prc5_g <- rbind(prc5_g,temp_res[[i]][(2*cl_genes_num[1,i]+1):(3*cl_genes_num[1,i]),])
                        prc25_g <- rbind(prc25_g,temp_res[[i]][(3*cl_genes_num[1,i]+1):(4*cl_genes_num[1,i]),])
                        prc50_g <- rbind(prc50_g,temp_res[[i]][(4*cl_genes_num[1,i]+1):(5*cl_genes_num[1,i]),])
                        prc75_g <- rbind(prc75_g,temp_res[[i]][(5*cl_genes_num[1,i]+1):(6*cl_genes_num[1,i]),])
                        prc95_g <- rbind(prc95_g,temp_res[[i]][(6*cl_genes_num[1,i]+1):(7*cl_genes_num[1,i]),])

                        
                         expr_t<-rbind(expr_t,temp_res[[i]][(cl_genes_num[1,i]*7+1):(cl_genes_num[1,i]*7+cl_isoform_num[i]),])
                        se_t<-rbind(se_t,temp_res[[i]][(1+ cl_genes_num[1,i]*7+cl_isoform_num[i]):(cl_genes_num[1,i]*7+2*cl_isoform_num[i]),])
                        prc5_t<-rbind(prc5_t,temp_res[[i]][(1+ cl_genes_num[1,i]*7+2*cl_isoform_num[i]):(cl_genes_num[1,i]*7+3*cl_isoform_num[i]),])
                        prc25_t<-rbind(prc25_t,temp_res[[i]][(1+ cl_genes_num[1,i]*7+3*cl_isoform_num[i]):(cl_genes_num[1,i]*7+4*cl_isoform_num[i]),])
                        prc50_t<-rbind(prc50_t,temp_res[[i]][(1+ cl_genes_num[1,i]*7+4*cl_isoform_num[i]):(cl_genes_num[1,i]*7+5*cl_isoform_num[i]),])
                        prc75_t<-rbind(prc75_t,temp_res[[i]][(1+ cl_genes_num[1,i]*7+5*cl_isoform_num[i]):(cl_genes_num[1,i]*7+6*cl_isoform_num[i]),])
                        prc95_t<-rbind(prc95_t,temp_res[[i]][(1+ cl_genes_num[1,i]*7+6*cl_isoform_num[i]):(cl_genes_num[1,i]*7+7*cl_isoform_num[i]),])
    
 		}
     
  }




  if(is.null(cl)){

      res_g<-res[1:(Total_genes*7),]
     res_t<-res[(Total_genes*7+1):((Total_genes+total_isoform)*7),]

     expr_g <- matrix(res_g[c(1:Total_genes),],Total_genes,chipnum);
     se_g <- matrix(res_g[c((Total_genes+1):(2*Total_genes)),],Total_genes,chipnum);
     prc5_g <- matrix(res_g[c((2*Total_genes+1):(3*Total_genes)),],Total_genes,chipnum)
     prc25_g <- matrix(res_g[c((3*Total_genes+1):(4*Total_genes)),],Total_genes,chipnum)
     prc50_g <- matrix(res_g[c((4*Total_genes+1):(5*Total_genes)),],Total_genes,chipnum)
     prc75_g <- matrix(res_g[c((5*Total_genes+1):(6*Total_genes)),],Total_genes,chipnum)
     prc95_g <- matrix(res_g[c((6*Total_genes+1):(7*Total_genes)),],Total_genes,chipnum)


     expr_t <- matrix(res_t[c(1:total_isoform),],total_isoform,chipnum);
     se_t <- matrix(res_t[c((total_isoform+1):(2*total_isoform)),],total_isoform,chipnum);
     prc5_t <- matrix(res_t[c((2*total_isoform+1):(3*total_isoform)),],total_isoform,chipnum)
     prc25_t <- matrix(res_t[c((3*total_isoform+1):(4*total_isoform)),],total_isoform,chipnum)
     prc50_t <- matrix(res_t[c((4*total_isoform+1):(5*total_isoform)),],total_isoform,chipnum)
     prc75_t <- matrix(res_t[c((5*total_isoform+1):(6*total_isoform)),],total_isoform,chipnum)
     prc95_t <- matrix(res_t[c((6*total_isoform+1):(7*total_isoform)),],total_isoform,chipnum)

    }
    

     expr_g_all<-cbind(expr_g_all,expr_g)
     se_g_all<-cbind(se_g_all,se_g)
     prc5_g_all<-cbind(prc5_g_all,prc5_g)
     prc25_g_all<-cbind(prc25_g_all,prc25_g)
     prc50_g_all<-cbind(prc50_g_all,prc50_g)
     prc75_g_all<-cbind(prc75_g_all,prc75_g)
     prc95_g_all<-cbind(prc95_g_all,prc95_g)

     expr_t_all<-cbind(expr_t_all,expr_t)
     se_t_all<-cbind(se_t_all,se_t)
     prc5_t_all<-cbind(prc5_t_all,prc5_t)
     prc25_t_all<-cbind(prc25_t_all,prc25_t)
     prc50_t_all<-cbind(prc50_t_all,prc50_t)
     prc75_t_all<-cbind(prc75_t_all,prc75_t)
     prc95_t_all<-cbind(prc95_t_all,prc95_t)

 }## done 

  if(!is.null(cl))
  {
	stopCluster(cl);
  }


 

 rm(pm);
 rm(Gene_T_INdex);
 

 rm(res); 
  gc();
   ###save result#####

   ###gene_names <- rownames(Probes_A_NUM);


   cat("\n")
   cat("Gene and transcirpt expression values are returned.");

   cat('Done.')
 

 
 
  if (gsnorm[1]=="mean")
   {
      expr_g_all <- as.data.frame(2^expr_g_all)
     
      chipm <- apply(expr_g_all,2,mean)
      chipm <- chipm/chipm[1]

      expr_g_all <- as.matrix(log2(expr_g_all))
      for (i in 1:chipnum)
      {
          expr_g_all[,i] <- expr_g_all[,i]-log2(chipm[i])
          se_g_all[,i]<- se_g[,i]/chipm[i]
          prc5_g_all[,i] <- prc5_g_all[,i]-log2(chipm[i])
          prc25_g_all[,i] <- prc25_g_all[,i]-log2(chipm[i])
          prc50_g_all[,i] <- prc50_g_all[,i]-log2(chipm[i])
          prc75_g_all[,i] <- prc75_g_all[,i]-log2(chipm[i])
          prc95_g_all[,i] <- prc95_g_all[,i]-log2(chipm[i])
    }
  }else if (gsnorm[1]=="median")
    {
       expr_g_all <- as.data.frame(2^expr_g_all)
    
       chipm <- apply(expr_g_all,2,median)
       chipm <- chipm/chipm[1]

        expr_g_all <- as.matrix(log2(expr_g_all))
        for (i in 1:chipnum)
        {
          expr_g_all[,i] <- expr_g_all[,i]-log2(chipm[i])
          se_g_all[,i]<- se_g[,i]/chipm[i]
          prc5_g_all[,i] <- prc5_g_all[,i]-log2(chipm[i])
          prc25_g_all[,i] <- prc25_g_all[,i]-log2(chipm[i])
          prc50_g_all[,i] <- prc50_g_all[,i]-log2(chipm[i])
          prc75_g_all[,i] <- prc75_g_all[,i]-log2(chipm[i])
          prc95_g_all[,i] <- prc95_g_all[,i]-log2(chipm[i])
        }
  }else if (gsnorm[1]=="meanlog")
  {
       chipm <- apply(expr_g_all,2,mean)
        chipm <- chipm-chipm[1]

       for (i in 1:chipnum)
       {
         expr_g_all[,i] <- expr_g_all[,i]-chipm[i]
         se_g_all[,i]<- se_g[,i]/chipm[i]
         prc5_g_all[,i] <- prc5_g_all[,i]-chipm[i]
         prc25_g_all[,i] <- prc25_g_all[,i]-chipm[i]
         prc50_g_all[,i] <- prc50_g_all[,i]-chipm[i]
         prc75_g_all[,i] <- prc75_g_all[,i]-chipm[i]
         prc95_g_all[,i] <- prc95_g_all[,i]-chipm[i]
      }
  }





  if (gsnorm[1]=="mean")
  {
    expr_t_all <- as.data.frame(2^expr_t_all)
    chipm <- apply(expr_t_all,2,mean)
    chipm <- chipm/chipm[1]

    expr_t_all <- as.matrix(log2(expr_t_all))
    for (i in 1:chipnum)
    {
      expr_t_all[,i] <- expr_t_all[,i]-log2(chipm[i])
      se_t_all[,i]<- se_t[,i]/chipm[i]
      prc5_t_all[,i] <- prc5_t_all[,i]-log2(chipm[i])
      prc25_t_all[,i] <- prc25_t_all[,i]-log2(chipm[i])
      prc50_t_all[,i] <- prc50_t_all[,i]-log2(chipm[i])
      prc75_t_all[,i] <- prc75_t_all[,i]-log2(chipm[i])
      prc95_t_all[,i] <- prc95_t_all[,i]-log2(chipm[i])
    }
  }
 
 else if (gsnorm[1]=="median")
  {
    expr_t_all <- as.data.frame(2^expr_t_all)
    
    chipm <- apply(expr_t_all,2,median)
    chipm <- chipm/chipm[1]

    expr_t_all <- as.matrix(log2(expr_t_all))
    for (i in 1:chipnum)
    {
      expr_t_all[,i] <- expr_t_all[,i]-log2(chipm[i])
      se_t_all[,i]<- se_t[,i]/chipm[i]
      prc5_t_all[,i] <- prc5_t_all[,i]-log2(chipm[i])
      prc25_t_all[,i] <- prc25_t_all[,i]-log2(chipm[i])
      prc50_t_all[,i] <- prc50_t_all[,i]-log2(chipm[i])
      prc75_t_all[,i] <- prc75_t_all[,i]-log2(chipm[i])
      prc95_t_all[,i] <- prc95_t_all[,i]-log2(chipm[i])
    }
  }
  else if (gsnorm[1]=="meanlog")
  {
    chipm <- apply(expr_t_all,2,mean)
    chipm <- chipm-chipm[1]

    for (i in 1:chipnum)
    {
      expr_t_all[,i] <- expr_t_all[,i]-chipm[i]
      se_t_all[,i]<- se_t[,i]/chipm[i]
      prc5_t_all[,i] <- prc5_t_all[,i]-chipm[i]
      prc25_t_all[,i] <- prc25_t_all[,i]-chipm[i]
      prc50_t_all[,i] <- prc50_t_all[,i]-chipm[i]
      prc75_t_all[,i] <- prc75_t_all[,i]-chipm[i]
      prc95_t_all[,i] <- prc95_t_all[,i]-chipm[i]
    }
  }

 
  rownames(expr_t_all) <- transcript_name[,1];

  colnames(expr_t_all) <- sampleNames(object)
  rownames(se_t_all) <- transcript_name[,1];
  colnames(se_t_all) <- sampleNames(object)
 
  rownames(prc5_t_all) <- transcript_name[,1];
  colnames(prc5_t_all) <- sampleNames(object)
  rownames(prc25_t_all) <- transcript_name[,1];
  colnames(prc25_t_all) <- sampleNames(object)
 
  rownames(prc50_t_all) <- transcript_name[,1];
  colnames(prc50_t_all) <- sampleNames(object);

  rownames(prc75_t_all) <- transcript_name[,1];
  colnames(prc75_t_all) <- sampleNames(object)
 
  rownames(prc95_t_all) <- transcript_name[,1];
  colnames(prc95_t_all) <- sampleNames(object)
  
 rm(transcript_name);
 gc();  


  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)
  return_exprReslt_t <- new(
		"exprReslt"
	,	exprs=log2((2^expr_t_all)+addConstant)
	,	se.exprs=se_t_all
  	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData(object)
		,	varMetadata=data.frame(labelDescription=varLabels(phenoData(object)))
		)
	# , notes=notes
	)
          prcfive(return_exprReslt_t) <- prc5_t_all
	prctwfive(return_exprReslt_t) <- prc25_t_all
	prcfifty(return_exprReslt_t) <- prc50_t_all
	prcsevfive(return_exprReslt_t) <- prc75_t_all
	prcninfive(return_exprReslt_t) <- prc95_t_all
  
  annotation(return_exprReslt_t) <- annotation(object)
  description(return_exprReslt_t) <- description(object)
  notes(return_exprReslt_t) <- notes(object)
 



  rownames(expr_g_all) <-rownames(Probes_A_NUM);
  colnames(expr_g_all) <- sampleNames(object)
  rownames(se_g_all) <- rownames(Probes_A_NUM);
  colnames(se_g_all) <- sampleNames(object)
  
  rownames(prc5_g_all) <- rownames(Probes_A_NUM);
  colnames(prc5_g_all) <- sampleNames(object)
  rownames(prc25_g_all) <- rownames(Probes_A_NUM);
  colnames(prc25_g_all) <- sampleNames(object)
  rownames(prc50_g_all) <- rownames(Probes_A_NUM);
  colnames(prc50_g_all) <- sampleNames(object)
  rownames(prc75_g_all) <- rownames(Probes_A_NUM);
  colnames(prc75_g_all) <- sampleNames(object)
  rownames(prc95_g_all) <- rownames(Probes_A_NUM);
  colnames(prc95_g_all) <- sampleNames(object)




  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)

  return_exprReslt_g <- new(
		"exprReslt"
	,	exprs=log2((2^expr_g_all)+addConstant)
	,	se.exprs=se_g_all
  	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData(object)
		,	varMetadata=data.frame(labelDescription=varLabels(phenoData(object)))
		)
	# , notes=notes
	)
          prcfive(return_exprReslt_g) <- prc5_g_all
         	prctwfive(return_exprReslt_g) <- prc25_g_all
	prcfifty(return_exprReslt_g) <- prc50_g_all
	prcsevfive(return_exprReslt_g) <- prc75_g_all
	prcninfive(return_exprReslt_g) <- prc95_g_all
  

  annotation(return_exprReslt_g) <- annotation(object)
  description(return_exprReslt_g) <- description(object)
  notes(return_exprReslt_g) <- notes(object)

  rm(object);
  gc();




  return_exprReslt<-list(gene=return_exprReslt_g,transcript=return_exprReslt_t)
  
  return (return_exprReslt)



}

