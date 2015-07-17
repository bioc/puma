


 gmhta<-function(
     object
     ,background=FALSE
     ,gsnorm=c("median", "none", "mean", "meanlog")
    ,savepar=FALSE
     ,eps=1.0e-6
    ,addConstant = 0
     ,cl=NULL
     ,BatchFold=10
)
{  
   
    # show(date())
   #  library(hta20hsaceviewtcdf)
     if(!is.null(cl))
	     clusterEvalQ(cl, library(puma))
     library(pumadata);          
     chipnum<-length(sampleNames(object));

         HTA_transcript_NO <- NULL
            HTA_transcript_name <-    NULL  
              HTA_probes_transcripts <- NULL
             HTA_Location <-NULL
            rm(HTA_transcript_NO);
            rm(HTA_probes_transcripts);
            rm(HTA_transcript_name);
           rm(HTA_Location);

###load corresponding between gene and transcirpt and probes and alpha num of every gene
    
   
      cat('HTA type is Human');      
   

         data(HTA_transcript_NO); 
         Gene_T_NO =HTA_transcript_NO ;
               
         data(HTA_probes_transcripts);
         Probes_A_NUM = HTA_probes_transcripts;
       
         data(HTA_transcript_name)
         transcript_name <- HTA_transcript_name;          
         data(HTA_Location) 
         location=HTA_Location  

    
      
         

    cat('\n');
 ###find pm of every gene######




    
    All_index<-Gene_T_NO[,2]*2680+Gene_T_NO[,1]+1    ###pos_x and pos_y to indices
    
    Gene_T_INdex = cbind( Gene_T_NO[,3],All_index);  ###table of gene transcript index   

   
    unique_index  <- unique(All_index);
    
   if(background == TRUE)
   {
    for (i in c(1:chipnum)){
     
      m<-min(pm(object)[,i])
      pm(object)[,i]<-pm(object)[,i]-m+1
      
    }
  }

    pm_g <- oligo::pm(object,target='probeset');                           ##pm of cel 
   

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
   total_isoform<-dim(HTA_transcript_name)[1];
   #rm(pm);
  

   
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
  ####calculate pm_index and gt_index##########
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
       ################################################################


  
    cl_genes_num<-matrix(0,1,ncol=BatchFold*length(cl))
    cl_genes_num_temp = trunc(Total_genes/(BatchFold*length(cl)));
   
  
     for (j in c(1:(BatchFold*length(cl)-1))){
          cl_genes_num[1,j]<-cl_genes_num_temp
      }
     cl_genes_num[1,BatchFold*length(cl)]<-Total_genes-cl_genes_num_temp*(BatchFold*length(cl)-1)


       ################################################################




         
        
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


  if(!is.null(cl))
  {
	stopCluster(cl);
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
      expr_temp <- as.data.frame(2^expr_g)
      expr_temp <-as.matrix(expr_temp)
      expr_temp[apply(is.infinite(expr_temp), FUN = any, 1), ] = exp(700)
      chipm <- apply(expr_temp,2,mean)
      chipm <- chipm/chipm[1]

      #expr_g <- as.matrix(log2(expr_g))
      for (i in 1:chipnum)
      {
          expr_g[,i] <- expr_g[,i]-log2(chipm[i])
          se_g[,i]<- se_g[,i]/chipm[i]
          prc5_g[,i] <- prc5_g[,i]-log2(chipm[i])
          prc25_g[,i] <- prc25_g[,i]-log2(chipm[i])
          prc50_g[,i] <- prc50_g[,i]-log2(chipm[i])
          prc75_g[,i] <- prc75_g[,i]-log2(chipm[i])
          prc95_g[,i] <- prc95_g[,i]-log2(chipm[i])
 
    } 
    rm(expr_temp);

  }else if (gsnorm[1]=="median")
    {
       expr_g <- as.data.frame(2^expr_g)
    
       chipm <- apply(expr_g,2,median)
       chipm <- chipm/chipm[1]

        expr_g <- as.matrix(log2(expr_g))
        for (i in 1:chipnum)
        {
          expr_g[,i] <- expr_g[,i]-log2(chipm[i])
          se_g[,i]<- se_g[,i]/chipm[i]
          prc5_g[,i] <- prc5_g[,i]-log2(chipm[i])
          prc25_g[,i] <- prc25_g[,i]-log2(chipm[i])
          prc50_g[,i] <- prc50_g[,i]-log2(chipm[i])
          prc75_g[,i] <- prc75_g[,i]-log2(chipm[i])
          prc95_g[,i] <- prc95_g[,i]-log2(chipm[i])
        }
  }else if (gsnorm[1]=="meanlog")
  {
       chipm <- apply(expr_g,2,mean)
        chipm <- chipm-chipm[1]

       for (i in 1:chipnum)
       {
        if(i!=1)
        {
          expr_g[,i] <- expr_g[,i]-chipm[i]
          se_g[,i]<- se_g[,i]/chipm[i]
          prc5_g[,i] <- prc5_g[,i]-chipm[i]
          prc25_g[,i] <- prc25_g[,i]-chipm[i]
          prc50_g[,i] <- prc50_g[,i]-chipm[i]
          prc75_g[,i] <- prc75_g[,i]-chipm[i]
          prc95_g[,i] <- prc95_g[,i]-chipm[i]
        }
         
      }
  }





  if (gsnorm[1]=="mean")
  {
    expr_temp <- as.data.frame(2^expr_t)
    expr_temp <-as.matrix(expr_temp)
    expr_temp[apply(is.infinite(expr_temp), FUN = any, 1), ] = exp(700)
    chipm <- apply(expr_temp,2,mean)
    chipm <- chipm/chipm[1]

    #expr_t <- as.matrix(log2(expr_t))
    for (i in 1:chipnum)
    {
      expr_t[,i] <- expr_t[,i]-log2(chipm[i])
      se_t[,i]<- se_t[,i]/chipm[i]
      prc5_t[,i] <- prc5_t[,i]-log2(chipm[i])
      prc25_t[,i] <- prc25_t[,i]-log2(chipm[i])
      prc50_t[,i] <- prc50_t[,i]-log2(chipm[i])
      prc75_t[,i] <- prc75_t[,i]-log2(chipm[i])
      prc95_t[,i] <- prc95_t[,i]-log2(chipm[i])
    }
    rm(expr_temp);

  }else if (gsnorm[1]=="median")
  {
    expr_t <- as.data.frame(2^expr_t)
    
    chipm <- apply(expr_t,2,median)
    chipm <- chipm/chipm[1]

    expr_t <- as.matrix(log2(expr_t))
    for (i in 1:chipnum)
    {
      expr_t[,i] <- expr_t[,i]-log2(chipm[i])
      se_t[,i]<- se_t[,i]/chipm[i]
      prc5_t[,i] <- prc5_t[,i]-log2(chipm[i])
      prc25_t[,i] <- prc25_t[,i]-log2(chipm[i])
      prc50_t[,i] <- prc50_t[,i]-log2(chipm[i])
      prc75_t[,i] <- prc75_t[,i]-log2(chipm[i])
      prc95_t[,i] <- prc95_t[,i]-log2(chipm[i])
    }
  }else if (gsnorm[1]=="meanlog")
  {
    chipm <- apply(expr_t,2,mean)
    chipm <- chipm-chipm[1]

    for (i in 1:chipnum)
    {
      if(i!=1)
      {
        expr_t[,i] <- expr_t[,i]-chipm[i]
        se_t[,i]<- se_t[,i]/chipm[i]
        prc5_t[,i] <- prc5_t[,i]-chipm[i]
        prc25_t[,i] <- prc25_t[,i]-chipm[i]
        prc50_t[,i] <- prc50_t[,i]-chipm[i]
        prc75_t[,i] <- prc75_t[,i]-chipm[i]
        prc95_t[,i] <- prc95_t[,i]-chipm[i]
      }
      
    }
  }

 
  rownames(expr_t) <- transcript_name[,1];

  colnames(expr_t) <- sampleNames(object)
  rownames(se_t) <- transcript_name[,1];
  colnames(se_t) <- sampleNames(object)
 
  rownames(prc5_t) <- transcript_name[,1];
  colnames(prc5_t) <- sampleNames(object)
  rownames(prc25_t) <- transcript_name[,1];
  colnames(prc25_t) <- sampleNames(object)
 
  rownames(prc50_t) <- transcript_name[,1];
  colnames(prc50_t) <- sampleNames(object);

  rownames(prc75_t) <- transcript_name[,1];
  colnames(prc75_t) <- sampleNames(object)
 
  rownames(prc95_t) <- transcript_name[,1];
  colnames(prc95_t) <- sampleNames(object)
  
 rm(transcript_name);
 gc();  


  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)
  return_exprReslt_t <- new(
		"exprReslt"
	,	exprs=log2((2^expr_t)+addConstant)
	,	se.exprs=se_t
  	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData(object)
		,	varMetadata=data.frame(labelDescription=varLabels(phenoData(object)))
		)
	# , notes=notes
	)
          prcfive(return_exprReslt_t) <- prc5_t
	prctwfive(return_exprReslt_t) <- prc25_t
	prcfifty(return_exprReslt_t) <- prc50_t
	prcsevfive(return_exprReslt_t) <- prc75_t
	prcninfive(return_exprReslt_t) <- prc95_t
  
  annotation(return_exprReslt_t) <- annotation(object)
  description(return_exprReslt_t) <- description(object)
  notes(return_exprReslt_t) <- notes(object)
 



  rownames(expr_g) <-rownames(Probes_A_NUM);
  colnames(expr_g) <- sampleNames(object)
  rownames(se_g) <- rownames(Probes_A_NUM);
  colnames(se_g) <- sampleNames(object)
  
  rownames(prc5_g) <- rownames(Probes_A_NUM);
  colnames(prc5_g) <- sampleNames(object)
  rownames(prc25_g) <- rownames(Probes_A_NUM);
  colnames(prc25_g) <- sampleNames(object)
  rownames(prc50_g) <- rownames(Probes_A_NUM);
  colnames(prc50_g) <- sampleNames(object)
  rownames(prc75_g) <- rownames(Probes_A_NUM);
  colnames(prc75_g) <- sampleNames(object)
  rownames(prc95_g) <- rownames(Probes_A_NUM);
  colnames(prc95_g) <- sampleNames(object)




  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)

  return_exprReslt_g <- new(
		"exprReslt"
	,	exprs=log2((2^expr_g)+addConstant)
	,	se.exprs=se_g
  	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData(object)
		,	varMetadata=data.frame(labelDescription=varLabels(phenoData(object)))
		)
	# , notes=notes
	)
          prcfive(return_exprReslt_g) <- prc5_g
         	prctwfive(return_exprReslt_g) <- prc25_g
	prcfifty(return_exprReslt_g) <- prc50_g
	prcsevfive(return_exprReslt_g) <- prc75_g
	prcninfive(return_exprReslt_g) <- prc95_g
  

  annotation(return_exprReslt_g) <- annotation(object)
  description(return_exprReslt_g) <- description(object)
  notes(return_exprReslt_g) <- notes(object)

  rm(object);
  gc();

 

 # show(date())
  return_exprReslt<-list(gene=return_exprReslt_g,transcript=return_exprReslt_t)
  #save(return_exprReslt,file="return_exprReslt.RDATA");
  return (return_exprReslt)



}
