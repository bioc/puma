


 gmoExon<-function(
     object
     ,exontype = c("Human", "Mouse", "Rat")
     ,GT = c("gene", "transcript")
     ,background=FALSE
     ,gsnorm=c("median", "none", "mean", "meanlog")
     ,savepar=FALSE
     ,eps=1.0e-6
     ,addConstant = 0
)
{  
     library(pumadata);          
     chipnum<-length(sampleNames(object));

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

            data(Human_Location) 
            location=Human_Location
            
         
       }else if(exontype[1]=="Mouse")
         {  
           cat('Exon type is Mouse');      

            data(Mouse_transcript_NO); 
            Gene_T_NO =Mouse_transcript_NO ;

            data(Mouse_probes_transcripts);
            Probes_A_NUM = Mouse_probes_transcripts;

            data(Mouse_transcript_name)
            transcript_name <- Mouse_transcript_name;
          
            data(Mouse_Location) 
            location=Mouse_Location      

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




    All_index <- xy2indices(Gene_T_NO[,1],Gene_T_NO[,2],abatch=object);    ###pos_x and pos_y to indices
    
    Gene_T_INdex = cbind( Gene_T_NO[,3],All_index);  ###table of gene transcript index   

   
    unique_index  <- unique(All_index);
    
   if(background == TRUE)
   {
    for (i in c(1:chipnum)){
     
      m<-min(pm(object)[,i])
      pm(object)[,i]<-pm(object)[,i]-m+1
      
    }
  }

    pm_g <- pm(object);                           ##pm of cel 
   

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
  res <-
  	.Call(
  	  "gme_c"
         , pm
         , Gene_T_INdex
         , unique_probes_everygene
         , total_probes_everygene
         , alpha_num_everygene
         , Total_genes
         , savepar
         , eps
         , PACKAGE="puma"

  	 )
  cat('Done.')

  cat("Expression values calculating");
  g <- c();
  g_NO <- c();
  for(i in c(1:Total_genes))
    {
          g_NO <- rep(i,alpha_num_everygene[i]*chipnum+2);
          g <- c(g,g_NO);
          
            
    }
  
   Alpha_cd <- cbind(g,res);
 
   
   rm(res);
   gc();
      
   lk=0;
   lg=0;
   lt=0;
   isoform_means<-c()
   isoform_std<-c()
   gene_means<-c()
   gene_std<-c()
   prctiles <- 0.01*c(5, 25, 50, 75, 95);
   gene_pr5<-c();
   gene_pr25<-c();
   gene_pr50<-c();
   gene_pr75<-c();
   gene_pr95<-c();
   isoform_pr5<-c();
   isoform_pr25<-c();
   isoform_pr50<-c();
   isoform_pr75<-c();
   isoform_pr95<-c();


###exprs and std of transcript ########## 
    
   
   for(g in (1:Total_genes)) 
     {
         
             isoform_cur_means <- array();
             isoform_cur_std <- array();
            
            
 
             cond <- alpha_num_everygene[g]*chipnum;
     

             cur_index <- c();               
             cur_pm <- c();                   ###pm of current processed gene
             index_start <- 1+lg;
             index_end <- lg+total_probes_everygene[g];


             cur_index <-Gene_T_INdex[index_start:index_end,];     ###index of current processed gene

             if(index_end-index_start==0)
              {
                cur_index=t(cur_index)
               }     
     
             pm_start <- 1+lt;
             pm_end <- unique_probes_everygene[g]+lt;
             cur_pm <- pm[pm_start:pm_end,];
             if(pm_end-pm_start==0)
             {
               cur_pm=t(cur_pm);
             } 

        
        
              lg <- total_probes_everygene[g]+lg;           ###next gene 
              lt <- unique_probes_everygene[g]+lt;
     

              cur_pos_alpha_cd <- array(0,dim=c(1,cond+2));
              cur_pos_alpha_cd <- which(Alpha_cd[,1]==g);    ###posite of alpha ,c and d for current processed gene
              cur_pos_alpha <-cur_pos_alpha_cd[1:cond];          ###posite of alpha for current processed gene
  
              cur_alpha <- Alpha_cd[cur_pos_alpha,2];       ###alpha data for current processed gene
    
              cur_pos_c <- max(cur_pos_alpha)+1;       ###c and d  data for current processed gene
              cur_pos_d <- max(cur_pos_alpha)+2;
              cur_c <- Alpha_cd[cur_pos_c,2];
              cur_d <- Alpha_cd[cur_pos_d,2];






         metrix_probe_alpha <- array(0,dim=c(unique_probes_everygene[g], alpha_num_everygene[g]));  ### corresponding alpha for every probe in one chip##
  
         mutlimetrix_probe_alpha <- array(0,dim=c(unique_probes_everygene[g],cond));       ###corresponding alpha for every probe in multi chips##
  
        for(i in c(1:unique_probes_everygene[g]))
        {

            
             
                 for(j in c(1:total_probes_everygene[g]))
                 {
                      
                     if(cur_index[j,2]==cur_pm[i,1])
                        metrix_probe_alpha[i,cur_index[j,1]] <- 1;
                 }
      
               mutlimetrix_probe_alpha[i,] <- rep(metrix_probe_alpha[i,],chipnum);
             
           
        } 


    ###calculate transcript values ########

        D1 <- array(0,dim=c(unique_probes_everygene[g],cond));
        h_isoform <- array(0,dim=c(cond,cond));

        for(j in c(1:unique_probes_everygene[g]))
          {
                 ###covariance metrix of alpha for current gene###
              MC <- c();
              h1 <- array();
              real_alpha <- cur_alpha*mutlimetrix_probe_alpha[j,];    ###the real alpha for current gene
      
              temp <- array(0,dim=c(alpha_num_everygene[g],alpha_num_everygene[g]));
              for(mal in c(1:alpha_num_everygene[g]))
              {
                   if(metrix_probe_alpha[j,mal]==1)
                        temp[mal,] <- metrix_probe_alpha[j,];

              }
     

              MC1 <- c();
               for(kt in c(1:chipnum))
               {
                  MC1 <- cbind(MC1,temp);
               }
          
               for(kt in c(1:chipnum))
               {
                  MC <- rbind(MC,MC1);
               }
      

              q <- sum(real_alpha)+cur_c;
              w_end <- chipnum+1;
              w <- cur_d+sum(cur_pm[j,2:w_end]);
       
              AM <- array(real_alpha,dim=c(alpha_num_everygene[g],chipnum));
              H <- array(trigamma(q),dim=c(cond,cond));          ###covariance metrix of alpha for current gene###
      
              k_end <- chipnum-1;
              for(k in c(0:k_end))
              {   
                  almid_start <- 1+k;
                  total_alpha_onechip <- sum(AM[,almid_start]);
                  dif1_alpha_start <- 2+k;
                  dif1_alpha <- log(cur_pm[j, dif1_alpha_start])-log(w)-digamma(total_alpha_onechip)+digamma(q);   ###First Derivative of alpha
               
               
                  dif2_alpha <- trigamma(q)-trigamma(total_alpha_onechip);
                  Dif1 <- rep(dif1_alpha,alpha_num_everygene[g]);
              
                  D1_start <- 1+k*alpha_num_everygene[g];
                  D1_end <- (k+1)*alpha_num_everygene[g];
                  D1[j,D1_start:D1_end] <- Dif1;
                  H[D1_start:D1_end,D1_start:D1_end] <- dif2_alpha;         ###second derivative of alpha
              }
              h1 <- H*MC;
              h_isoform <- h_isoform+h1;
         }
  
        Q <- mutlimetrix_probe_alpha*D1;
        G <- array(0,dim=c(cond,1));
  
        for(jk in c(1:cond))
         {
            G[jk,1] <- sum(Q[,jk]);             ###the  devirative of current gene ###
         }


    
         MapH <- 1/(-diag(h_isoform));
         MU <- MapH%*%G+cur_alpha;
       
                

          
        TC_muy = array(0,dim=c(1,cond));
        TC_sigy = array(0,dim = c(1,cond));
        TC_kk_isoform = array(0,dim = c(1,cond));


      
      
        for(i in c(1:cond))
        {
             sig_isoform <- MapH[i];
            
             if(cur_alpha[i]-(1e-06)>0.00001)
             {
                   mu_isoform <- cur_alpha[i];
             }else
                 {
                      mu_isoform = MU[i];
                 }
                   


             kk_isoform<-2/(1+(pnorm(mu_isoform/sqrt(sig_isoform))/0.5-1));
        
             muy_isoform = kk_isoform*(sqrt(sig_isoform)*exp(-mu_isoform^2/(2*sig_isoform))/sqrt(2*pi)
                            +mu_isoform*(1+(pnorm(mu_isoform/sqrt(sig_isoform))/0.5-1))/2);  
             
             sigy_isoform = kk_isoform*((sig_isoform+(mu_isoform-muy_isoform)^2)*(1+(pnorm(mu_isoform/sqrt(sig_isoform))/0.5-1))/2
                            +sqrt(sig_isoform/(2*pi))*exp(-mu_isoform^2/(2*sig_isoform))*(mu_isoform-2*muy_isoform)); 
           
             TC_muy[i] = muy_isoform;
             TC_sigy[i] = sigy_isoform;  
             TC_kk_isoform[i]=kk_isoform;
         }

         Mu_TC = matrix(TC_muy,alpha_num_everygene[g],chipnum);
         Sig_TC = matrix(TC_sigy,alpha_num_everygene[g],chipnum);
         kk_TC = matrix(TC_kk_isoform,alpha_num_everygene[g],chipnum);

         
         
        if(GT[1] == 'gene')
        {

           e_m_gene <- array(0,dim=c(1,chipnum));
           std_gene <- array(0,dim=c(1,chipnum));

           pr5_gene <- array(0,dim=c(1,chipnum));
           pr25_gene <- array(0,dim=c(1,chipnum));
           pr50_gene <- array(0,dim=c(1,chipnum));
           pr75_gene <- array(0,dim=c(1,chipnum));
           pr95_gene <- array(0,dim=c(1,chipnum));




         for(i in c(1:chipnum))
            {
            
                mu_gene = sum(Mu_TC[,i]);
                sig_gene = sum(Sig_TC[,i]);
              

                kk_gene <-2/(1+(pnorm(mu_gene/sqrt(sig_gene))/0.5-1));
  
                e_m_gene[i] = digamma(mu_gene)+log(cur_d)-digamma(cur_c)+psigamma(mu_gene,2)*sig_gene/2;
                std_gene[i] = sqrt(psigamma(mu_gene,1)^2*sig_gene);
        ##calculate standard deviation#### 

             ##5%
	   
                gene5_temp_qnorm = 1-2.0*(1.0-prctiles[1])/kk_gene;
                gene5_temp = mu_gene+sqrt(2.0*sig_gene)*(-qnorm((1-gene5_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
           
                pr5_gene[i]= (digamma(gene5_temp)+log(cur_d)-digamma(cur_c))/log(2.0);
		
            ###25%
                gene25_temp_qnorm = 1-2.0*(1.0-prctiles[2])/kk_gene;
                gene25_temp = mu_gene+sqrt(2.0*sig_gene)*(-qnorm((1-gene25_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                pr25_gene[i]= (digamma(gene25_temp)+log(cur_d)-digamma(cur_c))/log(2.0);
           
            ###50
                gene50_temp_qnorm = 1-2.0*(1.0-prctiles[3])/kk_gene;
                gene50_temp = mu_gene+sqrt(2.0*sig_gene)*(-qnorm((1-gene50_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                pr50_gene[i]= (digamma(gene50_temp)+log(cur_d)-digamma(cur_c))/log(2.0);
            
            ###75
                gene75_temp_qnorm = 1-2.0*(1.0-prctiles[4])/kk_gene;
                gene75_temp = mu_gene+sqrt(2.0*sig_gene)*(-qnorm((1-gene75_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                pr75_gene[i]= (digamma(gene75_temp)+log(cur_d)-digamma(cur_c))/log(2.0);
            
            ####95
                gene95_temp_qnorm = 1-2.0*(1.0-prctiles[5])/kk_gene;
                gene95_temp = mu_gene+sqrt(2.0*sig_gene)*(-qnorm((1-gene95_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                pr95_gene[i]= (digamma(gene95_temp)+log(cur_d)-digamma(cur_c))/log(2.0);           

              }
              
               gene_means <-rbind(gene_means,e_m_gene);
               gene_std <-rbind(gene_std,std_gene);  
               gene_pr5 <- rbind(gene_pr5,pr5_gene); 
               gene_pr25 <- rbind(gene_pr25,pr25_gene); 
               gene_pr50 <- rbind(gene_pr50,pr50_gene); 
               gene_pr75 <- rbind(gene_pr75,pr75_gene);  
               gene_pr95 <- rbind(gene_pr95,pr95_gene);
          }



           
          
           
        if(GT[1] == 'transcript')
          {
              e_m_isoform <- matrix(0,alpha_num_everygene[g],chipnum);
              std_isoform <- matrix(0,alpha_num_everygene[g],chipnum);
              pr5_isoform <-matrix(0,alpha_num_everygene[g],chipnum);
              pr25_isoform <-matrix(0,alpha_num_everygene[g],chipnum);
              pr50_isoform <-matrix(0,alpha_num_everygene[g],chipnum); 
              pr75_isoform <-matrix(0,alpha_num_everygene[g],chipnum);
              pr95_isoform <-matrix(0,alpha_num_everygene[g],chipnum);

            for(i in c(1:alpha_num_everygene[g]))
             {
              for(j in c(1:chipnum))
               {
                
                 kk_isoform = kk_TC[i,j];
                 e_m_isoform[i,j] = digamma(Mu_TC[i,j])+log(cur_d)-digamma(cur_c)+psigamma(Mu_TC[i,j],2)*Sig_TC[i,j]/2;
                 std_isoform[i,j] = sqrt(psigamma(Mu_TC[i,j],1)^2*Sig_TC[i,j]);
                 
                 isoform5_temp_qnorm = 1-2.0*(1.0-prctiles[1])/kk_isoform;
                 isoform5_temp = Mu_TC[i,j]+sqrt(2.0*Sig_TC[i,j])*(-qnorm((1-isoform5_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                 pr5_isoform[i,j] = (digamma(isoform5_temp)+log(cur_d)-digamma(cur_c))/log(2.0);
  
             ####25
                 isoform25_temp_qnorm = 1-2.0*(1.0-prctiles[2])/kk_isoform;
                 isoform25_temp = Mu_TC[i,j]+sqrt(2.0*Sig_TC[i,j])*(-qnorm((1-isoform25_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                 pr25_isoform[i,j] = (digamma(isoform25_temp)+log(cur_d)-digamma(cur_c))/log(2.0);


            ####50
                isoform50_temp_qnorm = 1-2.0*(1.0-prctiles[3])/kk_isoform;
                isoform50_temp = Mu_TC[i,j]+sqrt(2.0*Sig_TC[i,j])*(-qnorm((1-isoform50_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                pr50_isoform[i,j] = (digamma(isoform50_temp)+log(cur_d)-digamma(cur_c))/log(2.0);

            ####75
                isoform75_temp_qnorm = 1-2.0*(1.0-prctiles[4])/kk_isoform;
                isoform75_temp = Mu_TC[i,j]+sqrt(2.0*Sig_TC[i,j])*(-qnorm((1-isoform75_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
                pr75_isoform[i,j] = (digamma(isoform75_temp)+log(cur_d)-digamma(cur_c))/log(2.0);

            ###95 
               isoform95_temp_qnorm = 1-2.0*(1.0-prctiles[5])/kk_isoform;
               isoform95_temp = Mu_TC[i,j]+sqrt(2.0*Sig_TC[i,j])*(-qnorm((1-isoform95_temp_qnorm)/2,0,1,TRUE,FALSE)/sqrt(2));
               pr95_isoform[i,j] = (digamma(isoform95_temp)+log(cur_d)-digamma(cur_c))/log(2.0);
         
               }
            }
         
          isoform_means <-rbind(isoform_means,e_m_isoform);
  
          isoform_std <-rbind(isoform_std,std_isoform);

          isoform_pr5 <-rbind(isoform_pr5,pr5_isoform);
          isoform_pr25 <-rbind(isoform_pr25,pr25_isoform);
          isoform_pr50 <-rbind(isoform_pr50,pr50_isoform);
          isoform_pr75 <-rbind(isoform_pr75,pr75_isoform);
          isoform_pr95 <-rbind(isoform_pr95,pr95_isoform);

       }
    
     if(g%%1000==0) 
      {
      cat('.');
     
      }       
    
     
 }  
   
      cat('Done.')
     


 rm(pm);
 rm(Gene_T_INdex);
 rm(Alpha_cd);
 gc();
   ###save result#####

   ###gene_names <- rownames(Probes_A_NUM);

if(GT[1]=="gene")
 { 
   cat("\n")
   cat("Gene expression values are returned.");

   expr <- gene_means
   se <- gene_std
   
   prc5 <- gene_pr5;
   prc25 <- gene_pr25;
   prc50 <- gene_pr50;
   prc75 <- gene_pr75;
   prc95 <- gene_pr95;
 
   rm(gene_means)
   rm(gene_std)
   rm(gene_pr5)
   rm(gene_pr50)
   rm(gene_pr75)
   rm(gene_pr95)
   gc();
 

  if (gsnorm[1]=="mean")
   {
      expr <- as.data.frame(2^expr)
     
      chipm <- apply(expr,2,mean)
      chipm <- chipm/chipm[1]

       expr <- as.matrix(log2(expr))
      for (i in 1:chipnum)
      {
          expr[,i] <- expr[,i]-log2(chipm[i])
          prc5[,i] <- prc5[,i]-log2(chipm[i])
          prc25[,i] <- prc25[,i]-log2(chipm[i])
          prc50[,i] <- prc50[,i]-log2(chipm[i])
          prc75[,i] <- prc75[,i]-log2(chipm[i])
          prc95[,i] <- prc95[,i]-log2(chipm[i])
    }
  }else if (gsnorm[1]=="median")
    {
       expr <- as.data.frame(2^expr)
    
       chipm <- apply(expr,2,median)
       chipm <- chipm/chipm[1]

        expr <- as.matrix(log2(expr))
        for (i in 1:chipnum)
        {
          expr[,i] <- expr[,i]-log2(chipm[i])
          prc5[,i] <- prc5[,i]-log2(chipm[i])
          prc25[,i] <- prc25[,i]-log2(chipm[i])
          prc50[,i] <- prc50[,i]-log2(chipm[i])
          prc75[,i] <- prc75[,i]-log2(chipm[i])
          prc95[,i] <- prc95[,i]-log2(chipm[i])
        }
  }else if (gsnorm[1]=="meanlog")
  {
       chipm <- apply(expr,2,mean)
        chipm <- chipm-chipm[1]

       for (i in 1:chipnum)
       {
         expr[,i] <- expr[,i]-chipm[i]
         prc5[,i] <- prc5[,i]-chipm[i]
         prc25[,i] <- prc25[,i]-chipm[i]
         prc50[,i] <- prc50[,i]-chipm[i]
         prc75[,i] <- prc75[,i]-chipm[i]
         prc95[,i] <- prc95[,i]-chipm[i]
      }
  }




  rownames(expr) <-rownames(Probes_A_NUM);
  colnames(expr) <- sampleNames(object)
  rownames(se) <- rownames(Probes_A_NUM);
  colnames(se) <- sampleNames(object)
  
  rownames(prc5) <- rownames(Probes_A_NUM);
  colnames(prc5) <- sampleNames(object)
  rownames(prc25) <- rownames(Probes_A_NUM);
  colnames(prc25) <- sampleNames(object)
  rownames(prc50) <- rownames(Probes_A_NUM);
  colnames(prc50) <- sampleNames(object)
  rownames(prc75) <- rownames(Probes_A_NUM);
  colnames(prc75) <- sampleNames(object)
  rownames(prc95) <- rownames(Probes_A_NUM);
  colnames(prc95) <- sampleNames(object)

  rm(Probes_A_NUM);
  gc();

  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)

  return_exprReslt <- new(
		"exprReslt"
	,	exprs=log2((2^expr)+addConstant)
	,	se.exprs=se
  	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData(object)
		,	varMetadata=data.frame(labelDescription=varLabels(phenoData(object)))
		)
	# , notes=notes
	)
          prcfive(return_exprReslt) <- prc5
         	prctwfive(return_exprReslt) <- prc25
	prcfifty(return_exprReslt) <- prc50
	prcsevfive(return_exprReslt) <- prc75
	prcninfive(return_exprReslt) <- prc95
  

  annotation(return_exprReslt) <- annotation(object)
  description(return_exprReslt) <- description(object)
  notes(return_exprReslt) <- notes(object)
  rm(object);
  gc();
  return_exprReslt


} else if(GT[1] == "transcript")
{ 
  cat("\n") 
  cat("Transcript expression values are returned.")

  

   expr <- isoform_means
   se <- isoform_std
 
   prc5 <- isoform_pr5;
   prc25 <- isoform_pr25;
   prc50 <- isoform_pr50;
   prc75 <- isoform_pr75;
   prc95 <- isoform_pr95;
   
   rm(isoform_means)
   rm(isoform_std)
   rm(isoform_pr5)
   rm(isoform_pr25)
   rm(isoform_pr95)
   rm(isoform_pr75)
   gc();

  if (gsnorm[1]=="mean")
  {
    expr <- as.data.frame(2^expr)
    chipm <- apply(expr,2,mean)
    chipm <- chipm/chipm[1]

    expr <- as.matrix(log2(expr))
    for (i in 1:chipnum)
    {
      expr[,i] <- expr[,i]-log2(chipm[i])
      prc5[,i] <- prc5[,i]-log2(chipm[i])
      prc25[,i] <- prc25[,i]-log2(chipm[i])
      prc50[,i] <- prc50[,i]-log2(chipm[i])
      prc75[,i] <- prc75[,i]-log2(chipm[i])
      prc95[,i] <- prc95[,i]-log2(chipm[i])
    }
  }
 
 else if (gsnorm[1]=="median")
  {
    expr <- as.data.frame(2^expr)
    
    chipm <- apply(expr,2,median)
    chipm <- chipm/chipm[1]

    expr <- as.matrix(log2(expr))
    for (i in 1:chipnum)
    {
      expr[,i] <- expr[,i]-log2(chipm[i])
      prc5[,i] <- prc5[,i]-log2(chipm[i])
      prc25[,i] <- prc25[,i]-log2(chipm[i])
      prc50[,i] <- prc50[,i]-log2(chipm[i])
      prc75[,i] <- prc75[,i]-log2(chipm[i])
      prc95[,i] <- prc95[,i]-log2(chipm[i])
    }
  }
  else if (gsnorm[1]=="meanlog")
  {
    chipm <- apply(expr,2,mean)
    chipm <- chipm-chipm[1]

    for (i in 1:chipnum)
    {
      expr[,i] <- expr[,i]-chipm[i]
      prc5[,i] <- prc5[,i]-chipm[i]
      prc25[,i] <- prc25[,i]-chipm[i]
      prc50[,i] <- prc50[,i]-chipm[i]
      prc75[,i] <- prc75[,i]-chipm[i]
      prc95[,i] <- prc95[,i]-chipm[i]
    }
  }

 
  rownames(expr) <- transcript_name[,1];

  colnames(expr) <- sampleNames(object)
  rownames(se) <- transcript_name[,1];
  colnames(se) <- sampleNames(object)
 
  rownames(prc5) <- transcript_name[,1];
  colnames(prc5) <- sampleNames(object)
  rownames(prc25) <- transcript_name[,1];
  colnames(prc25) <- sampleNames(object)
 
  rownames(prc50) <- transcript_name[,1];
  colnames(prc50) <- sampleNames(object);

  rownames(prc75) <- transcript_name[,1];
  colnames(prc75) <- sampleNames(object)
 
  rownames(prc95) <- transcript_name[,1];
  colnames(prc95) <- sampleNames(object)
  
 rm(transcript_name);
 gc();  


  phenodata <- phenoData(object)
  annotation <- annotation(object)
  description <- description(object)
  return_exprReslt <- new(
		"exprReslt"
	,	exprs=log2((2^expr)+addConstant)
	,	se.exprs=se
  	,	phenoData = new(
			"AnnotatedDataFrame"
		,	data=pData(object)
		,	varMetadata=data.frame(labelDescription=varLabels(phenoData(object)))
		)
	# , notes=notes
	)
          prcfive(return_exprReslt) <- prc5
	prctwfive(return_exprReslt) <- prc25
	prcfifty(return_exprReslt) <- prc50
	prcsevfive(return_exprReslt) <- prc75
	prcninfive(return_exprReslt) <- prc95
  
  annotation(return_exprReslt) <- annotation(object)
  description(return_exprReslt) <- description(object)
  notes(return_exprReslt) <- notes(object)
  rm(object);
  gc();
  return_exprReslt 

}

}

