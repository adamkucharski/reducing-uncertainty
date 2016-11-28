# Code to accompany:
# Reducing uncertainty about flavivirus infections. Lancet Infect Dis. 2016
# (c) Kucharski AJ, Riley S

# - - - - - - - - - -
# Model code

generate_timeseries<-function(strains,inf_years,n_part=20,p.inf.in=0.2,sd.val.in=0.5){

  inf.n=length(inf_years)
  #Set per year incidence, to create correlation between participant infection histories
  time_series=NULL
  for(ii in 1:strains){
    if(length(p.inf.in)>1){p.inf=p.inf.in[ii]}else{p.inf}
    if(length(sd.val.in)>1){sd.val=sd.val.in[ii]}else{sd.val=sd.val.in}
    log.sd=sd.val
    attack.yr=sapply(rlnorm(inf.n,meanlog=log(p.inf)-log.sd^2/2,sdlog=log.sd),function(x){min(x,1)})
    time_series=cbind(time_series,attack.yr)
  }
  time_series
}

c.nume<-function(x,bdy=0.95){
  bp1=c(median(x),quantile(x,0.5-bdy/2),quantile(x,0.5+bdy/2))
  as.numeric(bp1)
}

c.text<-function(x,sigF=2){
  if(length(x)>3){
    bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975)),sigF)
  }else{
    bp1=signif(x,sigF)
  }
  paste(bp1[1]," (95% CrI: ",bp1[2],"--",bp1[3],")",sep="")
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Simulate serological data

simulate_sera_data<-function(strains,inf.years.sera=c(1980:2015),time.series.in=NULL,
                             theta=c(error=0.05,sigma=0.4),p.inf.in=0.2,sd.val.in=1,seedi=1,
                             roundv=F,dmatrix.in=NULL,zikv.attack=0.5,
                             per_sample=50){ # ii=participant | jj=test year
  
  set.seed(seedK)
  
  # theta guide:
  # error = probability negative given homotypic infection (i.e. 1 - sensitivity)
  # sigma = probability positive given heterotypic infection (i.e. 1 - specificity)
  # - - - - - - 
  # Set year of birth - sample uniformly
  #age.yr=sort(sample(1:50,n_part,replace = TRUE))
  inf.n=length(inf.years.sera)
  age.yr=c(0:(inf.n-1))
  n_part=inf.n

  # Generate annual attack rates
  if(is.null(time.series.in)){
    time.series=generate_timeseries(strains,inf.years.sera,n_part,p.inf.in,sd.val.in)
    
    # circulation in last year only for final strain, with attack rate = zikv.attack
    time.series[,strains]=0
    time.series[c((inf.n)),strains]=zikv.attack
    
  }else{
    time.series=time.series.in
  }
  
  # - - - - - - 
  # Simulate random infection history for each participant
  
  historytabSim=NULL
  for(ii in 1:n_part){
    # Identify which years individual was alive for (ensure new borns can be infected)
    mask.alive=c(max(1,inf.n-(age.yr[ii])):inf.n)
    # Calculate probability of infection with each strain
    if(length(mask.alive)>1){
      prob.inf=1-apply(1-time.series[mask.alive,], 2, prod)
    }else{
      prob.inf=1-sapply(1-time.series[mask.alive,], prod)
    }
    historytabSim=rbind(historytabSim,prob.inf)
  }

  # Define cross reaction matrix
  if(is.null(dmatrix.in)){
    dmatrix=diag(strains)*(1-theta[["error"]])+(1-diag(strains))*theta[["sigma"]] # default cross-reaction structure
  }else{
    dmatrix=dmatrix.in
  }
  
  # - - - - - - 
  # Simulate assay results for each participant
  
  test.list.sera=NULL
  
  for(ii in 1:n_part){
    prob.inf=historytabSim[ii,]
  
    # Probability test positive to each strain
    prob.positive = 1-apply(1-t(dmatrix*prob.inf),1,prod)
    test.list.sera = rbind(test.list.sera,prob.positive)
  }

  # - - - - - - 
  # Simulate binomial results for each participant
  
  n_sample=matrix(per_sample,nrow=strains,ncol=length(inf.years.sera))
  titre.prob=t(test.list.sera)
  
  prob_inf0=t(time.series) # Fix initial conditions to correct history
  prob_inf0=as.matrix(prob_inf0[,c(length(prob_inf0[1,]):1)])
  
  titre.data=NULL
  for(jj in 1:strains){
    tstrain=NULL
    for(kk in 1:length(inf.years.sera)){
      rand=rbinom(1,size=n_sample[jj,kk],prob=titre.prob[jj,kk])/per_sample
      tstrain=cbind(tstrain,rand)
    }
    titre.data=rbind(titre.data,tstrain)
  }
  
  #test.list.sera
  
  # - - - - - - 
  # Export data
  save(test.list.sera,titre.data,n_sample,prob_inf0,historytabSim,inf.years.sera,strains,n_part,time.series,age.yr,theta,file=paste("R_datasets/Sero_sim_",seedi,"_",strains,".RData",sep=""))

}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot simulated serological data

f.y<-function(x){x} #{rev(x)} # swap order of ages?

plot_simulated_sera_data<-function(seedi,strains){

  load(paste("R_datasets/Sero_sim_",seedi,"_",strains,".RData",sep=""))
  
  par(mar = c(3,3,1,1))
  par(mgp=c(1.8,0.6,0))
  par(mfrow=c(1,3))
  
  col.list=list(col1=rgb(0.9,0.6,0),col2=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
  label.age=seq(0,length(age.yr),2)
  lw.1=1.5
  
  # Plot attack rates
  plot(inf.years.sera,time.series[,1],type="l",ylim=c(0,1.01),col=col.list$col1,xlab="year",ylab="attack rate",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(inf.years.sera,time.series[,2],type="l",col=col.list$col2,lwd=lw.1)
  lines(inf.years.sera,time.series[,3],type="l",col=col.list$col3,lwd=lw.1)
  lines(inf.years.sera,time.series[,4],type="l",col=col.list$col4,lwd=lw.1)
  lines(inf.years.sera,time.series[,5],type="l",col=col.list$col5,lwd=lw.1)
  #grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  title(main=LETTERS[1],adj=0)

  # Plot probability of infection by age
  plot(age.yr,f.y(historytabSim[,1]),ylim=c(0,1.01),xlim=c(0,length(historytabSim[,1])),col=col.list$col1,xlab="age in 2015",ylab="probability ever infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(age.yr,f.y(historytabSim[,2]),col=col.list$col2,lwd=lw.1)
  lines(age.yr,f.y(historytabSim[,3]),col=col.list$col3,lwd=lw.1)
  lines(age.yr,f.y(historytabSim[,4]),col=col.list$col4,lwd=lw.1)
  lines(age.yr,f.y(historytabSim[,5]),col=col.list$col5,lwd=lw.1)
  #axis(1,at=label.age,labels=f.y(label.age))
  #grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  title(main=LETTERS[2],adj=0)
  
  # Plot probability of seropositivity by age
  plot(age.yr,f.y(test.list.sera[,1]),ylim=c(0,1.01),xlim=c(0,length(test.list.sera[,1])),col=col.list$col1,xlab="age in 2015",ylab="probability seropositive",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(age.yr,f.y(test.list.sera[,2]),col=col.list$col2,lwd=lw.1)
  lines(age.yr,f.y(test.list.sera[,3]),col=col.list$col3,lwd=lw.1)
  lines(age.yr,f.y(test.list.sera[,4]),col=col.list$col4,lwd=lw.1)
  lines(age.yr,f.y(test.list.sera[,5]),col=col.list$col5,lwd=lw.1)
  #axis(1,at=label.age,labels=f.y(label.age))
  title(main=LETTERS[3],adj=0)
  #grid(ny = NULL, nx = 0, col = rgb(0.9,0.9,0.9), lty = "solid")
  
  dev.copy(pdf,paste("plot_simulations/serology_plot.pdf",sep=""),width=8,height=2)
  dev.off()
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# MCMC Functions

SampleTheta<-function(theta_initial,m,covartheta,nparam,n_part){
  
  # sample from multivariate normal distribution - no adaptive sampling
  theta_star = as.numeric(exp(mvrnorm(1,log(theta_initial), Sigma=covartheta)))
  names(theta_star)=names(theta_initial)
  
  if(!is.na(match("error",names(theta_star)))){ # Check whether fitting theta or attack vector
    mu1=min(2-theta_star[["sigma"]],theta_star[["sigma"]])
    theta_star[["sigma"]]=ifelse(mu1<0,theta_initial[["sigma"]],mu1)
    
    mu1=min(2-theta_star[["error"]],theta_star[["error"]])
    theta_star[["error"]]=ifelse(mu1<0,theta_initial[["error"]],mu1)
  }else{
    #print("check")
    # Ensure FoI is below 3 (i.e. <95% probability infection per year)
    theta_check=sapply(theta_star,function(x){min(6-x,x)})
    theta_star=theta_check*(as.numeric(theta_check)>0) + theta_initial*(as.numeric(theta_check)<0)
  }
  
  theta_star
}

ComputeProbability<-function(marg_likelihood,marg_likelihood_star){
  # Flat priors on theta => symmetric update probability
  calc.lik = exp(marg_likelihood_star-marg_likelihood)
  min(1, calc.lik)
}


LikelihoodTitres<-function(titre.data,dmatrix,forcetab_star,inf.n,strains,n_sample){
  
  lik=NULL

  for(ii in 1:inf.n){ # Iterate across years
    p.inf=NULL
    for(jj in 1:strains){ # Iterate across strains
      p.inf=c(p.inf,1-exp(-sum(forcetab_star[jj,1:ii]))) # Infected in this period with strain j
    }
    # Calculate serological likelihood
    prob.positive = 1-apply(1-t(dmatrix*p.inf),1,prod)
    lik=rbind(lik,dbinom(round(n_sample[,ii]*titre.data[,ii]),size=n_sample[,ii],prob=prob.positive,log=T))
  }
  
  sum(lik)
  
  #colSums(lik)
}

# DEBUGGING
#dmatrix = diag(strains)*(1-tail(theta_post[["error"]],1))+(1-diag(strains))*tail(theta_post[["sigma"]],1) # default cross-reaction structure
#dmatrix = diag(strains)*(1-0.000001)+(1-diag(strains))*0.000001 # default cross-reaction structure
##post_force[5,2]=0.2
#forcetab_star=post_force
#LikelihoodTitres(titre.data,dmatrix,forcetab_star,inf.n,strains,n_sample) # DEBUG

#plot(llA[,4])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Run MCMC function

run_mcmc<-function(
  titre.data, # Note that this starts present and goes back into past (i.e. ordered increasing with age). Data as proportions
  n_sample, # Total samples in each age group and each strain
  inf_years,
  theta,
  strains=5,
  prob_inf=NULL, # initial conditions
  force_constrain, # matrix to constrain circulation years - Note that this is increasing age
  switch1=2,
  runs,
  seedi=1,
  pmask=NULL
){
  
  # Check inputs same size
  if(sum(dim(titre.data))!=sum(dim(n_sample))){return("titre.data and n_sample need to be same size")}
  
  if(sum(pmask=="sigma")>0){ theta[["sigma"]] = 0 } # Fix sigma 0 if hidden
  if(sum(pmask=="error")>0){ theta[["error"]] = 0 } # Fix error 0 if hidden
  
  inf.n=length(inf_years)
  # Preallocate memory
  nparam=length(theta); npcov=rep(1,nparam); npcov[match(pmask,names(theta))]=0 # mask specified parameters
  cov_matrix_theta0 = diag(npcov)
  cov_matrix_force0 = diag(rep(1,inf.n))
  
  thetatab=matrix(NA,nrow=(runs+1),ncol=length(theta)); colnames(thetatab)=names(theta)
  thetatab[1,]=theta

  if(!is.null(prob_inf)){
    # Convert attack rate to FoI
    #diff=as.matrix(prob_inf)-as.matrix(cbind(rep(0,strains),prob_inf[,1:(length(prob_inf[1,])-1)]))
    forcetab=-log(1-prob_inf)+0.001 # Add to ensure mixing if zero in simulation
  }else{
    forcetab=matrix(0.1,nrow=strains,ncol=inf.n)
    forcetab[strains,1]=1 # Increase FOI for ZIKV
  }
  forcetab=force_constrain*forcetab
  forcetabCollect=forcetab

  # Preallocate matrices
  likelihoodtab=matrix(-Inf,nrow=(runs+1))
  accepttabT=NULL
  accepttabH=NULL
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Run MCMC
  
  for (m in 1:runs){
    
    # Adaptive covariance matrix
    if(m==1){
      epsilon0=0.01
      varpart_prob0=0.01
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_force=varpart_prob0*cov_matrix_force0
    }else{
      epsilon0=max(1e-5,min(1,exp(log(epsilon0)+(accept_rateT-0.234)*0.999^m)))
      varpart_prob0=max(1e-5,min(1,exp(log(varpart_prob0)+(accept_rateH-0.234)*0.999^m))) # force of infection sampling
      
      cov_matrix_theta=epsilon0*cov_matrix_theta0
      cov_matrix_force=varpart_prob0*cov_matrix_force0  # ***DEBUGGING***
    }
    
    # - - - - - - - - - - - - - - - -
    # Resample parameters
    
    if(m %% switch1==0 | m==1){ # m==1 condition as have to calculate all liks on first step
      theta_star = SampleTheta(thetatab[m,], m,cov_matrix_theta,nparam=sum(cov_matrix_theta0)) #resample theta
      forcetab_star=forcetab
    }else{
      forcetab_star = t(apply(forcetab,1,function(x){SampleTheta(x, m,cov_matrix_force,nparam=sum(cov_matrix_force))})) #resample theta)
      theta_star = thetatab[m,]
    }
    
    dmatrix = diag(strains)*(1-theta_star[["error"]])+(1-diag(strains))*theta_star[["sigma"]] # default cross-reaction structure
    
    # - - - - - - - - - - - - - - - -
    # LIKELIHOOD function - Only calculate for updated history
    
    lik_val = LikelihoodTitres(titre.data,dmatrix,forcetab_star,inf.n,strains,n_sample)
    # - - - - - - - - - - - - - - - -
    # Metropolis Hastings step
    
    output_prob = ComputeProbability(sum(likelihoodtab[m,]),sum(lik_val)) 
    
    #if(is.na(output_prob)){print(c(m,sum(likelihoodtab[m,]),sum(lik_val),theta_star))} # Print likelihood (For DEBUG)}
    if(is.na(output_prob) & m==1){stop('check initial parameter values')}
    
    if(runif(1) < output_prob){
      thetatab[m+1,] = theta_star
      if(m %% switch1!=0){forcetab = forcetab_star} # Only change if resampled
      likelihoodtab[m+1] = lik_val
      if(m %% switch1==0){accepttabT=c(accepttabT,1)}
      if(m %% switch1!=0){accepttabH=c(accepttabH,1)}
      
    }else{
      thetatab[m+1,] = thetatab[m,]
      likelihoodtab[m+1] = likelihoodtab[m]
      if(m %% switch1==0){accepttabT=c(accepttabT,0)}
      if(m %% switch1!=0){accepttabH=c(accepttabH,0)}
    }
    
    
    if(m<max(100)){
      accept_rateT=0.234 # target acceptance rate for theta
      accept_rateH=0.234 # target acceptance rate for infection history
    }else{
      accept_rateT=sum(accepttabT)/length(accepttabT)
      accept_rateH=sum(accepttabH)/length(accepttabH)
    }
    
    if(m %% min(runs,10) ==0){ # Save force of infection every 10 runs
      forcetabCollect=rbind(forcetabCollect,forcetab)
    }
    
    if(m %% min(runs,1000) ==0){
      print(c(m,accept_rateT,varpart_prob0,likelihoodtab[m]))
      save(likelihoodtab,thetatab,inf_years,inf.n,strains,titre.data,forcetab,n_sample,forcetabCollect,switch1,file=paste("posterior_runs/outputR_f",seedi,".RData",sep=""))
    }
    
  } #End runs loop
  
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# RUN INFERENCE

inference_model<-function(seedK,strains,runsMCMC,scenario,per_sample,switch0=5){
  
  # Convert simulated samples into probabilities
  set.seed(seedK)
  
  load(paste("R_datasets/Sero_sim_",seedK,"_",strains,".RData",sep=""))
  
  # Impose ICs and run MCMC
  
  if(scenario==1){
    pmask0=c("sigma","error")
    force_constrain=matrix(1,nrow=strains,ncol=length(inf.years.sera))
  }
  
  if(scenario==2){
    pmask0=NULL
    force_constrain=matrix(1,nrow=strains,ncol=length(inf.years.sera))
  }
  
  if(scenario==3){
    pmask0=NULL
    force_constrain=matrix(1,nrow=strains,ncol=length(inf.years.sera))
    force_constrain[strains,2:length(inf.years.sera)]=0 # Constrain ZIKV to final year only
  }
  
  if(scenario==4){
    pmask0=NULL
    force_constrain=matrix(1,nrow=strains,ncol=length(inf.years.sera))
    force_constrain=(prob_inf0>0.01) # Constrain any year <1%
  }
  
  
  # TO ADD: FIT STRAINS
  run_mcmc(
    titre.data, # Note that this starts present and goes back into past (i.e. ordered with age). Data as proportions
    n_sample, # Total samples in each age group and each strain
    inf.years.sera,
    theta=c(error=0.05,sigma=0.2),
    strains,
    prob_inf=NULL, #prob_inf0, # initial conditions
    force_constrain, # matrix to constrain circulation years - Note that this is increasing age
    runs=runsMCMC,
    switch1=switch0, # How often sample FOI vs theta
    seedi=paste(per_sample,"_",scenario,"_",seedK,sep=""),
    pmask=pmask0 #NULL #
  )
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot dot distribution for Figure 2

plot.performance<-function(per_sample,age_out,strains,scenarioN=4,runs){
  
  # per_sample=per_sample0; age_out=20; strains=5; scenarioN=4; runs=seedRuns

  store.prop = array(NA,dim=c(runs,scenarioN))
  
  for(scenario in 1:scenarioN){
    
    if(scenario==0){
      # Leave out the base case for now...
      #for(seedK in 1:runs){
      #  load(paste("R_datasets/Sero_sim_",seedK,"_",strains,".RData",sep=""))
      #  store.prop[seedK,scenario+1] = as.numeric(prob_inf0[5,age_out])
      #}
      
    }else{
      
      # Estimate median probability for each inference run
      
      for(seedK in 1:runs){
        load(paste("R_datasets/Sero_sim_",seedK,"_",strains,".RData",sep=""))
        load(paste("posterior_runs/outputR_f",paste(per_sample,"_",scenario,"_",seedK,sep=""),".RData",sep=""))
        
        siml_force=t(time.series); 
        nblock=length(forcetabCollect)/(length(inf_years)*strains) # get blocks
        post_force=forcetabCollect[((nblock-1)*strains+1):(nblock*strains),]
        store.val=array(NA,dim=c(3,strains,inf.n))
        store.valCI=array(NA,dim=c(3,strains,inf.n))
        
        # Calculate median etc.
        for(ii in 1:strains){
          minN=round(0.2*(nblock-1))
          post_forceA=forcetabCollect[c(minN:(nblock-1))*strains+ii,]
          store.val[,ii,] = apply(post_forceA,2,function(x){c.nume(x)})
        }
        
        lik=NULL
        
        # Store 95% CI
        for(kk in 1:3){
          for(ii in 1:(inf.n)){ # Iterate across years
            p.inf=NULL
            for(jj in 1:strains){ # Iterate across strains
              p.inf=c(p.inf,1-exp(-sum(store.val[kk,jj,1:(ii)]))) # Infected in this period with strain j
            }
            # Convert to probability of infection
            store.valCI[kk, , ii] = p.inf
          }
        }
        
        store.prop[seedK,scenario] = store.valCI[1,5,age_out]
        
      } # end seed loop

    } # end scenario check
    
  } # end scenario loop
  
  
  par(mfrow=c(1,1))
  par(mar = c(3,3,1,1))
  par(mgp=c(1.8,0.6,0))
  
  plot_table=store.prop
  dim(plot_table)=NULL
  widX=0.15
  
  # Plot simulated data first
  
  # calculation median and IQR
  range=t(apply(store.prop,2,function(x){quantile(x,c(0.25,0.5,0.75))}))
  
  x_axis=rep((2*runif(runs)-1)*widX,scenarioN)+rep(1:scenarioN, each = runs) # Add scatter
  
  plot(x_axis,plot_table,xlim=c(0.5,4.5),ylim=c(0.4,1),cex=0.5,pch=19,ylab=paste("probability ",age_out,"yo previously infected with ZIKV",sep=""), xlab="model")
  
  # median
  for(ii in 1:4){
    lines(c(ii-widX*2,ii+widX*2),c(range[ii,1],range[ii,1]),col="blue",lwd=2,lty=2)
    lines(c(ii-widX*2,ii+widX*2),c(range[ii,2],range[ii,2]),col="blue",lwd=2)
    lines(c(ii-widX*2,ii+widX*2),c(range[ii,3],range[ii,3]),col="blue",lwd=2,lty=2)
  }
  
  lines(c(0,5),c(0.5,0.5),col="red",lwd=2)
  
  dev.copy(pdf,paste("plot_simulations/performance.pdf",sep=""),width=5,height=5)
  dev.off()
  
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Plot posteriors

plot.posteriors<-function(per_sample,strains,scenario,seedK){
  
  # per_sample = 10 ; strains = 5 ;  scenario = 1 ; seedK = 1
  
  load(paste("R_datasets/Sero_sim_",seedK,"_",strains,".RData",sep=""))
  load(paste("posterior_runs/outputR_f",paste(per_sample,"_",scenario,"_",seedK,sep=""),".RData",sep=""))
  
  par(mfrow=c(2,4))
  par(mar = c(3,3,1,1))
  par(mgp=c(1.8,0.6,0))
  
  colA=rgb(0.8,0.8,0.8)
  col.list=list(col2=rgb(0.9,0.6,0),col1=rgb(0.2,0,0.8),col3=rgb(0.1,0.6,0.2),col4=rgb(1,0.4,1),col5=rgb(0.8,0,0.2))
  alphc=0.2
  col.listF=list(col2=rgb(0.9,0.6,0,alphc),col1=rgb(0.2,0,0.8,alphc),col3=rgb(0.1,0.6,0.2,alphc),col4=rgb(1,0.4,1,alphc),col5=rgb(0.8,0,0.2,alphc))
  lw.1=1.5
  age.yr=c(0:(inf.n-1))
  label.age=seq(0,length(age.yr),2)
  
  # Plot posteriors
  maxlik=(likelihoodtab==max(likelihoodtab))
  run2=length(likelihoodtab)
  run1=round(0.2*run2)
  plot(likelihoodtab[run1:run2],type="l",ylab="log likelihood",xlab="iteration")
  theta_post=data.frame(thetatab[run1:run2,])
  hist(1-theta_post[["sigma"]],main="",col=colA,xlab="specificity",prob=TRUE,xlim=c(0.5,1)); abline(v=1-theta.serology[["sigma"]],col="red")
  hist(1-theta_post[["error"]],main="",col=colA,xlab="sensitivity",prob=TRUE,xlim=c(0.5,1)); abline(v=1-theta.serology[["error"]],col="red")
  
  # - - - 
  # Plot infection curves -  last sample from MCMC against observed data
  
  siml_force=t(time.series); 
  
  nblock=length(forcetabCollect)/(length(inf_years)*strains) # get blocks
  post_force=forcetabCollect[((nblock-1)*strains+1):(nblock*strains),]
  #post_force=1-exp(-post_force) # Convert FoI to probability
  #post_force=post_force[,c(length(inf_years):1)] # Swap arrangement to match
  
  store.val=array(NA,dim=c(3,strains,inf.n))
  store.valCI=array(NA,dim=c(3,strains,inf.n))
  
  # Calculate median etc.
  for(ii in 1:strains){
    post_forceA=forcetabCollect[c(0:(nblock-1))*strains+ii,]
    store.val[,ii,] = apply(post_forceA,2,function(x){c.nume(x)})
  }
  
  lik=NULL
  
  # Store 95% CI
  
  for(kk in 1:3){
    for(ii in 1:(inf.n)){ # Iterate across years
      p.inf=NULL
      for(jj in 1:strains){ # Iterate across strains
        p.inf=c(p.inf,1-exp(-sum(store.val[kk,jj,1:(ii)]))) # Infected in this period with strain j
      }
      # Convert to probability of infection
      store.valCI[kk, , ii] = p.inf
    }
  }
  
  #for(ii in 1:(inf.n)){ # Iterate across years
  #  p.inf=NULL
  #  for(jj in 1:strains){ # Iterate across strains
  #    p.inf=c(p.inf,1-exp(-sum(post_force[jj,1:(ii)]))) # Infected in this period with strain j
  #  }
  #  # Convert to probability of infection
  #  prob.positive = p.inf # 1-apply(1-t(dmatrix*p.inf),1,prod)
  #  lik=rbind(lik,prob.positive)
  #}
  
  
  plot(age.yr,f.y(historytabSim[,1]),ylim=c(0,1.01),col=col.list$col1,xlab="",ylab="proportion infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(f.y(store.valCI[1,1,]),col=col.list$col1,lty=2,lwd=lw.1)
  polygon(c(c(1:inf.n),rev(c(1:inf.n))),c(f.y(store.valCI[2,1,]),rev(f.y(store.valCI[3,1,]))),lty=0,col=col.listF$col1)
  #axis(1,at=label.age,labels=f.y(label.age))
  
  plot(age.yr,f.y(historytabSim[,2]),ylim=c(0,1.01),col=col.list$col2,xlab="",ylab="proportion infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(f.y(store.valCI[1,2,]),col=col.list$col2,lty=2,lwd=lw.1)
  polygon(c(c(1:inf.n),rev(c(1:inf.n))),c(f.y(store.valCI[2,2,]),rev(f.y(store.valCI[3,2,]))),lty=0,col=col.listF$col2)
  #axis(1,at=label.age,labels=f.y(label.age))
  
  plot(age.yr,f.y(historytabSim[,3]),ylim=c(0,1.01),col=col.list$col3,xlab="",ylab="proportion infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(f.y(store.valCI[1,3,]),col=col.list$col3,lty=2,lwd=lw.1)
  polygon(c(c(1:inf.n),rev(c(1:inf.n))),c(f.y(store.valCI[2,3,]),rev(f.y(store.valCI[3,3,]))),lty=0,col=col.listF$col3)
  #axis(1,at=label.age,labels=f.y(label.age))
  
  plot(age.yr,f.y(historytabSim[,4]),ylim=c(0,1.01),col=col.list$col4,xlab="",ylab="proportion infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(f.y(store.valCI[1,4,]),col=col.list$col4,lty=2,lwd=lw.1)
  polygon(c(c(1:inf.n),rev(c(1:inf.n))),c(f.y(store.valCI[2,4,]),rev(f.y(store.valCI[3,4,]))),lty=0,col=col.listF$col4)
  #axis(1,at=label.age,labels=f.y(label.age))
  
  plot(age.yr,f.y(historytabSim[,5]),ylim=c(0,1.01),col=col.list$col5,xlab="",ylab="proportion infected",type="l",bty="n",xaxs="i",yaxs="i",lwd=lw.1)
  lines(f.y(store.valCI[1,5,]),col=col.list$col5,lty=2,lwd=lw.1)
  polygon(c(c(1:inf.n),rev(c(1:inf.n))),c(f.y(store.valCI[2,5,]),rev(f.y(store.valCI[3,5,]))),lty=0,col=col.listF$col5)
  #axis(1,at=label.age,labels=f.y(label.age))
  dev.copy(pdf,paste("plot_simulations/posterior_plot",paste(per_sample0,"_",scenario,"_",seedK,sep=""),".pdf",sep=""),width=10,height=6)
  dev.off()
  
  post.tab=cbind(  
    c(c.text(1-theta_post[["error"]]),
      c.text(1-theta_post[["sigma"]]),
      c.text(x=c(mean(store.valCI[1,5,]),mean(store.valCI[2,5,]),mean(store.valCI[3,5,]))), # ZIKV attack rate
      c.text(x=c(mean(store.valCI[1,5,21]),mean(store.valCI[2,5,21]),mean(store.valCI[3,5,21]))) # ZIKV attack rate
  ))

  
  write.csv(post.tab,paste("plot_simulations/parameter_estimates",scenario,"_seed",seedK,".csv",sep=""))
  

}

