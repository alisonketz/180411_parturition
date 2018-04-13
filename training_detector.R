training = function(d,eps,pw=126){
    
    #d is dataframe of features, 
    #the first column is the id
    #the second column is the julian day of observation
    
    nCovs = dim(d)[2]-3
    id = unique(d[,1])
    nInd = length(id)
    ci = 4:dim(d)[2]

    if(length(eps)!=nCovs){cat("length(epsilon) != number covariates, try again \n");return}
    
    #For individual case, run anomaly dection on single individual
    if(nInd==1){
        d_temp=d[d[,2]>pw,]
        n_temp=dim(d_temp)[1]
        detect_quant = quantile(d[d[,2]<pw,ci],eps,na.rm=TRUE)
        if(nCovs==1){#nInd = 1, ncovs=1
            ind_mean = mean(d[d[,2]<pw,ci],na.rm=TRUE)
            ind_sd =  sd(d[d[,2]<pw,ci],na.rm=TRUE)
            
            threshold_density = dnorm(detect_quant,ind_mean,ind_sd)
            threshold_p = pnorm(detect_quant,mean=ind_mean,sd=ind_sd)
            lower_prob=pnorm(d_temp[,ci],mean=ind_mean,sd=ind_sd)
            
            detect_density = dnorm(d_temp[,ci],ind_mean,ind_sd)
            
            hits_indx = which(is.finite(detect_density) & detect_density<=threshold_density & lower_prob <= threshold_p)
            outs_indx = which(!(1:n_temp %in% hits_indx))
            alarm = list(rbind(d_temp[hits_indx,]))
            
            results_prob=rep(NA,n_temp)
            results_prob[hits_indx]=1-(lower_prob[hits_indx]/threshold_p)
            results_prob[outs_indx] =1-(1-lower_prob[outs_indx])/(1-threshold_p)

        }#end nInd=1,nCovs=1
        else {# if nInd=1 nCovs >1
                detect_quant=rep(NA,nCovs)
                for(i in 1:nCovs){
                  detect_quant[i] = quantile(d[d[,2]<pw,i+2],probs=eps[i],na.rm=TRUE)
                }
                ind_mean = apply(d[d[,2]<pw,ci],2,mean,na.rm=TRUE)
                ind_sd = apply(d[d[,2]<pw,ci],2,sd,na.rm=TRUE)
                Sigma=diag(ind_sd)
                threshold_density = dmvnorm(detect_quant,ind_mean,Sigma)
                threshold_p = pmvnorm(lower=rep(-Inf,nCovs),upper=detect_quant,mean=ind_mean,sigma=Sigma)[1]

                detect_density = rep(NA,n_temp)
                lower_prob=rep(NA,n_temp)
                for(i in 1:n_temp){
                    up=as.numeric(d_temp[i,ci])
                    up.na=is.na(up)
                    detect_density[i] = dmvnorm(up[!up.na],ind_mean[!up.na],diag(ind_sd[!up.na]))
                    lower_prob[i]=pmvnorm(lower=rep(-Inf,nCovs-sum(up.na)),upper = up[!up.na],mean=ind_mean[!up.na],sigma=diag(ind_sd[!up.na]))[1]
                }
                
                hits_indx = which(is.finite(detect_density) & detect_density<=threshold_density & lower_prob <= threshold_p)
                outs_indx = which(!(1:n_temp %in% hits_indx))
                alarm=list(rbind(d_temp[hits_indx,]))
                
                results_prob=rep(NA,n_temp)
                results_prob[hits_indx]=1-(lower_prob[hits_indx]/threshold_p)
                results_prob[outs_indx] =1-(1-lower_prob[outs_indx])/(1-threshold_p)
        }#end else nCovs>1,nInd = 1
    }#end nInd =1
    else{ # nInd>1       
        n_temp=rep(NA,nInd)
        for(j in 1:nInd){
            d_temp1=d[d[,1]==id[j],]
            d_temp=d_temp1[d_temp1[,2]>pw,]
            n_temp[j]=dim(d_temp)[1]
        }
        n_temp_max=max(n_temp)
        
        alarm=rep(list(),nInd)
        hits_indx=rep(list(),nInd)
        outs_indx=rep(list(),nInd)
        threshold_p=rep(NA,nInd)
        threshold_density=rep(NA,nInd)
        detect_density=matrix(NA,nr=n_temp_max,nc=nInd)
        lower_prob=matrix(NA,nr=n_temp_max,nc=nInd)
        results_prob=matrix(NA,nr=n_temp_max,nc=nInd)
        
        if(nCovs==1){ #nInd>1, nCovs = 1
            detect_quant=rep(NA,nInd)
            ind_mean=rep(NA,nInd)
            ind_sd=rep(NA,nInd)
            for(j in 1:nInd){
                d_temp1=d[d[,1]==id[j],]

                ind_mean[j] = mean(d_temp1[d_temp1[,2]<=pw,ci],na.rm=TRUE)
                ind_sd[j] = sd(d_temp1[d_temp1[,2]<=pw,ci],na.rm=TRUE)
                
                detect_quant[j]=quantile(as.numeric(d_temp1[d_temp1[,2]<=pw,ci]),eps,na.rm=TRUE)
                threshold_density[j] = dnorm(detect_quant[j],ind_mean[j],ind_sd[j])
                threshold_p[j] = pnorm(detect_quant[j],mean=ind_mean[j],sd=ind_sd[j])
                
                d_temp=d_temp1[d_temp1[,2]>pw,]
                
                #probability of hit/nothit
                for(i in 1:n_temp[j]){
                    lower_prob[i,j]=pnorm(as.numeric(d_temp[i,ci]),mean=ind_mean[j],ind_sd[j])
                    detect_density[i,j] = dnorm(as.numeric(d_temp[i,ci]),ind_mean[j],ind_sd[j])
                }
                
                
                hits_indx[[j]] = which(is.finite(detect_density[,j]) & detect_density[,j]<=threshold_density[j] & lower_prob[,j] <= threshold_p[j])
                outs_indx[[j]] = which(!(1:n_temp[j] %in% hits_indx[[j]]))
                
                alarm[[j]]=d_temp[hits_indx[[j]],]

                results_prob[hits_indx[[j]],j] = 1-(lower_prob[hits_indx[[j]],j]/threshold_p[j])
                results_prob[outs_indx[[j]],j] =1-(1-lower_prob[outs_indx[[j]],j])/(1-threshold_p[j])
                
            }#endfor nInd
        }#endif nCovs==1
        else{#nind>1,nCovs >1

            detect_quant=matrix(NA,nr=nCovs,nc=nInd)
            
            ind_mean=matrix(NA,nr=nCovs,nc=nInd)
            ind_sd=matrix(NA,nr=nCovs,nc=nInd)

            for(j in 1:nInd){
                d_temp1=d[d[,1]==id[j],]
                ind_mean[,j] = apply(d_temp1[d_temp1[,2]<pw,ci],2,mean,na.rm=TRUE)
                ind_sd[,j] = apply(d_temp1[d_temp1[,2]<pw,ci],2,sd,na.rm=TRUE)

                detect_quant[,j]=quantile(d_temp1[d_temp1[,2]<pw,ci],eps,na.rm=TRUE)
                threshold_density[j] = dmvnorm(detect_quant[,j],ind_mean[,j],diag(ind_sd[,j]))
                threshold_p[j] = pmvnorm(lower=rep(-Inf,nCovs),upper=detect_quant[,j],mean=ind_mean[,j],sigma=diag(ind_sd[,j]))[1]
                
                d_temp=d_temp1[d_temp1[,2]>pw,]
                
                for(i in 1:n_temp[j]){
                    up=as.numeric(d_temp[i,ci])
                    up.na=is.na(up)
                    detect_density[i,j] = dmvnorm(up[!up.na],mean=ind_mean[!up.na,j],sigma=diag(ind_sd[!up.na,j]))
                    lower_prob[i,j]=pmvnorm(lower=rep(-Inf,nCovs-sum(up.na)),upper = up[!up.na],mean=ind_mean[!up.na,j],sigma=diag(ind_sd[!up.na,j]))[1]
                }
                
                hits_indx[[j]] = which(is.finite(detect_density[,j]) & detect_density[,j]<=threshold_density[j] & lower_prob[,j] <= threshold_p[j])
                outs_indx[[j]] = which(!(1:n_temp[j] %in% hits_indx[[j]]))
                
                alarm[[j]]=d_temp[hits_indx[[j]],]
                
                results_prob[hits_indx[[j]],j] = 1-(lower_prob[hits_indx[[j]],j]/threshold_p[j])
                results_prob[outs_indx[[j]],j] =1-(1-lower_prob[outs_indx[[j]],j])/(1-threshold_p[j])
                
                }#endfor
            }#end else
        }#end group nInd>1 else

        return(list(alarm=alarm,
                    detect_density=detect_density,
                    threshold_density=threshold_density,
                    detect_quant=detect_quant,
                    results_prob=results_prob,
                    hits_indx=hits_indx,
                    outs_indx=outs_indx,
                    lower_prob=lower_prob))
    
}
