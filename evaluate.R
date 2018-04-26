
evaluate= function(alarm,possible.hits,nInd,vitdropday){
    if(nInd==1){

        d=alarm[[1]][alarm[[1]][,2]<=vitdropday,]

        # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
        tp = sum(d[,3],na.rm=TRUE)
        
        # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
        fp = sum(d[,3]==0,na.rm=TRUE)
        
        # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
        fn = possible.hits - tp
        
        tp.pop = tp
        fp.pop = fp
        fn.pop = fn
        
        # precision and recall
        out.prec = ifelse((tp.pop+fp.pop)==0,0,tp.pop/(tp.pop+fp.pop))
        out.recall = ifelse((tp.pop+fn.pop)==0,0,tp.pop/(tp.pop+fn.pop))
        
        # F1 value
        out.F1 =ifelse((out.prec+out.recall) == 0,0,(2*out.prec*out.recall)/(out.prec+out.recall))
    }
    else{
        d=rep(list(),nInd)
        
        for(j in 1:nInd){
            d[[j]]=alarm[[j]][alarm[[j]][,2]<=vitdropday[j],]
        }
    
        tp = rep(NA,nInd)
        fp = rep(NA,nInd)
        fn = rep(NA,nInd)
        
        for(j in 1:nInd){

            # TRUE POSITIVE:  groundtruth data  says it's an anomaly and so algorithm does.
            tp[j] = sum(d[[j]][,3],na.rm=TRUE)
            
            # FALSE POSITIVE:  groundtruth data says it's not an anomaly, but algorithm says anomaly.
            fp[j] = sum(d[[j]][,3]==0,na.rm=TRUE)
            
            # FALSE NEGATIVE: groundtruth data says it's an anomaly, but algorithm says not anomaly.
            fn[j] = possible.hits[j] - tp[j]
            
        }
        
        tp.pop = sum(tp)
        fp.pop = sum(fp)
        fn.pop = sum(fn)
        
        # precision and recall
        # F1 value
        out.prec =ifelse(tp.pop+fp.pop == 0,0,tp.pop/(tp.pop+fp.pop))
        out.recall =ifelse(tp.pop+fp.pop == 0,0,tp.pop/(tp.pop+fn.pop))
        out.F1 =ifelse((out.prec+out.recall) == 0,0,(2*out.prec*out.recall)/(out.prec+out.recall))

    }#endelse
    
    return(list(tp=tp,fp=fp,fn=fn,out.prec=out.prec,out.recall=out.recall,out.F1=out.F1,tp.pop=tp.pop,fp.pop=fp.pop,fn.pop=fn.pop))
    
}
