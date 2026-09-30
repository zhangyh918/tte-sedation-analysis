# Unstabilized inverse-probability-of-censoring weights paired with marginal
# means and a strategy-only cause-specific Cox model.
suppressPackageStartupMessages({library(data.table);library(survival)})
setDTthreads(1)
ccw_basevars <- c('age','sex','admission_type','diabetes','ckd','liver_disease','chf','copd',
  'pf_ratio','fio2_value','peep_value','map_value','heart_rate','spo2_value','vasopressor_at_t0',
  'initial_sedative','lactate','creatinine','bilirubin','platelets')
ccw_lagvars <- c('lag_rass','lag_propofol','lag_midazolam','lag_dex')

# Exact type-7 quantiles of literal integer-frequency replication, without
# materializing repeated nuisance rows. Tested against stats::quantile(rep()).
frequency_quantile <- function(x,frequency,probs=c(.01,.99)) {
  stopifnot(length(x)==length(frequency),all(is.finite(x)),
    all(is.finite(frequency)),all(frequency>=0),all(frequency==floor(frequency)))
  keep<-frequency>0;x<-x[keep];frequency<-frequency[keep]
  stopifnot(length(x)>0,all(probs>=0 & probs<=1))
  o<-order(x);x<-x[o];f<-frequency[o];cs<-cumsum(f);n<-sum(f)
  at<-function(k)x[findInterval(k-1,cs)+1L]
  h<-(n-1)*probs+1;j<-floor(h);g<-h-j
  at(j)+(at(pmin(j+1,n))-at(j))*g
}

ccw_prepare_design <- function(d,predictors) {
  df<-as.data.frame(d);f<-df$frequency;recipe<-list();missing<-list()
  stopifnot(all(f>0))
  for(v in predictors) {
    z<-df[[v]];if(is.logical(z))z<-as.integer(z)
    if(is.numeric(z)) {
      bad<-!is.finite(z)
      med<-if(any(!bad))frequency_quantile(z[!bad],f[!bad],.5) else 0
      z[bad]<-med;mu<-sum(z*f)/sum(f)
      ss<-if(sum(f)>1)sqrt(sum(f*(z-mu)^2)/(sum(f)-1)) else 1
      if(!is.finite(ss)||ss==0)ss<-1
      df[[v]]<-(z-mu)/ss;recipe[[v]]<-list(kind='numeric',replace=med,mean=mu,sd=ss)
    } else {
      z<-as.character(z);bad<-is.na(z)|z==''
      observed<-data.table(value=z[!bad],f=f[!bad])[,.(f=sum(f)),by=value]
      setorder(observed,-f,value)
      mode<-if(nrow(observed))observed$value[1] else 'Missing'
      z[bad]<-mode;df[[v]]<-factor(z)
      recipe[[v]]<-list(kind='factor',replace=mode,levels=levels(df[[v]]))
    }
    missing[[v]]<-data.table(variable=v,imputed_rows=sum(f[bad]),all_missing=all(bad),
      total_rows=sum(f))
  }
  df$node_grace<-factor(paste(df$node_id,df$counter_before,sep='_'))
  list(data=df,recipe=recipe,missing=rbindlist(missing))
}

ccw_fit_denominator <- function(d,rass_only=FALSE) {
  setorder(d,stay_id,node_id)
  stopifnot(all(d$lag_information_end_h<=d$node_start_h+1e-7),
    all(d$assessment_h>d$node_start_h),!anyDuplicated(d[,.(stay_id,node_id)]))
  forced<-rass_only & d$counter_before==0L
  stopifnot(!any(forced & d$uncensored==0L))
  prep<-ccw_prepare_design(d,c(ccw_basevars,ccw_lagvars))
  df<-droplevels(prep$data[!forced,,drop=FALSE])
  if(nrow(df)==0 || length(unique(df$uncensored))<2)
    stop('Non-structural retention model lacks both outcome classes')
  candidate<-c('node_grace',ccw_basevars,ccw_lagvars)
  use<-candidate[vapply(df[,candidate,drop=FALSE],function(x)length(unique(x))>1,logical(1))]
  X<-model.matrix(reformulate(use),df);warnings<-character()
  fit<-withCallingHandlers(glm.fit(X,df$uncensored,weights=df$frequency,
    family=binomial(),control=glm.control(maxit=100,epsilon=1e-10)),
    warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
  if(!isTRUE(fit$converged))stop('Denominator model failed to converge')
  if(fit$rank!=ncol(X))stop('Denominator model is rank deficient')
  p<-rep(1,nrow(d));p[!forced]<-fit$fitted.values
  if(any(!is.finite(p)|p<=0|p>1))stop('Invalid conditional retention probability; no silent flooring')
  q<-d[,.(stay_id,node_id,assessment_h,node_start_h,lag_information_end_h,counter_before,
    uncensored,frequency,strategy,target)]
  q[,`:=`(p_den=p,structural_one=forced)]
  q[,raw_weight:=exp(cumsum(-log(p_den))),by=stay_id]
  if(any(!is.finite(q$raw_weight)))stop('Cumulative weight overflow; review support')
  diagnostics<-data.table(rows=sum(d$frequency),unique_rows=nrow(d),converged=fit$converged,
    rank=fit$rank,columns=ncol(X),iterations=fit$iter,structural_one=sum(d$frequency[forced]),
    p_below_1e6=sum(d$frequency[p<1e-6]),p_below_1e6_retained=sum(d$frequency[p<1e-6 & d$uncensored==1]),
    p_above_one_minus_1e6=sum(d$frequency[p>1-1e-6 & !forced]),
    min_probability=min(p),max_raw_all=max(q$raw_weight),
    max_raw_retained=max(q$raw_weight[q$uncensored==1]),
    weighted_brier=weighted.mean((d$uncensored-p)^2,d$frequency),
    warnings=paste(unique(warnings),collapse=' | '),formula=paste('uncensored ~',paste(use,collapse=' + ')))
  list(probabilities=q,diagnostics=diagnostics,recipe=prep$recipe,missing=prep$missing,
    coefficients=fit$coefficients,columns=colnames(X))
}

ccw_weight_data <- function(h,counts=NULL) {
  b<-copy(h$baseline)
  stopifnot(!anyDuplicated(b$stay_id),!anyDuplicated(b$subject_id),length(ccw_basevars)==20L,
    h$partial_policy=='fraction50')
  if(is.null(counts))counts<-data.table(stay_id=b$stay_id,frequency=1L)
  stopifnot(!anyDuplicated(counts$stay_id),all(counts$frequency>=0),
    all(counts$frequency==floor(counts$frequency)),all(counts$stay_id %in% b$stay_id))
  counts<-counts[frequency>0];b<-merge(b,counts,by='stay_id')
  d<-merge(h$decisions,b[,c('stay_id','frequency',ccw_basevars),with=FALSE],by='stay_id')
  prob<-list();dg<-list();recipes<-list();missing<-list();k<-0L
  for(target_name in c('vfd','ext'))for(strategy_name in c('SLA','DTLS')) {
    take<-d$target==target_name & d$strategy==strategy_name
    fit<-ccw_fit_denominator(copy(d[take]),isTRUE(h$scenario$rass_only))
    k<-k+1L;key<-paste(target_name,strategy_name,sep='_')
    prob[[k]]<-fit$probabilities
    fit$diagnostics[,`:=`(target=target_name,strategy=strategy_name)];dg[[k]]<-fit$diagnostics
    recipes[[key]]<-list(recipe=fit$recipe,coefficients=fit$coefficients,columns=fit$columns)
    fit$missing[,`:=`(target=target_name,strategy=strategy_name)];missing[[k]]<-fit$missing
  }
  p<-rbindlist(prob)
  limits<-p[,{
    z<-frequency_quantile(raw_weight,frequency,c(.01,.99))
    list(lower=z[1],upper=z[2],model_rows=sum(frequency))
  },by=target]
  p<-merge(p,limits[,.(target,lower,upper)],by='target')
  p[,weight:=pmax(lower,pmin(upper,raw_weight))]
  retained<-p[uncensored==1L,.(stay_id,strategy,target,time=assessment_h,weight,raw_weight)]
  stopifnot(!anyDuplicated(retained[,.(stay_id,strategy,target,time)]))
  # The already-audited risk split starts at 0 or at a retained update.
  risk<-merge(copy(h$risk),counts,by='stay_id')
  risk<-merge(risk,retained,by.x=c('stay_id','strategy','target','start'),
    by.y=c('stay_id','strategy','target','time'),all.x=TRUE,sort=FALSE)
  if(any(is.na(risk$weight)&risk$start!=0))stop('Missing weight at a nonzero risk boundary')
  risk[start==0,`:=`(weight=1,raw_weight=1)]
  term<-merge(h$status[target=='vfd' & endpoint_observed],counts,by='stay_id')
  term<-merge(term,retained,by.x=c('stay_id','strategy','target','last_retained_decision_h'),
    by.y=c('stay_id','strategy','target','time'),all.x=TRUE,sort=FALSE)
  if(any(is.na(term$weight)&term$last_retained_decision_h!=0))stop('Missing terminal weight')
  term[last_retained_decision_h==0,`:=`(weight=1,raw_weight=1)]
  stopifnot(all(is.finite(risk$weight)),all(risk$weight>0),all(risk$start<risk$stop),
    all(term$last_retained_decision_h<term$endpoint_h),
    all(term$observed_vfd>=0 & term$observed_vfd<=28),
    all(term$observed_vfd[term$endpoint_type=='death']==0))
  # Any trigger is the last decision; it never becomes a carried risk weight.
  check<-p[order(stay_id,node_id),.(triggers=sum(uncensored==0),
    trigger_last=all(which(uncensored==0)==.N)),by=.(stay_id,strategy,target)]
  stopifnot(all(check$triggers<=1L),all(check$trigger_last))
  cal<-p[,.(rows=sum(frequency),observed=weighted.mean(uncensored,frequency),
    predicted=weighted.mean(p_den,frequency),retained=sum(frequency[uncensored==1]),
    structural_one=sum(frequency[structural_one])),by=.(target,strategy,node_id,counter_before)]
  terminal_diag<-term[,{
    w<-weight;f<-frequency;q<-frequency_quantile(w,f,c(.5,.99))
    list(observed_clones=sum(f),weighted_total=sum(w*f),kish_ess=sum(w*f)^2/sum(w^2*f),
      median_weight=q[1],p99_weight=q[2],max_weight=max(w),max_raw_weight=max(raw_weight),
      truncated_clones=sum(f[w!=raw_weight]),deaths=sum(f[endpoint_type=='death']),
      no_update_clones=sum(f[last_retained_decision_h==0]))
  },by=strategy]
  event_diag<-risk[target=='ext' & extubation==1,.(events=sum(frequency),
    event_weight_ess=sum(weight*frequency)^2/sum(weight^2*frequency),
    max_event_weight=max(weight)),by=strategy]
  support<-p[,.(decisions=sum(frequency),retained=sum(frequency[uncensored==1]),
    censored=sum(frequency[uncensored==0]),min_p_retained=if(any(uncensored==1))min(p_den[uncensored==1]) else NA_real_,
    max_raw_retained=if(any(uncensored==1))max(raw_weight[uncensored==1]) else NA_real_,
    max_weight_retained=if(any(uncensored==1))max(weight[uncensored==1]) else NA_real_),
    by=.(target,strategy,node_id)]
  list(risk=risk,terminal=term,probabilities=p,limits=limits,diagnostics=rbindlist(dg),
    recipes=recipes,missing=rbindlist(missing),calibration=cal,support=support,
    terminal_diagnostics=terminal_diag,event_diagnostics=event_diag,
    patients=sum(counts$frequency),numerator='one',scenario=h$scenario$scenario,
    clinical_effect_ready=FALSE)
}

ccw_estimate <- function(w,draw=NULL,check_ph=FALSE) {
  means<-w$terminal[,.(mean_vfd=weighted.mean(observed_vfd,weight*frequency)),by=strategy]
  stopifnot(setequal(means$strategy,c('SLA','DTLS')))
  r<-copy(w$risk[target=='ext'])
  if(!is.null(draw)) {
    # Literal replication preserves Efron ties; multiplying case weights by
    # bootstrap frequency does not in general reproduce that likelihood.
    r$frequency<-NULL
    r<-merge(r,draw,by='stay_id',allow.cartesian=TRUE)
    r[,cluster_id:=bootstrap_id]
  } else {
    stopifnot(all(r$frequency==1));r[,cluster_id:=stay_id]
  }
  r[,strategy:=factor(strategy,levels=c('SLA','DTLS'))]
  warnings<-character()
  fit<-withCallingHandlers(if(is.null(draw))coxph(Surv(start/24,stop/24,extubation)~strategy,
    data=r,weights=weight,cluster=cluster_id,robust=TRUE,ties='efron',
    x=check_ph,model=check_ph,control=coxph.control(iter.max=50)) else
    coxph(Surv(start/24,stop/24,extubation)~strategy,data=r,weights=weight,
      robust=FALSE,ties='efron',control=coxph.control(iter.max=50)),
    warning=function(z){warnings<<-c(warnings,conditionMessage(z));invokeRestart('muffleWarning')})
  beta<-unname(coef(fit)[1]);if(!is.finite(beta))stop('Nonfinite Cox estimate')
  ans<-data.table(scenario=w$scenario,patients=w$patients,
    mean_sla=means$mean_vfd[means$strategy=='SLA'],mean_dtls=means$mean_vfd[means$strategy=='DTLS'],
    vfd_difference=means$mean_vfd[means$strategy=='DTLS']-means$mean_vfd[means$strategy=='SLA'],
    log_hr=beta,hr=exp(beta),cox_warnings=paste(unique(warnings),collapse=' | '))
  ph<-NULL
  if(is.null(draw)) {
    se<-sqrt(vcov(fit)[1,1]);ans[,`:=`(cox_robust_se=se,
      cox_fixed_nuisance_lower=exp(beta-qnorm(.975)*se),
      cox_fixed_nuisance_upper=exp(beta+qnorm(.975)*se),
      cox_fixed_nuisance_p=2*pnorm(-abs(beta/se)))]
    if(check_ph)ph<-tryCatch(as.data.table(cox.zph(fit,transform='km')$table,keep.rownames='term'),
      error=function(e)data.table(term='ERROR',message=conditionMessage(e)))
  }
  list(result=ans,ph=ph,fit=if(is.null(draw))fit else NULL)
}
