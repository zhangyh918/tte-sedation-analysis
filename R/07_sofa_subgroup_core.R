# Modified-SOFA tertile subgroup analysis.
# The main IPCW models are unchanged, and SOFA values are not imputed.
ccw_sofa_estimate <- function(w,sofa,draw=NULL) {
  r<-merge(copy(w$risk[target=='ext']),sofa[,.(stay_id,sofa_group)],by='stay_id')
  r<-r[!is.na(sofa_group)]
  if(!is.null(draw)) {
    r$frequency<-NULL;r<-merge(r,draw,by='stay_id',allow.cartesian=TRUE)
    r[,cluster_id:=bootstrap_id]
  } else r[,cluster_id:=stay_id]
  r[,`:=`(strategy=factor(strategy,levels=c('SLA','DTLS')),
    sofa_group=factor(sofa_group,levels=c('T1','T2','T3')))]
  v<-merge(copy(w$terminal),sofa[,.(stay_id,sofa_group)],by='stay_id')[!is.na(sofa_group)]
  warns<-character()
  cfit<-function(dat,interaction=FALSE) {
    form<-if(interaction)Surv(start/24,stop/24,extubation)~strategy*sofa_group else
      Surv(start/24,stop/24,extubation)~strategy
    withCallingHandlers(if(is.null(draw))coxph(form,data=dat,weights=weight,
      cluster=cluster_id,robust=TRUE,ties='efron') else
      coxph(form,data=dat,weights=weight,robust=FALSE,ties='efron'),
      warning=function(z){warns<<-c(warns,conditionMessage(z));invokeRestart('muffleWarning')})
  }
  rows<-list()
  for(g in c('T1','T2','T3')) {
    fit<-cfit(r[r$sofa_group==g]);b<-unname(coef(fit)[1])
    if(!is.finite(b))stop('Nonfinite SOFA subgroup coefficient')
    m<-v[v$sofa_group==g,.(vfd=weighted.mean(observed_vfd,weight*frequency)),by=strategy]
    if(nrow(m)!=2)stop('Missing SOFA subgroup strategy')
    row<-data.table(group=g,log_hr=b,hr=exp(b),
      vfd_difference=m$vfd[m$strategy=='DTLS']-m$vfd[m$strategy=='SLA'])
    if(is.null(draw)) {
      se<-sqrt(vcov(fit)[1,1]);row[,`:=`(fixed_nuisance_lower=exp(b-qnorm(.975)*se),
        fixed_nuisance_upper=exp(b+qnorm(.975)*se))]
    }
    rows[[g]]<-row
  }
  fi<-cfit(r,TRUE);ii<-grep(':',names(coef(fi)))
  stopifnot(length(ii)==2L)
  co<-coef(fi)[ii];if(any(!is.finite(co)))stop('Nonfinite SOFA interaction')
  inter<-data.table(T2=unname(co[1]),T3=unname(co[2]))
  if(is.null(draw)) {
    stat<-drop(t(co)%*%solve(vcov(fi)[ii,ii,drop=FALSE],co))
    inter[,`:=`(fixed_nuisance_wald=stat,df=2L,fixed_nuisance_p=pchisq(stat,2,lower.tail=FALSE))]
  }
  list(groups=rbindlist(rows),interaction=inter,warnings=paste(unique(warns),collapse=' | '))
}
