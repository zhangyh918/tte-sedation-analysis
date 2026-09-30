# Sensitivity analyses for alternative infusion-coverage fractions.
repo_root <- Sys.getenv("TTE_REPO_ROOT", unset = normalizePath(".", mustWork = TRUE))
source(file.path(repo_root, "R", "08_sensitivity_core.R"))
source(file.path(repo_root, "R", "06_ipcw_estimation_core.R"))

coverage_policy <- function(fraction){
  stopifnot(length(fraction)==1L,is.finite(fraction),fraction %in% c(.25,.5,.75))
  paste0('fraction',as.integer(100*fraction))
}
coverage_nodes <- function(nodes,required_fraction){
  policy<-coverage_policy(required_fraction);d<-data.table::copy(nodes)
  stopifnot(all(is.finite(d$mv_h)),all(d$mv_h>=0 & d$mv_h<=12+1e-7),
    all(is.finite(d$coverage_h)),all(d$coverage_h>=0 & d$coverage_h<=d$mv_h+1e-6))
  qualifies<-d$mv_h>0 & d$coverage_h>=required_fraction*d$mv_h-1e-6
  d[,`:=`(coverage_required_fraction=required_fraction,
    state_fraction=data.table::fcase(mv_h==0,'paused',is.na(rass_carried),'missing',
      rass_carried>=-2,'light',qualifies,'deep',default='low_without'))]
  d
}

coverage_pause_resume_trace <- function(nodes,strategy=c('SLA','DTLS'),transition_h=48,
  horizon_h=144,endpoint_h=Inf,icu_out_h=NA_real_,partial_policy='fraction50',rass_only=FALSE){
  strategy<-match.arg(strategy);partial_policy<-match.arg(partial_policy,c('fraction25','fraction50','fraction75'))
  d<-as.data.frame(nodes);d<-d[order(d$node_id),,drop=FALSE]
  stopifnot(all(c('state_fraction','coverage_required_fraction','assessment_h','mv_h') %in% names(d)),
    all(d$coverage_required_fraction==as.numeric(sub('fraction','',partial_policy))/100))
  counter<-0L;censor_h<-Inf;records<-list();reason<-'fixed_horizon'
  for(i in seq_len(nrow(d))){
    a<-d[i,]
    if(a$node_start_h>=horizon_h)break
    if(!is.na(icu_out_h) && a$node_start_h>=icu_out_h){reason<-'index_icu_ended';break}
    t<-min(a$assessment_h,horizon_h)
    if(endpoint_h<=t){reason<-'terminal_event_first';break}
    if(a$mv_h==0)state<-'paused' else if(rass_only){
      state<-if(is.na(a$rass_carried))'missing' else if(a$rass_carried>=-2)'light' else 'deep'
    }else state<-a$state_fraction
    expected<-if(strategy=='DTLS' && a$node_start_h<transition_h)'deep' else 'light'
    before<-counter
    if(!state %in% c('paused','missing'))counter<-if(state==expected)0L else counter+1L
    trigger<-state=='low_without' || counter>1L
    records[[length(records)+1L]]<-data.frame(node_id=a$node_id,assessment_h=t,
      state=state,expected=expected,counter_before=before,counter_after=counter,
      censor_now=trigger,paused=state=='paused',partial_mv_node=a$partial_mv_node)
    if(trigger){censor_h<-t;reason<-'artificial_censor';break}
  }
  empty<-data.frame(node_id=integer(),assessment_h=numeric(),state=character(),expected=character(),
    counter_before=integer(),counter_after=integer(),censor_now=logical(),paused=logical(),partial_mv_node=logical())
  list(rows=if(length(records))do.call(rbind,records) else empty,censor_h=censor_h,
    stop_reason=reason,partial_policy=partial_policy,clinical_effect_ready=FALSE)
}

# Exact expression substitutions, with match counts. No estimator is retyped;
# the weight function changes only its permitted policy-label guard.
coverage_replace <- function(expr,needle,replacement,expected=1L){
  hits<-0L
  walk<-function(x){
    if(identical(x,needle)){hits<<-hits+1L;return(replacement)}
    if(is.call(x)||is.pairlist(x)||is.expression(x))for(i in seq_along(x))x[i]<-list(walk(x[[i]]))
    x
  }
  result<-walk(expr);stopifnot(hits==expected);result
}
coverage_horizon_history<-sync_horizon_history
formals(coverage_horizon_history)<-c(formals(sync_horizon_history),alist(required_fraction=.5))
body(coverage_horizon_history)<-coverage_replace(body(sync_horizon_history),as.name('pause_resume_trace'),as.name('coverage_pause_resume_trace'))
body(coverage_horizon_history)<-coverage_replace(body(coverage_horizon_history),'fraction50',quote(coverage_policy(required_fraction)))

coverage_build_histories<-function(b,nodes,horizon,scenario,required_fraction){
  stopifnot(all(nodes$coverage_required_fraction==required_fraction))
  nn<-split(as.data.frame(nodes),nodes$stay_id)
  all<-lapply(seq_len(nrow(b)),function(i){
    rr<-list();j<-0L
    for(s in c('SLA','DTLS'))for(t in c('vfd','ext')){
      j<-j+1L
      z<-coverage_horizon_history(nn[[as.character(b$stay_id[i])]],strategy=s,target=t,
        death_h=b$death_h[i],ext_h=b$ext_h[i],vfd=b$vfd[i],icu_out_h=b$icu_out_h[i],
        horizon_h=horizon,required_fraction=required_fraction)
      for(k in c('trace','decisions','status','risk')){
        z[[k]]$stay_id<-rep(b$stay_id[i],nrow(z[[k]]))
        z[[k]]$strategy<-rep(s,nrow(z[[k]]))
        z[[k]]$target<-rep(t,nrow(z[[k]]))
      }
      rr[[j]]<-z
    }
    if(i%%1000L==0L){cat(scenario,'histories',i,'/',nrow(b),'\n');flush.console()}
    rr
  })
  flat<-unlist(all,recursive=FALSE)
  out<-list(baseline=data.table::as.data.table(b),
    scenario=data.table::data.table(scenario=scenario,transition=48,rass_only=FALSE),
    horizon_h=horizon,partial_policy=coverage_policy(required_fraction),
    coverage_required_fraction=required_fraction,clinical_effect_ready=FALSE)
  for(k in c('trace','decisions','status','risk'))out[[k]]<-data.table::rbindlist(lapply(flat,`[[`,k),fill=TRUE)
  stopifnot(nrow(out$status)==4*nrow(b),!anyDuplicated(out$baseline$stay_id),
    all(out$decisions$assessment_h<=horizon),all(out$risk$stop>out$risk$start))
  out
}
coverage_weight_data<-ccw_weight_data
body(coverage_weight_data)<-coverage_replace(body(ccw_weight_data),
  quote(h$partial_policy=='fraction50'),
  quote(h$partial_policy %in% c('fraction25','fraction50','fraction75')))
