# Event-order and adherence core. This module has no file or database access
# and does not fit statistical models. Entry rules, exposure construction and
# weight estimation are implemented in the other pipeline modules.
# Full-node RASS/coverage cannot be known before the node ends: a deviation
# derived from that full node is evaluated at node_end, not retrospectively
# at node_start. Events at or before that assessment retain priority.

event_order_history <- function(nodes, strategy=c("SLA","DTLS"),
                                transition_hours=48, grace_nodes=1L,
                                horizon_hours=144, outcome_hours=Inf,
                                outcome_type="none", intervention_stop_hours=Inf,
                                rass_only=FALSE) {
  strategy <- match.arg(strategy)
  needed <- c("node_start","node_end","state")
  stopifnot(all(needed %in% names(nodes)), horizon_hours>0,
            transition_hours>=0, grace_nodes>=0,
            !is.na(outcome_hours),outcome_hours>=0,
            !is.na(intervention_stop_hours),intervention_stop_hours>=0)
  nd <- nodes[order(nodes$node_start,nodes$node_end),needed,drop=FALSE]
  stopifnot(all(nd$node_end>nd$node_start),all(nd$node_start>=0),
            !anyDuplicated(nd$node_start),
            all(nd$state %in% c("light","deep","low_without","missing","not_applicable")))
  if(nrow(nd)>1L) stopifnot(all(head(nd$node_end,-1)<=tail(nd$node_start,-1)))
  if(rass_only) nd$state[nd$state=="low_without"] <- "deep"
  censor_time <- Inf
  censor_reason <- NA_character_
  consecutive <- 0L
  decisions <- vector("list",nrow(nd))
  n_decisions <- 0L
  treatment_end <- min(horizon_hours,intervention_stop_hours)
  for(i in seq_len(nrow(nd))) {
    t <- nd$node_end[i]
    # A terminal outcome or intervention stop takes precedence even if it
    # coincides with a full-node assessment. Never fill future RASS after death.
    if(outcome_hours<=t || intervention_stop_hours<=t || t>horizon_hours) break
    state <- nd$state[i]
    expected <- if(strategy=="DTLS" && nd$node_start[i]<transition_hours) "deep" else "light"
    adherent <- if(state=="missing") NA else state==expected
    reason <- NA_character_
    if(state=="not_applicable") {
      # This value must be assigned from observed non-MV time, not future success.
      consecutive <- 0L
      adherent <- NA
    } else if(state=="low_without") {
      reason <- "low_RASS_without_qualifying_infusion"
    } else if(!is.na(adherent)) {
      consecutive <- if(adherent) 0L else consecutive+1L
      if(consecutive>grace_nodes) reason <- "ordinary_deviation_after_grace"
    }
    n_decisions <- n_decisions+1L
    decisions[[n_decisions]] <- data.frame(node_start=nd$node_start[i],decision_time=t,
      state=state,expected=expected,adherent=adherent,
      consecutive_deviation=consecutive,censor_now=!is.na(reason))
    if(!is.na(reason)) {censor_time <- t; censor_reason <- reason; break}
  }
  decisions <- if(n_decisions) do.call(rbind,decisions[seq_len(n_decisions)]) else
    data.frame(node_start=numeric(),decision_time=numeric(),state=character(),expected=character(),
               adherent=logical(),consecutive_deviation=integer(),censor_now=logical())
  endpoint_observed <- is.finite(outcome_hours) && outcome_hours<=censor_time
  list(strategy=strategy,censor_time=censor_time,censor_reason=censor_reason,
       endpoint_observed=endpoint_observed,
       endpoint_time=if(endpoint_observed) outcome_hours else NA_real_,
       endpoint_type=if(endpoint_observed) outcome_type else NA_character_,
       risk_followup_stop=min(outcome_hours,censor_time),
       intervention_end=min(treatment_end,outcome_hours,censor_time),
       decisions=decisions)
}

# Timeline for an event endpoint. At an artificial censor, pre-censor person-time
# remains in the data; only subsequent time is removed. This differs from
# deleting every clone that fails to remain adherent through the whole horizon.
make_risk_intervals <- function(history, followup_hours=672) {
  stopifnot(followup_hours>0)
  finish <- min(followup_hours,history$risk_followup_stop)
  if(finish<=0) return(data.frame(start=numeric(),stop=numeric(),event=integer()))
  knots <- sort(unique(c(0,history$decisions$decision_time[
    history$decisions$decision_time<finish],finish)))
  out <- data.frame(start=head(knots,-1),stop=tail(knots,-1),event=0L)
  if(history$endpoint_observed && history$endpoint_time<=followup_hours)
    out$event[nrow(out)] <- 1L
  out
}
