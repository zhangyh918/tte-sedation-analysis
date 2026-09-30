# Exact boundary comparisons at PostgreSQL timestamp precision (microseconds).
# Never compare rounded presentation values such as 1.999 days with 2 days.
seconds_from_hours <- function(x) round(x*3600,6)

is_successful_liberation <- function(end_h,next_start_h,death_h) {
  en <- seconds_from_hours(end_h)
  nx <- seconds_from_hours(next_start_h)
  de <- seconds_from_hours(death_h)
  !is.na(de) & en<=28*86400 & nx>en+48*3600 & de>en+48*3600
}

merge_episodes <- function(a,b) {
  ok <- is.finite(a)&is.finite(b)&b>a
  a<-a[ok];b<-b[ok]
  if(!length(a)) return(data.table::data.table(start_h=numeric(),end_h=numeric()))
  ii<-order(a,b);a<-a[ii];b<-b[ii]
  starts<-ends<-numeric(length(a));k<-1L;starts[k]<-a[1];ends[k]<-b[1]
  if(length(a)>1L) for(j in 2:length(a)) {
    if(a[j]<=ends[k]) ends[k]<-max(ends[k],b[j]) else {
      k<-k+1L;starts[k]<-a[j];ends[k]<-b[j]
    }
  }
  data.table::data.table(start_h=starts[seq_len(k)],end_h=ends[seq_len(k)])
}
