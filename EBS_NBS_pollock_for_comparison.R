pollockEBSwarm<-read_csv(file = here::here("data","POLL_EBS_NBS_shelf_warm.csv"))
spec_data <- pollockEBSwarm %>% 
  clean_names()
weight<-spec_data$weight #set up for nls(), get weight
length<-spec_data$length
nls_LW_fit<-nls(weight~alpha*length^beta,
                data=list(length,weight),
                start=list(alpha=.001,beta=5),
                trace=F, 
                control=list(maxiter = 500))
lw_results <- broom::tidy(nls_LW_fit)
