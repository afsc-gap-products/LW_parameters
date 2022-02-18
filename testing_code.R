species <- rep(seq(1:3), times = 10)
bio_data <- bind_cols(length = rnorm(n = 30, mean = 50, sd = 25), 
                      weight = rnorm(n = 30, mean = 100, sd = 35))

dat <- bind_cols(species = species, bio_data)

# test --------------------------------------------------------------------

# manually get mean by species:
dat %>% filter(species == 1) %>% summarize(mean_wt = mean(weight))
dat %>% filter(species == 2) %>% summarize(mean_wt = mean(weight))
dat %>% filter(species == 3) %>% summarize(mean_wt = mean(weight))
dat %>% summarize(mean_wt = mean(weight))

# create function and loop over function:
get_mean <- function(data, group){
  data %>% 
    dplyr::filter(species == group) %>% 
    summarise(mean_wt = mean(weight))
}

for(i in unique(species)){
  print(get_mean(dat, i))
}

# use purrr map() package
 results <- map(unique(dat$species), get_mean, data = dat)

 
 # get_mean2 <- function(data){
 #   data %>% 
 #     summarise(mean_wt = mean(weight))
 # }
 # results <- map(unique(dat$species), get_mean2, data = dat)

# #NOTES:
map(dat, mean) #takes mean by column...


# old method --------------------------------------------------------------

spec_test <- dat %>% dplyr::filter(species == 1)

#define columns
weight<-spec_test$weight 
length<-spec_test$length 

#run nonlinear least squares regression and print summary report
nls_LW_fit<-nls(weight~alpha*length^beta,data=list(length,weight),start=list(alpha=.001,beta=5),trace=T, control=list(maxiter = 500))
lw_results <- broom::tidy(nls_LW_fit)

# save results:
alpha_result <- lw_results$estimate[1]
beta_result <- lw_results$estimate[2]


# new method --------------------------------------------------------------

get_lw <- function(df, group){ #get l-w relationship from the data frame (df) by species (group)
  
  spec_test <- df %>% dplyr::filter(species == group)
  
  #define columns
  weight<-spec_test$weight 
  length<-spec_test$length 
  
  #run nonlinear least squares regression and print summary report
  nls_LW_fit<-nls(weight~alpha*length^beta,
                  data=list(length,weight),
                  start=list(alpha=.001,beta=5),
                  trace=T, 
                  control=list(maxiter = 500))
  
  lw_results <- broom::tidy(nls_LW_fit) #get tidy results
  print(c("species", group))
  print(lw_results)
  return(lw_results)
}

get_params <- map(unique(dat$species), get_lw, df = dat)

####################

get_params <- map(dat,
                  ~nls(weight~alpha*length^beta,
                       data = list(dat$length, 
                                   dat$weight),
                       start=list(alpha=.001,beta=5),
                       trace=T, 
                       control=list(maxiter = 500))) %>% 
  broom::tidy()


get_params <- map(dat,
                  ~nls(weight~alpha*length^beta,
                       data = list(dat[dat$species %in% species]$length, 
                                   dat[dat$species %in% species]$weight),
                       start=list(alpha=.001,beta=5),
                       trace=T, 
                       control=list(maxiter = 500))) %>% 
  broom::tidy()


get_params <- map(.x = list(dat[dat$species %in% species]$length, 
                            dat[dat$species %in% species]$weight),
                  ~nls(weight~alpha*length^beta,
                       data = .x,
                       start=list(alpha=.001,beta=5),
                       trace=T, 
                       control=list(maxiter = 500))) %>% 
  broom::tidy()
