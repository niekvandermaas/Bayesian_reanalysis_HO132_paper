set.seed(123)
library("rjags")
library("tidyverse")
library("tidybayes")
library("survminer")
library("survival")
library('bayesplot')
library("ggthemes")

#data 
data_150 <- read.csv("data_150_5.csv", sep="", stringsAsFactors=TRUE)
data_300 <- read.csv("data_300_5.csv", sep="", stringsAsFactors=TRUE)
data_450 <- read.csv("data_450_5.csv", sep="", stringsAsFactors=TRUE)
data_600 <- read.csv("data_600_5.csv", sep="", stringsAsFactors=TRUE)

#get HR
cph <- coxph(Surv(efs_i2,efs_i2_i) ~ group_2, data = data_600[data_600$trialnr %in% "HO132",])
summary(cph)

#plot the curves of EFS
fit1 <- survfit(Surv(efs_i2,efs_i2_i) ~ group, data = data_600)

ggsurvplot(
  fit1,
  data = data_600,
  pval = T
)

##############################################
#Bayesian interim analysis event free survival
##############################################
# The model specification
file.show("survival_model_cp.txt")

#make a list of all dataframes
interim_list <- list(data_150 = data_150 %>% arrange(trialnr),
                     data_300 = data_300 %>% arrange(trialnr),
                     data_450 = data_450 %>% arrange(trialnr),
                     data_600 = data_600 %>% arrange(trialnr))

# data preperations
censored <- list()
is_censored <- list()
efs <- list()
for (i in 1:length(interim_list)) {
#censored data
censored[[i]] <- ifelse(interim_list[[i]]["efs_i2_i"] == 1 |is.na(interim_list[[i]]["efs_i2_i"]), FALSE, TRUE)
is_censored[[i]] <- censored[[i]]*1
efs[[i]] <- interim_list[[i]]["efs_i2"]
efs[[i]][censored[[i]]] <- NA
}

t_cen <- list()
group <- list()
n1 <- list()
n2 <- list()
for (i in 1:length(censored)){
t_cen[[i]] <- rep(0, length(censored[[i]]))
t_cen[[i]][censored[[i]]] <- interim_list[[i]]["efs_i2"][censored[[i]]]
t_cen[[i]][is.na(t_cen[[i]])] <- 0
group[[i]] <- ifelse(interim_list[[i]]["group"]== "132 treatment", TRUE, FALSE)
n1[[i]] <- nrow(interim_list[[i]] %>% filter(trialnr == "HO102"))
n2[[i]] <- nrow(interim_list[[i]])
}

aml_data <- list()
mcmc_samples <- list()
results <- list()

N <- list("1" = 800,
          "2" = 800,
          "3" = 800,
          "4" = 800)
group_2 <- list("1" <- rep(c(TRUE, FALSE), 800),
                "2" <- rep(c(TRUE, FALSE), 800),
                "3" <- rep(c(TRUE, FALSE), 800),
                "4" <- rep(c(TRUE, FALSE), 800))

for(i in 1:length(interim_list)){
  # put the data together for JAGS
  aml_data[[i]] <- list(time = efs[[i]][,1] +0.001,
                        t_cen = t_cen[[i]] + 0.001,
                        is_censored = is_censored[[i]][,1],
                        n1 = n1[[i]],
                        n2 = n2[[i]],
                        group = group[[i]][,1],
                        N = N[[i]],
                        group_2 = group_2[[i]],
                        #"Spike and Slab" hyperparameters
                        l.slab= 1.0E-4, #S_l (lower bound of slab) 
                        u.slab= 2, #S_u (upper bound of slab)
                        p0= 0.5, #p_0 (prior probability of being in the slab)
                        spike=5000 #entry is Spike for logit_theta1
  )
}


#Run the model
for(i in 1:length(interim_list)){
#which paramaters to monitor
params_to_monitor <- c("hr_1", "p1", "p2", "p3")

# Running the model
model <- jags.model(file = "survival_model_cp.txt", aml_data[[i]], 
                    n.chains = 3, n.adapt= 1000)
update(model, 1000); # Burning 1000 samples
mcmc_samples[[i]] <- coda.samples(model, params_to_monitor, n.iter= 50000,
                                  inits = list(.RNG.name = "base::Wichmann-Hill",
                                               .RNG.seed = 2022,
                                               thin = 50))

results[[i]] <- mcmc_samples[[i]][[1]]
}

summarise_draws(mcmc_samples[[1]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 300 
summarise_draws(mcmc_samples[[2]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 450
summarise_draws(mcmc_samples[[3]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 600
summarise_draws(mcmc_samples[[4]], mean, sd, ~quantile(.x, probs = c(0.5, 0.025,0.975))) #n = 600

#function for auc posterior
calc_auc_post <- function(data, threshold) {
  prob <- sum(data<threshold)/length(data)
  prb <- paste0(sapply(prob, round, 2)*100,"%")
  return(prb)
}

# results 
hr_1 <- list()
prob <- list()
for (hr in 1:length(results)){
  hr_1[[hr]] <- results[[hr]][,1]
  
  prob[[hr]] <- calc_auc_post(hr_1[[hr]], 0.76)
}


# Make density plots
prob <- unlist(prob)
post <- data.frame(sample(hr_1[[1]], 1000), sample(hr_1[[2]], 1000), 
                   sample(hr_1[[3]], 1000), sample(hr_1[[4]], 1000))

colnames(post) <- c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3", "Interim analysis 4")

dat_text <- data.frame(
  label = prob,
  key   = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3", "Interim analysis 4")
)

# plot posteriors with 95% Bayesian credible intervals
p <- post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value),
               fill = "grey92", adjust=2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  xlab(NULL) +
  facet_wrap(~key, scales = "free", ncol = 2)

#extract values for fill
to_fill <- data_frame(
  x = ggplot_build(p)$data[[1]]$x,
  y = ggplot_build(p)$data[[1]]$y,
  key = factor(ggplot_build(p)$data[[1]]$PANEL, levels = c(1,2,3,4), 
               labels = c("Interim analysis 1", "Interim analysis 2", "Interim analysis 3", "Interim analysis 4")))


#plot HR posterio distribution
post %>% 
  gather() %>% 
  
  mutate(key = factor(key)) %>% 
  
  ggplot() +
  geom_density(aes(x = value), trim = TRUE,
               fill = "#86D2ED", size = 0.8, adjust=2, alpha = 0.3) +
  geom_vline(aes(xintercept = 0.76), linetype = "dashed") +
  stat_pointinterval(aes(x = value, y = 0),
                     point_interval = mode_hdi, .width = c(.95, .5)) +
  scale_y_continuous(NULL, breaks = NULL) +
  scale_x_continuous(limits = c(0, 2), breaks = c(0,0.5,0.76,1,1.5,2)) + 
  xlab(NULL) + geom_text(
    data    = dat_text,
    mapping = aes(x = 0.75, y = 0.5, label = label),
    hjust = 1, vjust = -1.5) +
  geom_area(data = to_fill[to_fill$x <= 0.76, ], 
            aes(x=x, y=y), fill = "#0C2074",
            alpha = 0.3) +
  facet_wrap(~key, scales = "fixed", ncol = 1, strip.position = "left") + 
  theme_fivethirtyeight() + 
  theme(
    strip.text.y.left = element_text(angle = 90),
    plot.background = element_rect(fill = " white"),
    panel.background = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white"),
    axis.title.x =  element_text()
  ) +
  xlab("Hazard ratio")