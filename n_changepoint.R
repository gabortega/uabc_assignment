# Clean the environment
rm(list = ls())

library(boot)
library(coda)
library(parallel)
data(coal)

parallel_cluster = makePSOCKcluster(detectCores(logical=TRUE)-1)

iterations = 10000
burn_in = 1000

x = 1851:1962
y = tabulate(floor(coal[[1]]))
coal_data = y[x]
N = length(coal_data)
prior_alpha = 2
prior_beta = 1

############################################################################################
############################################################################################
############################################################################################
############################################################################################

n_conditional_dist_sample = function(changepoint_vector, lambda_vector, data) {
  sample_vector = head(seq_along(data), -1)[(changepoint_vector[1]+1):(changepoint_vector[3]-1)]
  if (length(sample_vector) > 1) {
    all_prob = 
      sapply(sample_vector, 
             function(n) {
               partial_sum_1 = sum(data[(changepoint_vector[1]+1):n])
               partial_sum_2 = sum(data[(n+1):changepoint_vector[3]])
               return(exp((log(lambda_vector[1]) * partial_sum_1) + (-(n - changepoint_vector[1]) * lambda_vector[1]) +
                            (log(lambda_vector[2]) * partial_sum_2) + (-(changepoint_vector[3] - n) * lambda_vector[2])))
             })
    return(sample(sample_vector, 1, prob=(all_prob/sum(all_prob)), replace=TRUE))
  } else {
    return(sample_vector)
  }
}

lambda_conditional_dist_sample = function(changepoint_vector, data, alpha, beta) {
  partial_sum = sum(data[(changepoint_vector[1]+1):changepoint_vector[2]])
  return(rgamma(1, (alpha + partial_sum), (beta + changepoint_vector[2] - changepoint_vector[1])))
}

gibbs_iter = function(changepoint_vector, lambda_vector, data, alpha, beta, skip_changepoints=0, skip_lambdas=0) {
  changepoint_vector = c(0, changepoint_vector, N)
  for (i in 1:max_changepoints) {
    if (i > skip_changepoints) {
      selected_changepoints = changepoint_vector[i:(i+2)]
      selected_lambdas = lambda_vector[i:(i+1)]
      #
      changepoint_vector[(i+1)] = n_conditional_dist_sample(selected_changepoints, selected_lambdas, data)
    }
    if (i > skip_lambdas) {
      lambda_vector[i] =  lambda_conditional_dist_sample(changepoint_vector[i:(i+1)], data, alpha, beta)
    }
  }
  if ((max_changepoints+1) > skip_lambdas) {
    lambda_vector[(max_changepoints+1)] = lambda_conditional_dist_sample(tail(changepoint_vector, 2), data, alpha, beta)
  }
  #
  changepoint_vector = head(tail(changepoint_vector, -1), -1)
  #
  return(c(changepoint_vector, lambda_vector))
}

run_gibbs = function(run_iterations, data, alpha, beta, const_changepoint_vector=NULL, const_lambda_vector=NULL) {
  if (!is.null(const_changepoint_vector)) {
    changepoint_vector = const_changepoint_vector
  } else {
    changepoint_vector = sample(head(seq_along(data), -max_changepoints)[1:(N-max_changepoints)], 1)
  }
  if (length(changepoint_vector) < max_changepoints) {
    for (n in (max_changepoints-length(changepoint_vector)):1) {
      samples = head(seq_along(data), -n)[(tail(changepoint_vector, 1)+1):(N-n)]
      changepoint_vector = append(changepoint_vector, ifelse((length(samples) > 1), sample(samples, 1), samples))
    }
  }
  #
  if (!is.null(const_lambda_vector)) {
    lambda_vector = c(const_lambda_vector, rgamma((max_changepoints+1)-length(const_lambda_vector), alpha, beta))
  } else {
    lambda_vector = rgamma(max_changepoints+1, alpha, beta)
  }
  #
  if (!is.null(const_changepoint_vector)) {
    skip_changepoints = length(const_changepoint_vector)
  } else {
    skip_changepoints = 0
  }
  #
  if (!is.null(const_lambda_vector)) {
    skip_lambdas = length(const_lambda_vector)
  } else {
    skip_lambdas = 0
  }
  #
  samples = matrix(nrow=run_iterations, ncol=((max_changepoints*2)+1))
  samples[1,] = c(changepoint_vector, lambda_vector)
  for (i in 2:run_iterations) {
    changepoint_vector = samples[(i - 1), 1:max_changepoints]
    lambda_vector = samples[(i - 1), (max_changepoints+1):((max_changepoints*2)+1)]
    #
    samples[i,] = gibbs_iter(changepoint_vector, lambda_vector, data, alpha, beta, 
                             skip_changepoints=skip_changepoints,
                             skip_lambdas=skip_lambdas)
  }
  #
  colnames(samples) = c(
    sapply(seq_along(changepoint_vector), function(i) paste0("T_", i)),
    sapply(seq_along(lambda_vector), function(i) paste0("lambda_", i)))
  return(as.data.frame(samples))
}

############################################################################################
############################################################################################
############################################################################################
############################################################################################

log_likelihood = function(changepoint_vector, lambda_vector, data) {
  changepoint_vector = c(0, changepoint_vector, N)
  return(
    sum(sapply(1:N, function(i) {
      return(dpois(data[i], 
                   lambda_vector[tail(which(i>changepoint_vector), 1)], log = TRUE))
    }))
  )
}

log_prior = function(lambda_vector) {
  return(sum(dgamma(lambda_vector, prior_alpha, prior_beta, log=TRUE), log(1/choose(111, max_changepoints))))
}

tau_posterior = function(tau, changepoint_vector, lambda_vector, data) {
  sample_vector = head(seq_along(data), -1)[(changepoint_vector[1]+1):(changepoint_vector[3]-1)]
  if (length(which(sample_vector==tau))) {
    all_prob = sapply(sample_vector, 
                      function(n) {
                        partial_sum_1 = sum(data[(changepoint_vector[1]+1):n])
                        partial_sum_2 = sum(data[(n+1):changepoint_vector[3]])
                        return(exp((log(lambda_vector[1]) * partial_sum_1) + (-(n - changepoint_vector[1]) * lambda_vector[1]) +
                                     (log(lambda_vector[2]) * partial_sum_2) + (-(changepoint_vector[3] - n) * lambda_vector[2])))
                      })
    return((all_prob[which(sample_vector==tau)])/sum(all_prob))
  } else {
    return(0)
  }
}

lambda_posterior = function(changepoint_vector, lambda, data, alpha, beta){
  partial_sum = sum(data[(changepoint_vector[1]+1):changepoint_vector[2]])
  return(exp(dgamma(lambda, (alpha + partial_sum), (beta + changepoint_vector[2] - changepoint_vector[1]), log=TRUE)))
}

log_marginal_likelihood = function(data, postsample) {
  print("Approximating local maxima log likelihood for this model...")
  postsample_index = which.max(sapply(1:nrow(postsample), 
                                      function(i) {
                                        return(log_likelihood(unlist(postsample[i, 1:max_changepoints]), 
                                                              unlist(postsample[i, (max_changepoints+1):((max_changepoints*2)+1)]),
                                                              data))
                                      }))
  max_changepoint_vector = unlist(postsample[postsample_index, 1:max_changepoints])
  # print(max_changepoint_vector)
  max_lambda_vector = unlist(postsample[postsample_index, (max_changepoints+1):((max_changepoints*2)+1)])
  # print(max_lambda_vector)
  #
  print("Calculating likelihood from Gibbs sample...")
  likelihood = log_likelihood(max_changepoint_vector, max_lambda_vector, data)
  # print(likelihood)
  #
  print("Calculating prior from Gibbs sample...")
  prior = log_prior(max_lambda_vector)
  # print(prior)
  #
  posterior = 0
  reduced_postsample = postsample
  for (i in 1:max_changepoints) {
    print(paste0("Calculating tau_", i, " from Reduced Gibbs sample..."))
    posterior = posterior + log(mean(sapply(1:nrow(reduced_postsample),
                                            function(j) {
                                              changepoint_vector = c(0, unlist(reduced_postsample[j, 1:max_changepoints]), N)
                                              lambda_vector = unlist(reduced_postsample[j, (max_changepoints+1):((max_changepoints*2)+1)])
                                              return(
                                                tau_posterior(
                                                  max_changepoint_vector[i],
                                                  changepoint_vector[i:(i+2)],
                                                  lambda_vector[i:(i+1)], data))
                                            })))
    # print(posterior)
    if (i == 1) {
      lambda_vector = NULL
    } else {
      lambda_vector = max_lambda_vector[1:(i-1)]
    }
    print("Running Reduced Gibbs sampler...")
    reduced_gibbs_sample = run_gibbs(iterations, coal_data, prior_alpha, prior_beta, max_changepoint_vector[1:i], lambda_vector)
    reduced_postsample = reduced_gibbs_sample[min(iterations/10, burn_in):iterations,]
    #
    print(paste0("Calculating lambda_", i, " from Gibbs sample..."))
    posterior = posterior + log(mean(sapply(1:nrow(reduced_postsample), 
                                            function(j) {
                                              changepoint_vector = c(0, unlist(reduced_postsample[j, 1:max_changepoints]), N)
                                              return(
                                                lambda_posterior(changepoint_vector[i:(i+1)],
                                                                 max_lambda_vector[i], data, prior_alpha, prior_beta))
                                            })))
    # print(posterior)
    if (i < max_changepoints) {
      print("Running Reduced Gibbs sampler...")
      reduced_gibbs_sample = run_gibbs(iterations, coal_data, prior_alpha, prior_beta, max_changepoint_vector[1:i], max_lambda_vector[1:i])
      reduced_postsample = reduced_gibbs_sample[min(iterations/10, burn_in):iterations,]
    }
  }
  print(paste0("Calculating lambda_", (max_changepoints+1), " from Gibbs sample..."))
  posterior = posterior + log(lambda_posterior(c(max_changepoint_vector[max_changepoints], N),
                                               max_lambda_vector[(max_changepoints+1)], data, prior_alpha, prior_beta))
  # print(posterior)
  
  return(likelihood + prior - posterior)
}

############################################################################################
############################################################################################
############################################################################################
############################################################################################

write_data = function(table, filename) {
  curr_time = format(Sys.time(), "%Y_%m_%d_T%H_%M_%S")
  write.csv(as.data.frame(table),
            file=
              paste0("./data/", filename, "_", curr_time, ".csv"),
            row.names=FALSE)
}

modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

############################################################################################
############################################################################################
############################################################################################
############################################################################################

if (!dir.exists("./data")) {
  dir.create("./data", showWarnings=TRUE)
}

runs = 6

log_marginal_likelihoods = matrix(double((N-1)*(runs)), nrow=(N-1), ncol=(runs))
# log_marginal_likelihoods=as.matrix(read.csv("./log_marginal_likelihoods_100_chains.csv", header=TRUE))

starting_changepoint_iteration = 54
max_changepoint_iterations = 54
iterations_increment = 1
#
for (i in seq(starting_changepoint_iteration, max_changepoint_iterations, iterations_increment)) {
  gc(verbose = FALSE)
  #
  max_changepoints=i
  #
  clusterExport(cl=parallel_cluster,
                c(# Variables
                  "iterations", "burn_in", "x", "y", "coal_data", "N", "prior_alpha", "prior_beta", "max_changepoints",
                  # Methods
                  "n_conditional_dist_sample", "lambda_conditional_dist_sample", "gibbs_iter", "run_gibbs", 
                  "log_likelihood", "log_prior", "lambda_posterior", "tau_posterior", "log_marginal_likelihood"))
  #
  print(paste("Executing", runs, "parallel chains for", max_changepoints, "changepoint(s)..."))
  t1 = proc.time()
  
  # Parallel runs
  log_marginal_likelihoods[max_changepoints,] =
    parSapply(parallel_cluster, 1:runs,
              function(run_iteration) {
                print("#############################################################")
                
                print(paste("Running Complete Gibbs sampler for", max_changepoints, "changepoint(s)..."))
                
                # Single runs
                ###########################################################################
                # t1 = proc.time()
                # gibbs_sample = run_gibbs(iterations, coal_data, prior_alpha, prior_beta)
                # t2 = proc.time()
                # 
                # print(t2-t1)
                # 
                # # skipping burn-in
                # postsample = gibbs_sample[min(iterations/10, burn_in):iterations,]
                # 
                # # for (i in 1:(max_changepoints)) {
                # #   hist(gibbs_sample[,i], main=colnames(gibbs_sample)[i])
                # # }
                # 
                # # Draw from the posterior predictive distribution
                # predsample = matrix(0, nrow=length(postsample), ncol=N)
                # for (i in 1:length(postsample)) {
                #   changepoint_vector = c(0, unlist(postsample[i, 1:max_changepoints]), N)
                #   lambda_vector = unlist(postsample[i, (max_changepoints+1):((max_changepoints*2)+1)])
                #   predsample[i,] = cumsum(unlist(sapply(1:(max_changepoints+1), function(j) {
                #     rpois(changepoint_vector[[(j+1)]] - changepoint_vector[[j]], lambda_vector[[j]])
                #   })))
                # }
                # 
                # # Evaluate posterior predictive quantiles
                # predq = matrix(0,nrow=3,ncol=length(x))
                # for (i in 1:length(x)) {
                #   predq[,i] = quantile(predsample[,i],probs=c(0.025,0.5,0.975))
                # }
                # 
                # # t1 = proc.time()
                # # result = log_marginal_likelihood(coal_data, postsample)
                # # t2 = proc.time()
                # # print(paste("Log marginal likelihood:", result))
                # # print(t2-t1)
                # 
                # plot(x, cumsum(coal_data), type="l", xlab="Year", ylab="Cumulative number of disasters")
                # lines(x, predq[1,], col="green")
                # lines(x, predq[2,], col="red")
                # lines(x, predq[3,], col="green")
                # 
                # legend("topleft", legend=c("95% Prediction Interval", "Prediction Median", "Observations"),
                #        col=c("green","red","black"),
                #        lty = c(1, 1, 1))
                # 
                # # text(x = 1930, y = 65,
                # #      label=bquote(.("MAP(")~"T"[1]~.(paste0(")=",1851+modes(postsample[,"T_1"])))))
                # # text(x = 1930, y = 50,
                # #      label=bquote(.("MAP(")~"T"[2]~.(paste0(")=",1851+modes(postsample[,"T_2"])))))
                # # text(x = 1930, y = 35,
                # #      label=bquote(.("MAP(")~"T"[3]~.(paste0(")=",1851+modes(postsample[,"T_3"])))))
                # # text(x = 1930, y = 20,
                # #      label=bquote(.("MAP(")~"T"[4]~.(paste0(")=",1851+modes(postsample[,"T_4"])))))
                # # text(x = 1930, y = 5,
                # #      label=bquote(.("MAP(")~"T"[5]~.(paste0(")=",1851+modes(postsample[,"T_5"])))))
                # # 
                # # text(x = 1945, y = 65,
                # #      label=bquote(.("MAP(")~"T"[6]~.(paste0(")=",1851+modes(postsample[,"T_6"])))))
                # # text(x = 1945, y = 50,
                # #      label=bquote(.("MAP(")~"T"[7]~.(paste0(")=",1851+modes(postsample[,"T_7"])))))
                # # text(x = 1945, y = 35,
                # #      label=bquote(.("MAP(")~"T"[8]~.(paste0(")=",1851+modes(postsample[,"T_8"])))))
                # # text(x = 1945, y = 20,
                # #      label=bquote(.("MAP(")~"T"[9]~.(paste0(")=",1851+modes(postsample[,"T_9"])))))
                # # text(x = 1945, y = 5,
                # #      label=bquote(.("MAP(")~"T"[10]~.(paste0(")=",1851+modes(postsample[,"T_10"])))))
                # 
                # # plot(mcmc(postsample), smooth=TRUE, auto.layout=TRUE, ask=FALSE)
                ###########################################################################
                
                # Parallel runs
                ###########################################################################
                gibbs_sample = run_gibbs(iterations, coal_data, prior_alpha, prior_beta)
                postsample = gibbs_sample[min(iterations/10, burn_in):iterations,]
                return(log_marginal_likelihood(coal_data, postsample))
                ###########################################################################
              }
    )
  t2 = proc.time()
  print(t2-t1)
  write_data(log_marginal_likelihoods, "log_marginal_likelihoods")
}
stopCluster(parallel_cluster)
#
mean_log_marginal_likelihoods = double(N-1)
for (i in 1:nrow(log_marginal_likelihoods))
  mean_log_marginal_likelihoods[i]=mean(log_marginal_likelihoods[i,])

############################################################################################
############################################################################################
############################################################################################
############################################################################################
df = data.frame((log_marginal_likelihoods)[c(5,10,21,32,43),])

df$mean = apply(df, 1, mean)
df$stdv = apply(df, 1, sd)

ggplot(df[,c(101,102)], aes(x=1:5, y=mean)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv), width=.2,
                position=position_dodge(.9)) +
  labs(x="Number of changepoints considered", y="Log Marginal Likelihood for a model") +
  coord_cartesian(ylim=c(-174,-177))


# df = data.frame(exp(log_marginal_likelihoods)[1:5,])
# 
# df$median = apply(df, 1, median)
# df$conf_int_lower = apply(df, 1, function(x) {t.test(x)$conf.int[1]})
# df$conf_int_upper = apply(df, 1, function(x) {t.test(x)$conf.int[2]})
# 
# ggplot(df[,101:103], aes(x=c(5,10,21,32,43), y=median)) +
#   geom_bar(stat="identity", color="green", position=position_dodge(), aes(ymin=conf_int_lower, ymax=conf_int_upper)) +
#   geom_bar(stat="identity", color="green", position=position_dodge(), aes(ymin=conf_int_lower, ymax=conf_int_upper)) +
#   geom_bar(stat="identity", color="green", position=position_dodge(), aes(ymin=conf_int_lower, ymax=conf_int_upper)) +
#   geom_errorbar(aes(ymin=conf_int_lower, ymax=conf_int_upper), width=.2,
#                 position=position_dodge(.9)) +
#   labs(x="Number of changepoints considered", y="Marginal Likelihood for a model")

# log_marginal_likelihoods = rep(0,5)
# 
# marginal_likelihood_iter = function(changepoint_vector) {
#   changepoint_vector = c(0, changepoint_vector, N)
#   logsum = 0
#   for (i in 1:(max_changepoints+1)) {
#     partial_sum = sum(coal_data[(changepoint_vector[i]+1):changepoint_vector[i+1]])
#     length = changepoint_vector[(i+1)] - changepoint_vector[i]
#     logsum = logsum + ((-partial_sum - 2) * log(length + 1)) + lgamma(2 + partial_sum)
#   }
#   return(exp(logsum))
# }
# 
# marginal_likelihood = function (x, m)
# {
#   stopifnot(length(m) == 1L, is.numeric(m))
#   if (m < 0)
#     stop("m < 0", domain = NA)
#   n = length(seq_len(x))
#   if (n < m)
#     stop("n < m", domain = NA)
#   print("Calculating log marginal likelihood for this model...")
#   result = 0
#   m = as.integer(m)
#   e = 0
#   h = m
#   a = seq_len(m)
#   count = as.integer(round(choose(n, m)))
#   result = result + marginal_likelihood_iter(a)
#   if (m > 0) {
#     i = 2L
#     nmmp1 = n - m + 1L
#     while (a[1L] != nmmp1) {
#       if (e < n - h) {
#         h = 1L
#         e = a[m]
#         j = 1L
#       }
#       else {
#         e = a[m - h]
#         h = h + 1L
#         j = 1L:h
#       }
#       a[m - h + j] = e + j
#       result = result + marginal_likelihood_iter(a)
#       i = i + 1L
#     }
#   }
#   #
#   result = result / (prod(factorial(coal_data)) * choose(N-1, max_changepoints))
# }
# 
# max_changepoints=5
# 
# t1 = proc.time()
# log_marginal_likelihoods[max_changepoints] = log(marginal_likelihood(N-1, max_changepoints))
# t2 = proc.time()
# # print(paste("Log marginal likelihood:", result))
# print(t2-t1)
# 
# log_marginal_likelihoods = c(-176.4679 -175.6190 -175.3718 -175.2496 -175.2511)
# 
# df = data.frame((log_marginal_likelihoods))
# 
# ggplot(df, aes(x=1:5, y=log_marginal_likelihoods)) + 
#   geom_bar(stat="identity", position=position_dodge()) +
#   # geom_errorbar(aes(ymin=mean-stdv, ymax=mean+stdv), width=.2,
#   #               position=position_dodge(.9)) +
#   labs(x="Number of changepoints considered", y="Log Marginal Likelihood for a model") +
#   coord_cartesian(ylim=c(-174,-177))
