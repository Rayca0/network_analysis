set.seed(41)
sample_matrix <- matrix(rnorm(10000), nrow = 100, ncol = 100)

#1 Unweighted network (input: data matrix, output: adjacent matrix)
hard_threshold_network <- function(sample_matrix, tau) {
  # check for invalid inputs
  if (!is.numeric(tau)) {
    stop("Error: tau should be a numeric value")
  }
  if (missing(sample_matrix) || missing(tau)) {
    stop("Error: Both similarity_matrix and tau should be provided")
  }
  
  # Perform hard thresholding
  sample_matrix=t(sample_matrix)
  similarity_matrix<-abs(cor(sample_matrix))
  diag(similarity_matrix) <- 0
  pos_key<-similarity_matrix>=tau
  adjacent_matrix <- ifelse(pos_key, 1, 0)
  
  return(adjacent_matrix)
}

#2 Weighted network

# Soft Thresholding Function (input: data matrix, output: adjacent)
soft_threshold_network <- function(sample_matrix, method, tau = NULL, alpha = NULL, beta = NULL) {
  
  # check for invalid inputs
  if (missing(sample_matrix) || missing(method)) {
    stop("Error: Both sample_matrix and method should be provided")
  }
  
  sample_matrix=t(sample_matrix)
  similarity_matrix<-abs(cor(sample_matrix))
  diag(similarity_matrix) <- 0
  
  #1.sigmoid function
  if (method == "sigmoid") {
    if (is.null(tau) || is.null(alpha)) {
      stop("Both tau and alpha are required when using the 'sigmoid' method.")
    }
    adjacent_matrix <- 1 / (1 + exp(-alpha * (similarity_matrix - tau))) 
    
  } #2. power adjacency function
  else if (method == "power") {
    if (is.null(beta)) {
      stop("Beta is required when using the 'power' method.")
    }
    adjacent_matrix <- similarity_matrix^ beta  # Thresholding with the power of beta
  } 
  return(adjacent_matrix)

}

#function to create the log10(p(k)) vs log10(k) (input: adjacent matrix, output: plot and r^2)
scale_free_topology_plot <- function(adjacent_matrix) {
  degree <- rowSums(adjacent_matrix) 
  degree_freq <- table(degree)
  
  degrees <- as.numeric(names(degree_freq))
  freq <- as.numeric(degree_freq)
  total_nodes <- sum(freq)
  
  # Calculate probability distribution p(k)
  probabilities <- freq / total_nodes
  
  log_degrees <- log10(degrees)
  log_probabilities <- log10(probabilities)
  
  # Perform linear regression of log(p(k)) on log(k)
  linear_model <- lm(log_probabilities ~ log_degrees)
  
  # Calculate R-squared value and slope of the regression line
  rsquared <- round(summary(linear_model)$r.squared, 4)
  slope <- round(coef(linear_model)[2], 4)
  
  # Plot log-log plot with regression line, R-squared value, and slope
  plot<-plot(log_degrees, log_probabilities, xlab = "log10(k)", ylab = "log10(p(k))", 
             main = paste("Scale-Free Topology Plot\nR-squared =", rsquared, ", Slope =", slope), 
             col = "blue", pch = 16)
  
  # Add regression line
  abline(linear_model, col = "red")
  
  return(rsquared)
}

#network<-hard_threshold_network(sample_matrix,tau=0.1)
#rsquared<-scale_free_topology_plot(network)

