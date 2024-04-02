## Library aphid is written by Shaun Wilkinson
library("aphid")
states <- c("Begin", "Fair", "Loaded")
residues <- paste(1:6)
### Define transition probability matrix A
A <- matrix(c(0, 0, 0, 0.99, 0.95, 0.1, 0.01, 0.05, 0.9), nrow = 3)
dimnames(A) <- list(from = states, to = states)
### Define emission probability matrix E
E <- matrix(c(rep(1/6, 6), rep(1/10, 5), 1/2), nrow = 2, byrow = TRUE)
dimnames(E) <- list(states = states[-1], residues = residues)
### Create the HMM object
x <- structure(list(A = A, E = E), class = "HMM")
### Plot the model
plot(x, textexp = 1.5)
### Optionally add the transition probabilities as text
text(x = 0.02, y = 0.5, labels = "0.95")
text(x = 0.51, y = 0.5, labels = "0.90")
text(x = 0.5, y = 0.9, labels = "0.05")
text(x = 0.5, y = 0.1, labels = "0.10")


data(casino)
### The actual path is stored in the names attribute of the sequence
actual <- c("F", "L")[match(names(casino), c("Fair", "Loaded"))]
### Find the predicted path
vit1 <- Viterbi(x, casino)
predicted <- c("F", "L")[vit1$path + 1]
### Note the path element of the output Viterbi object is an integer vector
### the addition of 1 to the path converts from C/C++ to R's indexing style


casino.post <- posterior(x, casino)
plot(1:300, seq(0, 1, length.out = 300), type = "n", xlab = "Roll number",
     ylab = "Posterior probability of dice being fair")
starts <- which(c("L", actual) == "F" & c(actual, "F") == "L")
ends <- which(c("F", actual) == "L" & c(actual, "L") == "F") - 1
for(i in 1:6) rect(starts[i], 0, ends[i], 1, col = "grey", border = NA)
lines(1:300, casino.post[1, ])


y <- deriveHMM(list(casino), logspace = FALSE)
plot(y, textexp = 1.5)

### Optionally add the transition probabilities as text
text(x = 0.02, y = 0.5, labels = round(y$A["Fair", "Fair"], 2))
text(x = 0.51, y = 0.5, labels = round(y$A["Loaded", "Loaded"], 2))
text(x = 0.5, y = 0.9, labels = round(y$A["Fair", "Loaded"], 2))
text(x = 0.5, y = 0.1, labels = round(y$A["Loaded", "Fair"], 2))

data(globins)
globins

globins.PHMM <- derivePHMM(globins, residues = "AMINO", pseudocounts = "Laplace")
plot(globins.PHMM)

################################################################################
###################################################################### REFERENCE
################################################################################
###### https://cran.r-project.org/web/packages/aphid/vignettes/aphid-vignette.html
