# load file
nutrients <- read.csv("nutrients_vs_ph.csv", header = TRUE, sep = ",")

# plot
ph <- nutrients$ph

nutrients <- nutrients[, -1]

# drop na
nutrients <- na.omit(nutrients)


print(nutrients)
