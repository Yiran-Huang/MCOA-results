#Install three packages

packages <- c("ggplot2", "gridExtra", "glmnet")
install_if_missing <- function(pack) {
  if (!require(pack, character.only = TRUE)) {
    install.packages(pack, dependencies = TRUE)
    library(pack, character.only = TRUE)
  }
}

# Install the packages
lapply(packages, install_if_missing)

# Confirm installation
installed_packages <- installed.packages()
details=installed_packages[installed_packages[, "Package"] %in% packages, ]
print(rownames(details))