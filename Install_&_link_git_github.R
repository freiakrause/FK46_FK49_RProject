# Install git and link R
#https://rfortherestofus.com/2021/02/how-to-use-git-github-with-r
# check in R terminal "git --version" 
#if not recognized, git not installed
# install git from https://git-scm.com/install/windows
#after that restart R and check again "git --version"
install.packages("usethis")
library(usethis)
use_git_config(user.name = "freiakrause", user.email = "freia.krause@plus.ac.at")