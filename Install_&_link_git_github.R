# Install git and link R
#https://rfortherestofus.com/2021/02/how-to-use-git-github-with-r
# check in R terminal "git --version" 
# if this command not not recognized, git is not installed
# install git from https://git-scm.com/install/windows
# after that restart R and check again in R terminal "git --version"
install.packages("usethis")
library(usethis)
use_git_config(user.name = "freiakrause", user.email = "freia.krause@plus.ac.at")

# Now need to create git repositroy inside of R project. 
# For this work inside you R project ( when starting R studio, load the Rproject file and then you are working in the project
#library(usethis)
# this now creates git repository and asks to commit the files
use_git()
