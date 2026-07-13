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

#### Link local git/Rstudio with GitHub ------
# Create an account in github (using email adress specified abive for git)
# Now we need to create a personal access token to be able to link Rstudio and Github
# On the website opend by this code give PAT a name and copy it. DOnt loos it, it will never appear again
create_github_token()
#safe it in a save space not here
#token expires in october 2026

#now need to give token to r so that it can talk to github
install.packages("gitcreds")
library(gitcreds)
#enter the token when asked for password or token
gitcreds_set()

#### Link Rstudio with github ------
use_github()


#### Using GitHub -----
#When you are done working (pre day or when ever you did important changes) 
#generate a commit. Give a short not about the content of commit and then commit

#after committing, we can push the changes to github to save the new versions there

#when you start working the next day and are unsure of you local version of the script but now, you have committed and pushed the latest version
#pull from git hub
#very important when working with multiple people on one script
#little less important if you are alone on the script and know about the changes

