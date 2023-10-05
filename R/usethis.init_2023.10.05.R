library(usethis)
library(tidyverse)

# set up ssh
# see GitHub folder in /Documents

# set up .Rprofile
usethis::edit_r_profile()
# paste this into the .Rprofile
# options(
#   usethis.description = list(
#     "Authors@R" = utils::person(
#       "Matthew", "DAntuono",
#       emaiil = "mrd78@georgetown.edu",
#       role = c("aut", "cre"),
#       comment = c(ORCID = "000-0002-1837-7004")
#     )
#   ),
#   usethis.destdir = "/Users/fyendal/Documents/SSD/projects",
#   usethis.overwrite = TRUE
# )

# configure Git, introduce yourself
use_git_config(user.name = "MatthewDAntuono", user.email = "mrd78@georgetown.edu")

# make personal access token for github
usethis::create_github_token()

# set github token in rstudio
gitcreds::gitcreds_set()

# verify your token (personal access token = PAT)
gh_whoami()
usethis::git_sitrep()

# make new RStudio project

# initialize your RStudio project into a git repository
# this will initialize your local project, add, and commit initial files
# to your version history
usethis::use_git()


#### once done, proceed to set up your project ---------------------------------
# more things not done here: https://usethis.r-lib.org/

# make a .readme file
usethis::use_readme_md()

# modify description, add a licence
usethis::use_mit_license("Matthew DAntuono")


#### once set up, push to GitHub -----------------------------------------------

# use_github() will create a remote repository at your GitHub, then push your
# local branch to GitHub
usethis::use_github()

# THEN go and add and commit the new files created
# Git tab in top right window, check boxes, click Commit
# then add a title, note, and push to GitHub

# want to set this up as a package, makes folders cleaner:
path <- "~/Documents/SSD/projects/cognitive.seeds"
usethis::create_package(path)



#### clone from GitHub ---------------------------------------------------------

create_from_github(
  "git@github.com:MatthewDAntuono/basal.leaves.git",
  destdir = "/data/mrd",
  fork = FALSE,
  protocol = "ssh"
)

# alternative if above doesn't work (esp on quorra)
# git clone https://github.com/MatthewDAntuono/basal.leaves.git
