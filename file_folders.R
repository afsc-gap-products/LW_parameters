# Function for file folders
# Created by: Caitlin Allen Akselrud
# Contact: caitlin.allen_akselrud@noaa.gov
# Created: 2020-12-18
# Modified: 2021-01-18


# here are the names of the file folders you want to create:
# here is your function:
file_folders <- function(dirs = c("code", "data", "documentation", "figures", "functions", "output"))
{
  
  # for each file folder, do the following:
  for(i in 1:length(dirs)){
    # if the file folder already exists, do nothing
    if(dir.exists(dirs[i])==FALSE){
      # if the file folder does not exist, create it
      dir.create(dirs[i])
    }
  }
}

# create the folders by running the function:
file_folders()
