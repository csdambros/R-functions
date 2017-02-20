## Functions to send ssh ids to other nodes in the cluster
# Grants passwordless access within local network
# Required to spawn MPI processes across nodes

passSSHid<-function(passwd,IPfile="/etc/openmpi/openmpi-default-hostfile"){

# Create .ssh directory if it does not exist
system("mkdir ~/.ssh")
  
# Change folder permissions (only owner can access)
system("chmod -R 700 ~/.ssh")
  
# Generate ssh key
system('ssh-keygen -t rsa -f ~/.ssh/id_rsa -q -P ""')

# Remove any known hosts if they exist
system("rm ~/.ssh/known_hosts")

# Read file with host IPs
IPs<-read.csv(IPfile)

# 
for(i in 1:nrow(IPs)){

  #Skip if IP has a hashtag
  if(regexpr("#",IPs[i,1])[1]==1){cat('skipping',as.character(IPs[i,1]),"\n");next}

  # create code to be passed
  code<-paste0('sshpass -p "',passwd,'" ssh-copy-id -f ',IPs[i,1])

  # copy ssh id to all computers
  system(code)

  cat('IP',as.character(IPs[i,1]),'done',"\n")

}

}


# USE
#passSSHid(passwd = "XXX")

