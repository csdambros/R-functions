## Functions to send ssh ids to other nodes in the cluster
# Grants passwordless access within local network
# Required to spawn MPI processes across nodes

passSSHid<-function(passwd=NULL,IPfile="/etc/openmpi/openmpi-default-hostfile"){

if(is.null(passwd)){passwd <-.rs.askForPassword("Digite sua senha/Type your password")}

# Create .ssh directory if it does not exist

system("mkdir ~/.ssh",ignore.stdout = TRUE,ignore.stderr = TRUE)

# Change folder permissions (only owner can access)
system("chmod -R 700 ~/.ssh",ignore.stdout = TRUE)
  
# Generate ssh key
system('ssh-keygen -t rsa -f ~/.ssh/id_rsa -q -P ""',ignore.stdout = TRUE)

# Remove any known hosts if they exist
system("rm ~/.ssh/known_hosts",ignore.stdout = TRUE)

# Read file with host IPs
IPs<-read.csv(IPfile,header = FALSE)

# 
for(i in 1:nrow(IPs)){

  IPo<-sub(":.*","",IPs[i,1])
  #Skip if IP has a hashtag
  if(regexpr("#",IPs[i,1])[1]==1){cat('skipping Slave ',i,"\n");next}

  # create code to be passed
  code<-paste0('sshpass -p "',passwd,'" ssh-copy-id -f ',IPo)

  # copy ssh id to all computers
  system(code,ignore.stdout = TRUE,ignore.stderr = TRUE)

  cat('Slave',i,' authorized',"\n")

}

}


# USE
#passSSHid(passwd = "XXX")

