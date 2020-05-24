## Course 8 Informatica Project
* Christiaan Posthuma
* Rutger Kemperman
* Yaris van Thiel

## Download docker

### Windows
1.  Download [Docker Desktop for Windows](https://hub.docker.com/editions/community/docker-ce-desktop-windows/)
2. Start Docker Desktop
  
### Ubuntu
1. [Install](https://docs.docker.com/engine/install/ubuntu/) Docker (which is the same as the steps below)
1. Remove any old versions of docker 
    ```shell script
    sudo apt-get remove docker docker-engine docker.io containerd runc
    ```
2. Upade the ```apt``` package index and install packages to allow ```apt``` to use a repository over HTTPS:
    ```shell script
    sudo apt-get update
    
    sudo apt-get install \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg-agent \
        software-properties-common
    ```
3. Add Docker's official GPG key:
    ```shell script
     url -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    ```
4. Verify that you now have the key with the fingerprint ```9DC8 5822 9FC7 DD38 854A  E2D8 8D81 803C 0EBF CD88```
    ```shell script
    sudo apt-key fingerprint 0EBFCD88
    
    pub   rsa4096 2017-02-22 [SCEA]
          9DC8 5822 9FC7 DD38 854A  E2D8 8D81 803C 0EBF CD88
    uid           [ unknown] Docker Release (CE deb) <docker@docker.com>
    sub   rsa4096 2017-02-22 [S]
    ```
 5. Set up **stable** repository
    ```shell script
    sudo add-apt-repository \
       "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
       $(lsb_release -cs) \
       stable" 
    ``` 
6. Install Docker Engine
    ```shell script
    sudo apt-get install docker-ce docker-ce-cli containerd.io
   ```

### After installing docker
1. CD to directory 
    ```shell script
     cd bi8a-project-groep5
    ```
2. Run docker-compose file
    ```shell script
    docker-compose up
   ```
3. Go to ```localhost:5000``` in your browser

