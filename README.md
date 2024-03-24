## You can set up the analysis tools in an environment to ensure stability of dependencies
### Install mamba
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
```
 run the command below to install conda
 ```
bash Mambaforge-Linux-x86_64.sh
```
Close the terminal and reopen a new terminal
create an environment named covmap-bsl by running the following commands using the yaml file provided i.e. covmap-bsl.yml
```
mamba env create --file covmap-bsl.yml
```
Tools required are listed under the analysis file
