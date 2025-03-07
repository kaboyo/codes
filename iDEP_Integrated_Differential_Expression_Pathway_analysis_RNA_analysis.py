#https://github.com/gexijin/idepGolem?tab=readme-ov-file


```
#https://youtu.be/Hs5SamHHG9s
#Integrated Differential Expression & Pathway analysis (iDEP)
#https://github.com/gexijin/idepGolem?tab=readme-ov-file
#https://idepsite.wordpress.com/
#https://thinkr-open.github.io/golem/
#https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2486-6
#https://bioinformatics.sdstate.edu/idep/ : Online platform

#https://palomerolab.org/how-to/workflows/RNAseq/

#TUXEDO
#https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/Apps/TuxedoSuite_RNASeqTools.htm

#associated papers:
#https://www.nature.com/articles/nbt.2450#accession-codes


#other tools for RNASeq
-https://salmon.readthedocs.io/en/latest/salmon.html ; https://combine-lab.github.io/salmon/ 
-https://pmc.ncbi.nlm.nih.gov/articles/PMC9558114/: High-quality reference transcriptome construction improves RNA-seq quantification in Oryza sativa indica

```


*#Integrated Differential Expression & Pathway analysis (iDEP)*

#Original repo

#Description

#iDEP is a bioinformatics platform for analyzing gene expression data graphically, interactively, and reproducibly. 
# The input file is a gene-level expression matrix derived from RNA-Seq, microarray, proteomics, or other methods. 
# Hosted at South Dakota State Univerisity, iDEP is developed as an R package based on the Golem framework, by a small team 
# led by Dr. Steven Ge. See documentation and paper.

#License:(CC BY-NC 3.0) Non-commercial use. For local installation at private institutions, please contact us.

#Run iDEP locally on Windows, MacOS, or Linux***

#Following the redesign of our database, the updated iDEP can now be effortlessly executed on a 
# laptop or local server. Please be aware that iDEP undergoes frequent updates; if you are running the software locally, 
# we recommend updating to the most recent version on a monthly basis.

#Once installed as an R package, iDEP can be initiated from R using the iDEP::run_app() command. 
# However, this approach can be time-consuming, as it necessitates the installation of all 355 dependent R packages. 
# This process takes approximately an hour on Windows 10, and potentially even longer for Mac or Linux users due to potential 
# troubleshooting requirements for many of the R packages.


#A more efficient alternative is to utilize the iDEP Docker image available on DockerHub. 
# The only prerequisite is to install the Docker software on Linux or Docker Desktop on Windows or MacOS.

Please note that local installation is free for non-profit organizations only.

#System requirements
#Most of modern laptop can run iDEP locally, if it has more than 10GB storage and 4GB of memory.

Windows: Docker Desktop (~20 minutes, recommended)
Just follow this detailed video. No prior experience with Docker is needed.

#Install Docker Desktop.
Start Windows PowerShell as an Administrator. From Windows search bar, type PowerShell to find the Windows PowerShell app. And then select Run as an Administrator. For details, see here.
Enable Windows Subsystem for Linux 2 (WSL2).
Start the Docker app. From Windows search bar, type Docker, and then select Run as an Administrator. Accept the terms when asked. The Docker engine is now running in the background.
Pull the iDEP Docker image from DockerHub and start a container from PowerShell.
docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest 
You can now use iDEP locally by starting your web browser and enter localhost:3838 in the address bar.
Note that the Docker engine is now running in the backgroup, acting as a webserver. It works even if you restart your computer. To stop it, run these from Windows PowerShell:

docker stop idep 
docker rm idep
After stopping it, you can restart it by repeating Step 5, which also pulls the latest Docker image from DockerHub. Make sure you update your image at least on a monthly basis.

Windows: iDEP as an R package (~1 hour)
Install a recent version of R. If your R is a few years old, uninstall it, and manually delete the folder that contains all existing R packages.
Optional: Install IDE such as RStudio or VS Code.
Start R and run these commands. It takes about an hour to install the 355 dependencies! If there are issues with the installation of some of the individual packages, you need to resolve the issues and try to install them.
install.packages("devtools")
devtools::install_github("https://github.com/gexijin/idepGolem", upgrade = "never")
You might get the following warnings, which can be ignored.

WARNING: Rtools is required to build R packages, but is not currently installed.
package ‘KEGG.db’ is not available for this version of R
A version of this package for your version of R might be available elsewhere,
2: In i.p(...) : installation of package ‘GO.db’ had non-zero exit status
3: In i.p(...) :
installation of package ‘/tmp/RtmpCpiUsZ/file2326c2a50656f/PGSEA_1.60.0.tar.gz’ had non-zero exit status
Start iDEP from R command console.
idepGolem::run_app()
The benefit of this approach is that you always have the most recent version from GitHub. The next time you install iDEP, it will take much less time as the dependencies have been installed.

Windows: Developer mode (~1 hour)
With this method, users can customize iDEP by changing the source code. We also welcome contributions through Pull Requests on GitHub.

Install a recent version of R, such as R 4.30.
Optional: Install IDE such as RStudio or VS Code.
Obtain a copy of the source code. This can be done manually by downloading an zip file by clicking on the green Code button above. Alternatively, you can fork this repository, install GitHub Desktop, and clone this repository locally using an URL(https://github.com/gexijin/idepGolem.git).
Start R and install the golem package.
install.packages("golem")
Optional: Install all dependencies. We pretend to install the idepGolem package. This take about an hour.
install.packages("devtools")
devtools::install_github("https://github.com/gexijin/idepGolem", upgrade = "never")
remove.packages("idepGolem")
Open the downloaded reporsitory from R. The Shiny app can be started by running the run_dev.R script. If some R packages are missing, install them and try again. Note that the app is devided into 11 Shiny Modules.
MacOS: Docker container (~10 minutes on MacBook Air)
See video for more details. No prior experience with Docker is needed.

Download Docker Desktop follow instructions here. Make sure you choose the correct version based on your CPU type. For MacBook Air, I chose Apple Chip. If you use an Mac computer with Intel CPU, choose the Intel Chip instead.
Install Docker Desktop. First double-click the downloaded Docker.dmg file in the Downloads folder. Drag the Docker icon into the Application folder from the pop-up window. From Lunch pad, or the Application folder, click on the Docker icon. Click Open when asked. Accept the Terms and the Docker engine is running.
Start a Terminal window by clicking the Launchpad, and type terminal in the search field. Then click the Terminal app.
Pull the iDEP Docker image from DockerHub and start a container.
docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest 
You can now use iDEP locally by starting your web browser and enter localhost:3838 in the address bar.
Note that the Docker engine is now running in the backgroup, acting as a webserver. It works even if you restart your computer. To stop it, run these from the Terminal app:

docker stop idep 
docker rm idep
After stopping it, you can restart it by repeating Step 4, which also pulls the latest Docker image from DockerHub. Make sure you upgrade your iDEP image at least on a monthly basis.

Alternatively, you can also install iDEP as an R package or copy the iDEP code locally. The method is the same as above in the Windows section.

Linux: Docker container (~10 minutes)
Install Docker Engine follow instructions here. Many Linux systems have Docker installed by default. On Ubuntu, I used the following scripts.
curl -fsSL https://get.docker.com -o install_docker.sh
sudo sh install_docker.sh
Pull the iDEP Docker image from DockerHub and start a container.
sudo docker run --pull always -d --name idep -p 3838:3838 gexijin/idep:latest 
You can now use iDEP locally by starting your web browser and enter localhost:3838 in the address bar. If this is an server, make sure that the port 3838 is open. Then replace the localhost with the IP address.
Note that the Docker engine is now running in the backgroup, acting as a webserver. To stop it:

sudo docker stop idep 
sudo docker rm idep
After stopping it, you can restart it by repeating Step 2, which also pulls the latest iDEP image from DockerHub. We update it frequently, make sure you upgrade your image at least on a monthly basis.

