# Pull from centos7 image
FROM centos:centos7


# copy needed folders and files into container
COPY environment_linux.yml /usr/local/etc/environment_linux.yml
COPY .container-files /usr/local/etc/container-files
RUN chmod 555 -R /usr/local/etc/container-files

# Install needed linux packages silently
RUN yum install -y wget which tar curl mesa-libGL-devel

# Install python 3.7 most recent stable miniconda version 4.9.2
RUN curl -sSLO https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh && \
	bash Miniconda3-py37_4.9.2-Linux-x86_64.sh -b -p /usr/local/miniconda && \
	rm Miniconda3-py37_4.9.2-Linux-x86_64.sh

# Add miniconda to path and set other environment variables
ENV PATH="/usr/local/miniconda/bin:$PATH" \
    CPATH="/usr/local/miniconda/include:$CPATH" \
    LANG="C.UTF-8" \
    LC_ALL="C.UTF-8" \
    PYTHONNOUSERSITE=1

# Install simnibs from conda based on: https://simnibs.github.io/simnibs/build/html/installation/conda.html
RUN conda env create -f /usr/local/etc/environment_linux.yml
SHELL ["conda", "run", "-n", "simnibs_env", "/bin/bash", "-c", "pip", "install", "-f", "https://github.com/simnibs/simnibs/releases/v3.2.1", "simnibs"]
SHELL ["conda", "run", "-n", "simnibs_env", "/bin/bash", "-c", "pip", "install", "-f", "nibabel", "pybids", "numpy","nipype","pandas"]
RUN mkdir /usr/local/SimNIBS 
SHELL ["conda", "run", "-n", "simnibs_env", "/bin/bash", "-c", "postinstall_simnibs", "--copy-matlab", "--setup-links", "-d", "/usr/local/SimNIBS"]

# Place simnibs linux tar.gz into image and install silently
#RUN cd /tmp && \
#	wget https://github.com/simnibs/simnibs/releases/download/v3.2.1/simnibs_installer_linux.tar.gz && \
#	tar -xzvf simnibs_install_linux.tar.gz && \
#	cd simnibs_installer && \
#	./install -t /usr/local -s

ENTRYPOINT ["conda","run","--no-capture-output","-n","simnibs_env","python","/usr/local/etc/container-files/run.py"]
	
 
