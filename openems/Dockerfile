FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y autoremove
RUN apt-get update -y && apt-get upgrade -y
RUN apt-get install -y --fix-missing build-essential python3-pip cmake git libhdf5-dev libvtk7-dev libboost-all-dev libcgal-dev libtinyxml-dev qtbase5-dev libvtk7-qt-dev octave liboctave-dev parallel
RUN pip3 install matplotlib cython h5py setuptools numpy scipy
RUN git clone --recursive https://github.com/thliebig/openEMS-Project.git
WORKDIR "/openEMS-Project"
RUN bash -c './update_openEMS.sh /opt/openEMS --python'
RUN bash -c 'export PATH="/opt/openEMS:$PATH"'
