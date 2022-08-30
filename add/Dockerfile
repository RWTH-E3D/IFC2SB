FROM ubuntu:20.04 # ubuntu:latest

WORKDIR /bim2sim

RUN apt-get update
RUN DEBIAN_FRONTEND="noninteractive" apt-get -y install tzdata
RUN apt-get -y install unzip


################################################################################
## Dependencies
RUN apt-get install -y git cmake gcc g++ tcllib tklib tcl-dev tk-dev libfreetype6-dev libxt-dev libxmu-dev libxi-dev libgl1-mesa-dev libglu1-mesa-dev libfreeimage-dev libtbb-dev qt5-default libboost-all-dev libxml2-dev libomp-dev
################################################################################


################################################################################
# OpenCascade
WORKDIR /tmp
#RUN git clone https://github.com/Open-Cascade-SAS/OCCT && cd OCCT && mkdir build && cd build 
#WORKDIR /tmp/OCCT/build
RUN git clone -b master https://github.com/bimworld/occ.git && cd occ && unzip opencascade-7.5.3.zip -d /tmp/
RUN cd opencascade-7.5.3 && mkdir build && cd build 
WORKDIR /tmp/opencascade-7.5.3/build
RUN cmake .. -Wno-dev -DBUILD_DOC_Overview=OFF -DBUILD_Inspector=OFF -DBUILD_LIBRARY_TYPE=Shared -DBUILD_RESOURCES=OFF -DBUILD_SAMPLES_QT=OFF -DBUILD_USE_PCH=OFF -DBUILD_YACCLEX=OFF -DCMAKE_BUILD_TYPE=Release -DUSE_TBB=ON -DUSE_FREEIMAGE=ON
RUN make -j16 
RUN make install
################################################################################


################################################################################
# IfcOpenShell
WORKDIR /tmp
RUN git clone --branch v0.6.0 https://github.com/IfcOpenShell/IfcOpenShell.git  && cd IfcOpenShell && mkdir build && cd build
WORKDIR /tmp/IfcOpenShell/build
RUN  cmake ../cmake -Wno-dev -DOCC_INCLUDE_DIR=/usr/local/include/opencascade -DBUILD_CONVERT=OFF -DBUILD_DOCUMENTATION=OFF -DBUILD_EXAMPLES=OFF -DBUILD_GEOMSERVER=OFF -DBUILD_CONVERT=OFF -DBUILD_PACKAGE=OFF -DCOLLADA_SUPPORT=OFF -DGLTF_SUPPORT=OFF -DBUILD_IFCPYTHON=OFF -DCMAKE_BUILD_TYPE=Release -DENABLE_BUILD_OPTIMIZATIONS=ON -DIFCXML_SUPPORT=OFF -DUSE_MMAP=OFF -DUSE_VLD=OFF -DBUILD_SHARED_LIBS=ON
RUN make -j16
RUN make install
################################################################################


################################################################################
# Clipper Library
RUN apt-get -y install libpolyclipping-dev
################################################################################


################################################################################
ADD ifc2sb /tmp/ifc2sb/code
RUN cd code && mkdir build
WORKDIR /tmp/ifc2sb/code/build
RUN cmake ../cmake -DCLIPPER_INCLUDE_PATH="/usr/include/polyclipping"
RUN make -j16
RUN ./IFC2SB -h
################################################################################


################################################################################
# Clean
WORKDIR /tmp
RUN rm -r /tmp/IfcOpenShell
RUN rm -r /tmp/occ
RUN rm -r /tmp/opencascade-7.5.3
# RUN rm -r /tmp/triangle
################################################################################
