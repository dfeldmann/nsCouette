################################################################################
# This file is part of nsCouette -- A high-performance code for direct         #
# numerical simulations of turbulent Taylor-Couette flow                       #
#                                                                              #
# Copyright (C) 2019 Marc Avila, Bjoern Hof, Jose Manuel Lopez, Markus Rampp,  #
#                    Liang Shi, Alberto Vela-Martin, Daniel Feldmann.          #
#                                                                              #
# nsCouette is free software: you can redistribute it and/or modify it under   #
# the terms of the GNU General Public License as published by the Free         #
# Software Foundation, either version 3 of the License, or (at your option)    #
# any later version.                                                           #
#                                                                              #
# nsCouette is distributed in the hope that it will be useful, but WITHOUT ANY #
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS    #
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more        #
# details.                                                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# nsCouette. If not, see <http://www.gnu.org/licenses/>.                       #
################################################################################

#Makefile for the GPU version of nsCouette

NVCC = /opt/cuda/bin/nvcc   
LIBS = -lcufft -lcublas -lhdf5 -lhdf5_hl -lcurand -lcusparse
DEBUG = -g
GPU_SOURCES = $(wildcard *.cu)
GPU_OBJECTS = $(GPU_SOURCES:.cu=.o)



all: $(GPU_OBJECTS)
	$(NVCC) -o taylorC $(GPU_OBJECTS) $(LIBS)

$(GPU_OBJECTS): %.o: %.cu
	$(NVCC) -c   $< -o $@

clean:
	rm src/*.o taylorC
