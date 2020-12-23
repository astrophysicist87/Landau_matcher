# =======================================================================================
#  Makefile for hydrodynamic causality checks     Christopher Plumberg, December 23, 2020
# =======================================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			distclean	remove all binaries
##  

# Set compiler and flags
CC := g++
CFLAGS= -std=c++11 -lgsl -lgslcblas -lm

# Various directories and definitions
RM          =   rm -f
O           =   .o
LDFLAGS     =   
INCFLAGS    =   
SYSTEMFILES =   $(SRCGNU)

# --------------- Files involved ------------------

MAIN		 =	convert_IPGlasma
MAIN2		 =	convert_IPGlasma_for_MUSIC
MAIN3		 =	transpose_eps_u_pi

MAINSRC      =   convert_IPGlasma.cpp
MAIN2SRC	 =   convert_IPGlasma_for_MUSIC.cpp
MAIN3SRC	 =   transpose_eps_u_pi.cpp

INC		= 	

# -------------------------------------------------

TARGET		=	$(MAIN)
TARGET2		=	$(MAIN2)
TARGET3		=	$(MAIN3)

# --------------- Pattern rules -------------------

$(TARGET):
	$(CC) $(MAINSRC) -o $(TARGET) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

$(TARGET2):
	$(CC) $(MAIN2SRC) -o $(TARGET2) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

$(TARGET3):
	$(CC) $(MAIN3SRC) -o $(TARGET3) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

# -------------------------------------------------

.PHONY:		all help distclean

all:		$(TARGET) $(TARGET2) $(TARGET3)

help:
		@grep '^##' GNUmakefile

distclean:	
		-rm $(TARGET)
		-rm $(TARGET2)
		-rm $(TARGET3)

# --------------- Dependencies -------------------
convert_IPGlasma.cpp:         
convert_IPGlasma_for_MUSIC.cpp:         
transpose_eps_u_pi.cpp:         
