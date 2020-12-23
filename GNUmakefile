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

ifeq "$(MAIN)" ""
MAIN		 =	convert_IPGlasma
endif

ifeq "$(MAIN)" ""
MAIN2		 =	convert_IPGlasma_for_MUSIC
endif

MAINSRC      =   convert_IPGlasma.cpp
MAIN2SRC	 =   convert_IPGlasma_for_MUSIC.cpp

INC		= 	

# -------------------------------------------------

TARGET		=	$(MAIN)
TARGET2		=	$(MAIN2)

# --------------- Pattern rules -------------------

$(TARGET):
	$(CC) $(MAINSRC) -o $(TARGET) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

$(TARGET2):
	$(CC) $(MAINSRC2) -o $(TARGET2) $(CFLAGS) $(INCFLAGS)  $(LDFLAGS)

# -------------------------------------------------

.PHONY:		all help distclean

all:		mkobjdir $(TARGET) $(TARGET2)

help:
		@grep '^##' GNUmakefile

distclean:	
		-rm $(TARGET)
		-rm $(TARGET2)

# --------------- Dependencies -------------------
convert_IPGlasma.cpp:         
convert_IPGlasma_for_MUSIC.cpp:         
