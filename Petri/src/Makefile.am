noinst_PROGRAMS = petri

if WITH_LIBDDD_PATH
LIBDDD_CPP=-I$(LIBDDD_ROOT)/include
LIBDDD_PATH=$(LIBDDD_ROOT)/lib/libDDD.a
else
LIBDDD_LINK=-lDDD
endif


if WITH_LIBSPOT_PATH
LIBSPOT_CPP=-I$(LIBSPOT_ROOT)/include
LIBSPOT_PATH=$(LIBSPOT_ROOT)/lib/libspot.a $(LIBSPOT_ROOT)/lib/libbddx.a
else
LIBSPOT_LINK=-lspot -lbddx
endif

if WITH_LIBEXPAT_PATH
LIBEXPAT_CPP=-I$(LIBEXPAT_ROOT)/include
LIBEXPAT_PATH=$(LIBEXPAT_ROOT)/lib/libexpat.a
else
LIBEXPAT_LINK=-lexpat
endif

AM_CPPFLAGS = $(LIBSPOT_CPP) $(LIBDDD_CPP)  -g -O2 -Wall $(LIBEXPAT_CPP)
AM_CXXFLAGS = -g -O2 -Wall

LDADDFLAGS = $(LIBSPOT_PATH) $(LIBSPOT_LINK) $(LIBDDD_PATH)  $(LIBEXPAT_PATH) $(LIBEXPAT_LINK) 

petri_SOURCES=Petri.cpp SparseIntArray.cpp SparseIntArray.h MatrixCol.h SparsePetriNet.h Walker.h parse/PTNetHandler.h parse/PTNetLoader.h parse/PTNetLoader.cpp expr/parse/ExprHandler.h expr/parse/ExprLoader.h expr/parse/ExprLoader.cpp
petri_LDADD = $(LDADDFLAGS)

