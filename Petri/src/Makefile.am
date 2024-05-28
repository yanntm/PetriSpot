noinst_PROGRAMS = petri32 petri64 petri128

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

petri_SOURCES=Arithmetic.hpp Petri.cpp SparseArray.h SparseIntArray.h MatrixCol.h SparsePetriNet.h Walker.h PTNetHandler.h PTNetLoader.h SparseBoolArray.cpp SparseBoolArray.h InvariantMiddle.h InvariantCalculator.h FlowMatrix.h

# Rules for each program
petri32_SOURCES = $(petri_SOURCES)
petri32_CPPFLAGS = $(AM_CPPFLAGS) -DVAL=int
petri32_CXXFLAGS = $(AM_CXXFLAGS)
petri32_LDADD = $(LDADDFLAGS)

petri64_SOURCES = $(petri_SOURCES)
petri64_CPPFLAGS = $(AM_CPPFLAGS) -DVAL=long
petri64_CXXFLAGS = $(AM_CXXFLAGS)
petri64_LDADD = $(LDADDFLAGS)

petri128_SOURCES = $(petri_SOURCES)
petri128_CPPFLAGS = $(AM_CPPFLAGS) -DVAL="long long"
petri128_CXXFLAGS = $(AM_CXXFLAGS)
petri128_LDADD = $(LDADDFLAGS)


