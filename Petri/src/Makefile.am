noinst_PROGRAMS = petri32 petri64 petri128 kersconv

AM_LDFLAGS = $(STATICFLAGS)

if WITH_LIBEXPAT_PATH
LIBEXPAT_CPP=-I$(LIBEXPAT_ROOT)/include
LIBEXPAT_PATH=$(LIBEXPAT_ROOT)/lib/libexpat.a
else
LIBEXPAT_LINK=-lexpat
endif

AM_CPPFLAGS = -g -DNDEBUG -O2 -Wnarrowing -Wextra -Wall -I$(srcdir) $(LIBEXPAT_CPP)
AM_CXXFLAGS = -g -O2 -DNDEBUG -Wall -Wconversion -Wnarrowing -Wextra 

LDADDFLAGS = $(LIBEXPAT_PATH) $(LIBEXPAT_LINK)

petri_SOURCES=Petri.cpp core/Arithmetic.hpp core/SparseArray.h core/SparseBoolArray.h core/MatrixCol.h core/SparsePetriNet.h core/Rational.h core/InvariantHelpers.h core/Log.h parse/PTNetHandler.h parse/PTNetLoader.h invariants/InvariantMiddle.h invariants/InvariantCalculator.h invariants/InvariantsTrivial.h invariants/RowSigns.h invariants/RowSignDomination.h invariants/MixedSignsUniqueTable.h invariants/Heuristic.h io/FlowPrinter.h io/SparseMatrixIO.h io/MatrixExporter.h io/PNMLExport.h walk/Walker.h expr/Expression.h expr/Simplify.h expr/Property.h parse/mcc/PropertyHandler.h parse/mcc/PropertyLoader.h

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

kersconv_SOURCES = kersconv.cpp io/SparseMatrixIO.h core/MatrixCol.h core/SparseArray.h
kersconv_CPPFLAGS = $(AM_CPPFLAGS) -DVAL=long
kersconv_CXXFLAGS = $(AM_CXXFLAGS)
kersconv_LDADD =


