ROOTCONFIG   := root-config

ObjSuf        = o
SrcSuf        = cxx
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

CXX           = g++
CXXFLAGS      = -O -Wall -fPIC -g
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

DEPS = physics.h BaseGen.h BasePart.h ComptonGen.h Pi0PhotGen.h

.SUFFIXES: .$(SrcSuf)

# programs
TG  = EventGen$(ExeSuf)
TGO = EventGen.$(ObjSuf)
TGS = EventGen.$(SrcSuf)

OBJS	= $(TGO)

PROGRAMS = $(TG)

all:	$(PROGRAMS)

EventGen.$(SrcSuf): $(DEPS) Makefile
EventGen.$(ObjSuf): $(DEPS) Makefile

$(TG):	$(TGO)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(OutPutOpt)$@
		@echo "$@ done"

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<

clean:
		@echo "cleaning up..."
		@rm -f $(OBJS) core* out/*.root
