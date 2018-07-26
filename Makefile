PATH_TO_FASTJET = /Users/prekshanaik/Downloads/fastjet-install/bin/fastjet-config
PATH_TO_BOOST = /usr/lib/
PATH_TO_PYTHIA = /Users/prekshanaik/Documents/MEng_Project/pythia8235/bin/pythia8-config


CXX = g++
CXXFLAGS= -O3 -Wall -Woverloaded-virtual -g -std=c++11 -stdlib=libc++ -m64 

FASTINC = `$(PATH_TO_FASTJET) --cxxflags`
FASTLIB = `$(PATH_TO_FASTJET) --libs --plugins` -lRecursiveTools

BOOSTINC = `$(PATH_TO_BOOST)`
BOOSTLIB = `$(PATH_TO_BOOST) -lboost_filesystem -lboost_regex`

ROOTINC = `root-config --cflags --glibs`


PYTHIAINC = `$(PATH_TO_PYTHIA) --cxxflags`
PYTHIALIB = `$(PATH_TO_PYTHIA) --libs`


OBJDIR=src
EXECDIR=exec
BINDIR=bin
INCDIR=interface
INC= -I$(INCDIR)

_OBJ = InfoCalibratedJet InfoPFC Event Trigger Property Condition helpers
OBJ  = $(patsubst %,$(OBJDIR)/%,$(_OBJ:=.o))



# _EXEC=skim analyze turn_on convert_to_pristine analyze_data write
_EXEC=analyze_pfc analyze analyze_triggers
EXEC=$(patsubst %,$(EXECDIR)/%,$(_EXEC:=.o))
BIN=$(patsubst %,$(BINDIR)/%,$(_EXEC))


all: $(BIN)

$(OBJDIR)/%.o : $(OBJDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC) $(PYTHIAINC)  

$(EXECDIR)/%.o : $(EXECDIR)/%.cc
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(INC) $(FASTINC) $(PYTHIAINC)
	
$(BINDIR)/% : $(EXECDIR)/%.o $(OBJ)
	$(CXX) $< $(OBJ) -o $@ $(CXXFLAGS) $(FASTLIB) $(PYTHIALIB)

.PHONY: clean
.PRECIOUS: $(OBJ) $(EXEC)

clean:
	rm $(OBJ) $(EXEC) $(BIN)


