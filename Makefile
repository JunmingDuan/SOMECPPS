INCDIR = .
HEADER = $(wildcard $(INCDIR)/*.h)
BIN = $(patsubst %.cpp, %, $(wildcard *.cpp))

CXXFLAGS = -I$(INCDIR) -Wall

LDFLAGS = -lm -ldl -O3 -march=native

all: $(BIN)
$(BIN): % : %.cpp $(HEADER)
	$(CXX) -o $@ $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS)
  

.PHONY: clean
clean:
	-rm $(BIN)
