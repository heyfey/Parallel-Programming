NVFLAGS  := -std=c++11 -O3 -Xptxas="-v" -arch=sm_61 
LDFLAGS  := -lm
EXES     := seq hw4-2
CXXFLAGS := -fopenmp

alls: $(EXES)

clean:
	rm -f $(EXES)

seq: seq.cc
	g++ $(CXXFLAGS) -o $@ $?

hw4-2: hw4-2.cu
	nvcc $(NVFLAGS) $(LDFLAGS) -Xcompiler="$(CXXFLAGS)" -o $@ $?
