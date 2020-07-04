CC = g++-10 -std=c++20
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev` -lboost_program_options
CXXFLAGS = -lboost_program_options -Wall -g `pkg-config --cflags ibsimu-1.0.6dev` 


all: analysis.o simulation.o config.o config-setup.o beam.o
	$(CC) -o simulation config.o config-setup.o beam.o simulation.o $(LDFLAGS)
	$(CC) -o analysis config.o config-setup.o analysis.o $(LDFLAGS)
# alla fine 

config.o: config.cpp config.h  datatype.h
	$(CC) $(CXXFLAGS) -c -o config.o config.cpp

config-setup.o: config-setup.cpp config-setup.h beam.h datatype.h
	$(CC) $(CXXFLAGS) -c -o config-setup.o config-setup.cpp

beam.o: beam.cpp beam.h datatype.h
	$(CC) $(CXXFLAGS) -c -o beam.o beam.cpp

simulation.o: simulation.cpp datatype.h
	$(CC) $(CXXFLAGS) -c -o simulation.o simulation.cpp

analysis.o: analysis.cpp datatype.h
	$(CC) $(CXXFLAGS) -c -o analysis.o analysis.cpp

# ecr: ecr.o
# 	$(CC) -o ecr ecr.o $(LDFLAGS)
# ecr.o: ecr.cpp

# analysis: analysis.o
# 	$(CC) -o analysis analysis.o $(LDFLAGS)
# analysis.o: analysis.cpp

# clean:
# 	$(RM) *.o analysis ecr 

