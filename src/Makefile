CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -lboost_program_options -Wall -g `pkg-config --cflags ibsimu-1.0.6dev` 


all: analysis.o config.o simulation.o
	$(CC) -o simulation simulation.o config.o  $(LDFLAGS) -lboost_program_options
	$(CC) -o analysis analysis.o config.o  $(LDFLAGS) -lboost_program_options

config.o: config.cpp config.h
#	$(CC) -c config.cpp $(LDFLAGS) -lboost_program_options

simulation: config.o simulation.o

simulation.o: simulation.cpp config.h
#	$(CC) -c simulation.cpp $(LDFLAGS) -lboost_program_options

analysis: config.o analysis.o
	$(CC) -o analysis analysis.o config.o  $(LDFLAGS) -lboost_program_options

analysis.o: analysis.cpp config.h

# ecr: ecr.o
# 	$(CC) -o ecr ecr.o $(LDFLAGS)
# ecr.o: ecr.cpp

# analysis: analysis.o
# 	$(CC) -o analysis analysis.o $(LDFLAGS)
# analysis.o: analysis.cpp

# clean:
# 	$(RM) *.o analysis ecr 

