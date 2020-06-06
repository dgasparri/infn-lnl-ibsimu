CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -lboost_program_options -Wall -g `pkg-config --cflags ibsimu-1.0.6dev` 

config: config.o simulation.o
	$(CC) -o simulation simulation.o config.o  $(LDFLAGS) -lboost_program_options

config.o: config.cpp config.h
#	$(CC) -c config.cpp $(LDFLAGS) -lboost_program_options

simulation.o: simulation.cpp config.h
#	$(CC) -c simulation.cpp $(LDFLAGS) -lboost_program_options




# ecr: ecr.o
# 	$(CC) -o ecr ecr.o $(LDFLAGS)
# ecr.o: ecr.cpp

# analysis: analysis.o
# 	$(CC) -o analysis analysis.o $(LDFLAGS)
# analysis.o: analysis.cpp

# clean:
# 	$(RM) *.o analysis ecr 

