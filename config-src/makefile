CC=g++

#this one if for min-gw (Windows)
CFLAGS=-I.  -enable-auto-import -Wall -g -O0  

#this one is for gcc (unix)
#CFLAGS=-I. -Irandomc  -Wall -g -O0  
DEPS= 
OBJ=


#%.o: %.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)


eDiffusion: main.cpp config.cpp config.h log.h  
	$(CC) -o main.exe main.cpp config.cpp config.h log.h log.cpp $(CFLAGS) 
