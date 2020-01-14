CC=c++

#this one if for min-gw (Windows)
#CFLAGS=-I. -Irandomc -Istocc -IEasyBMP -Iconfig-src -enable-auto-import -static-libstdc++ -static-libgcc -Wall -g -O0  

#this one is for gcc (unix)
CFLAGS=-I. -Irandomc -Istocc -IEasyBMP -Iconfig-src -Wall -g -O0  



DEPS= eDiffusion.h
OBJ=eDiffusion.obj 


#%.o: %.cpp $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)


eDiffusion: eDiffusion.cpp ClassMonteCarlo.cpp eDiffusion.h
	$(CC) -o eDiffusion.exe eDiffusion.cpp ClassMonteCarlo.cpp ClassExciton.cpp ClassMedium.cpp ClassQuencher.cpp ClassBool3D.cpp randomc/mersenne.cpp stocc/stoc1.cpp stocc/userintf.cpp EasyBMP/EasyBMP.cpp config-src/config.cpp config-src/log.cpp $(CFLAGS) 
