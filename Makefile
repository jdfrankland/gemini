OBJECTS = Nucleus.o Mass.o Chart.o Yrast.o TlArray.o  LevelDensity.o Angle.o Nuclide.o LightP.o Evap.o AngleDist.o Random.o TlBarDist.o Scission.o  Weight.o   SigCharged.o SigBarDist.o Fus.o Gdr.o


ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
FOBJECTS:=$(patsubst %.f,%.o,$(wildcard *.f))
CFLAGS= -c -Wall -W -O2
COMPILER= c++ 


testDecay:  testDecay.o $(OBJECTS)
	$(COMPILER) -o testDecay testDecay.o $(OBJECTS)  

testTheWidth:  testTheWidth.o $(OBJECTS)
	$(COMPILER) -o testTheWidth testTheWidth.o $(OBJECTS)  

testFusion:  testFusion.o $(OBJECTS)
	$(COMPILER) -o testFusion testFusion.o $(OBJECTS)  

testGemini:: testGemini.o gemini.o $(OBJECTS)
	f77 -o testGemini testGemini.o gemini.o $(OBJECTS) -lstdc++

testWidth:: testWidth.o gemini.o $(OBJECTS)
	f77 -o testWidth testWidth.o gemini.o $(OBJECTS) -lstdc++

gemini.a: $(OBJECTS)
	ar rcs gemini.a $(OBJECTS)
	ranlib gemini.a

$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

$(FOBJECTS): %.o : %.f
	f77 $(CFLAGS) $< -o $@


clean:
	rm -f *.o

