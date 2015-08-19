OBJECTS = Nucleus.o Mass.o Chart.o Yrast.o TlArray.o  LevelDensity.o Angle.o Nuclide.o LightP.o Evap.o AngleDist.o Random.o TlBarDist.o Scission.o  SigBarDist.o SigCharged.o Gdr.o Fus.o Run.o Weight.o RunThick.o RunAlpha.o AlphaOM.o Scatter.o Waves.o PotPara.o Merson.o Coul.o Avrigeanu.o

OBJECTS0 = Nucleus.o Mass.o Chart.o Yrast.o TlArray.o  LevelDensity.o Angle.o Nuclide.o LightP.o Evap.o AngleDist.o Random.o TlBarDist.o Scission.o Gdr.o

ALLOBJECTS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
FOBJECTS:=$(patsubst %.f,%.o,$(wildcard *.f))
CFLAGS= -c -Wall -W -O2 -I$(shell root-config --incdir)
COMPILER= c++ 
LINKOPTION = $(shell root-config --libs)

testDecay:  testDecay.o $(OBJECTS)
	$(COMPILER) -o testDecay testDecay.o $(OBJECTS) $(LINKOPTION) 

testFusion:  testFusion.o $(OBJECTS)
	$(COMPILER) -o testFusion testFusion.o $(OBJECTS) $(LINKOPTION) 

fusion:  fusion.o $(OBJECTS)
	$(COMPILER) -o fusion fusion.o $(OBJECTS) $(LINKOPTION) 

fusionThick:  fusionThick.o $(OBJECTS)
	$(COMPILER) -o fusionThick fusionThick.o $(OBJECTS) $(LINKOPTION) 

testGemini:: testGemini.o gemini.o $(OBJECTS0)
	f77 -o testGemini testGemini.o gemini.o $(OBJECTS0) -lstdc++

gemini.a: $(OBJECTS)
	ar rcs gemini.a $(OBJECTS)
	ranlib gemini.a

$(ALLOBJECTS): %.o : %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@

$(FOBJECTS): %.o : %.f
	f77 $(CFLAGS) $< -o $@


clean:
	rm -f *.o

