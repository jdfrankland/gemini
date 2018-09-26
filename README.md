## Build & Install Gemini++ shared library & header files

This is a modified version of the Gemini++ source code taken from R. J. Charity's
[[website|https://wustl.app.box.com/s/ehsih9oc1j41loxgl4ox/folder/6560861989]] (see `upstream` branch for version details).
The modifications concern only the way the library is built and installed:

  * the original `Makefile` made a static library, but did not install it or the header files required for development;
  * we build a shared library, `libGemini.so` and install it and all required files in a standard architecture
  
To build and install:

  * clone this repository: 
  
      git clone https://github.com/jdfrankland/gemini
      
  * in the source directory, run `cmake` (version 2.8 or greater) giving the required installation directory:
  
      cmake . -DCMAKE_INSTALL_PREFIX=[path to installation]
      
  * to build and install, do
  
      make [-j x] install
      
    where `x` is the optional number of CPUs to use for a parallel build, if desired.
    
### Gemini++ References
R. J. Charity, "Systematic description of evaporation spectra for light and heavy compound nuclei", [[Physical Review C82, 014610 (2010)|https://doi.org/10.1103/PhysRevC.82.014610]]

D. Mancusi, R. J. Charity & J. Cugnon, "Unified description of fission in fusion and spallation reactions", [[Physical Review C82, 044610 (2010)|https://doi.org/10.1103/PhysRevC.82.044610]]
