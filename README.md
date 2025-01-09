## Build & Install Gemini++ shared library & header files

This is a modified version of the Gemini++ source code taken from R. J. Charity's
[[website|https://wustl.app.box.com/s/ehsih9oc1j41loxgl4ox/folder/6560861989]] (see `upstream` branch for version details).
The modifications concern only the way the library is built and installed:

  * the original `Makefile` made a static library, but did not install it or the header files required for development;
  * we build a shared library, `libGemini++.so` and install it and all required files in a standard architecture:
  
      - [path to installation]/lib/gemini++/libGemini++.so
      - [path to installation]/lib/cmake/gemini++/ *CMake configuration files*
      - [path to installation]/include/gemini++/   *all header files*
      - [path to installation]/share/gemini++/    *tbl/ and tl/ directories*
  
To build and install:

  * clone this repository: 
  
      - git clone https://github.com/jdfrankland/gemini.git
      
  * outside the source directory, configure the build using `cmake` (version 3.12 or greater),
    giving the required installation directory:
  
      - cmake -S [path to sources] -B build -DCMAKE_INSTALL_PREFIX=[path to installation]
      
  * to build and install, do
  
      - cmake --build build --target install --parallel
      
  * to use the library in another CMake project, do the following:
  
~~~{.cmake}  
      find_package(Gemini++)
      add_executable(myExec myExec.cpp) # your code using the Gemini++ library
      target_link_libararies(myExec PUBLIC Gemini++)
~~~
  
  * there is no need to define/modify any environment variables
      
### Gemini++ References
R. J. Charity, "Systematic description of evaporation spectra for light and heavy compound nuclei", [[Physical Review C82, 014610 (2010)|https://doi.org/10.1103/PhysRevC.82.014610]]

D. Mancusi, R. J. Charity & J. Cugnon, "Unified description of fission in fusion and spallation reactions", [[Physical Review C82, 044610 (2010)|https://doi.org/10.1103/PhysRevC.82.044610]]
