C     test program to show how to get decay width 
c     include light-particle evaporation and E1 and E2 gamma emission.
      implicit none


      real*4 fex,fj
      real*4 width
      integer*4 iz,ia
      real*4 logleveldensity
      real*4 leveldensity

c     include common bock for 
      include "gemini.inc"

      call geminiinit()

      iz=12 !compound nucleus Z
      ia=24 !compound nucleus A
      fex = 14. !excitation energy of compound nucleus in MeV
      fj = 4. !angular momentum of compound nucleus

      call decaywidth(iz,ia,fex,fj,width,logleveldensity)
      leveldensity = exp(logleveldensity);
      Write(*,*) "gamma= " ,width," MeV, rho= ", leveldensity," MeV-1"
      Write(*,*) "gamma*rho= ", width*leveldensity 
      end

      
