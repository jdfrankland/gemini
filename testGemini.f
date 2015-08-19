      implicit none
 
      real*4 fEx,fJ,vx,vy,vz
      real*4 thetaJ,phiJ
      integer*4 iz,ia,n,gemini,i,j

c     include common bock for 
      include "gemini.inc"

      call geminiinit()
      vx = 0.
      vy = 0.
      vz = 0.
      iZ=70
      iA=160
      fEx = 100.
      fJ = 40.
      thetaJ = 90.
      phiJ = 90.
      n = gemini(iz,ia,fex,fj,thetaJ,phiJ,vx,vy,vz)

c     if n=0, then gemini was not able to decay the event
      IF (n .GT. 0) THEN
         write (*,*) "number of final particles=",n
         DO i=1,n
           write(*,*) z(i),a(i),T(i),P(i),theta(i),phi(i)
         END DO
      END IF 

      do j=1,10
         n = gemini(iz,ia,fex,fj,thetaJ,phiJ,vx,vy,vz)
         IF (n .GT.  0) THEN
            write (*,*) "number of final particles=",n
            DO i=1,n
              write(*,*) z(i),a(i),T(i),P(i),theta(i),phi(i)
            END DO
        END IF 
      END DO
      end

      
