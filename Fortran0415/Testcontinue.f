***
      PROGRAM Testcontinue
***
*
* Evolves a population of binaries using input parameters 
* read from input file binaries.in (M1, M2, P, e, Z, Tmax). 
*
***
      implicit none
      INTEGER i, j, k
        do i = 1,3
            do j=4,6
                do k=7,9
                    
                    if(k.eq.8)THEN
                        CONTINUE
                    ENDIF
                    print*,i,j,k
                enddo
            ENDDO
        ENDDO
      STOP
      END
***

