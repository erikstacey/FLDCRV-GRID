      SUBROUTINE CCHISQ(MODEL, DATA, ERR, CHISQ, NPH)
c	Computes the chi squared given input arrays MODEL, DATA, ERR
c	and the integer lengths of the data stores in those arrays NPH
c	Outputs a float to Chisq
        REAL, DIMENSION(1000) :: MODEL, DATA, ERR
        INTEGER NPH
        REAL CHISQ
        CHISQ = 0.0
        DO I=1, NPH
            CHISQ = CHISQ + ((DATA(I)-MODEL(I))/(ERR(I)))**2
        END DO
        RETURN
      END