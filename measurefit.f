      SUBROUTINE CCHISQ(MODEL, DATA, ERR, CHISQ, NPH)
        REAL, DIMENSION(1000) :: MODEL, DATA, ERR
        INTEGER NPH
        REAL CHISQ
        CHISQ = 0.0
        DO I=1, NPH
            CHISQ = CHISQ + ((DATA(I)-MODEL(I))/(ERR(I)))**2
        END DO
        RETURN
      END