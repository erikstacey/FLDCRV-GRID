      PROGRAM FLDCURV_GRID
      INTEGER :: NINCL, NBINCL, NPHADJ, NDP0S
      REAL :: CHISQ, BCHISQ, NITER
	REAL :: BESTI, BESTB, BESTPH, BESTDP
c     These parameters define the grid. UL=upper limit, LL=lowerlimit
c     N* = number of values
	REAL :: LLI, ULI, LLB, ULB, LLPH, ULPH, LLDP0, ULDP0
      PARAMETER(NINCL=20, NBINCL=20, NPHADJ=10, NDP0S=10)
      PARAMETER(LLI=0, LLB=0, LLPH=-0.5, LLDP0=500)
      PARAMETER(ULI=90, ULB=180, ULPH=0.5, ULDP0=8000)
      REAL, DIMENSION(1000) :: INPPH, CRVMDL, HJDS, GRBLS, SIGBLS, GRNUL
	REAL, DIMENSION(1000) :: SIGNUL
      REAL, DIMENSION(NINCL) :: INCLS
      REAL, DIMENSION(NBINCL) :: BINCLS
      REAL, DIMENSION(NPHADJ) :: PHADJS
      REAL, DIMENSION(NDP0S) :: DP0S
c     NPH is the number of phases to compute over, not be confused with the phase adjustments
      INTEGER NPH
c	Timing
	REAL :: TSTART, TEND
      WRITE(6, *) "[FLDCURV] Setting up grid..."
c     Set up grid to compute over
      DO I=1, NINCL
            INCLS(I) = LLI + (I-1) * (ULI-LLI)/NINCL
      END DO
      DO I=1, NBINCL
            BINCLS(I) = LLB + (I-1) * (ULB-LLB)/NBINCL
      END DO

      DO I=1, NPHADJ
            PHADJS(I) = LLPH + (I-1) * (ULPH-LLPH) / NPHADJ
      END DO

      DO I=1, NDP0S
            DP0S(I) = LLDP0 + (I-1) * (ULDP0-LLDP0) / NDP0S
      END DO

      OPEN(2,FILE='./fldcurv.out',STATUS='unknown')
      OPEN(3,FILE='./grunhut_bls.txt',STATUS='OLD')
      

c	Set fixed pars
	BQ=0
	BOCT=0
	A=0
	EPS=0.3

c	Set Best chisq to 0, if bchisq is 0 always accept challenger in loop
	BCHISQ = 0
	BESTI=0
	BESTB=0
	BESTPH=0
	BESTDP=0

	WRITE(6, *) "[FLDCURV] Reading data from file..."
c	Read in data from txt file
      READ(3, *) NPH
	DO I=1, NPH
		READ(3, *) HJDS(I), INPPH(I), GRBLS(I), SIGBLS(I), GRNUL(I), SIGNUL(I)
	ENDDO

	WRITE(6, *) "[FLDCURV] Computing for ", NPH, " phases over "
	WRITE(6, *) NINCL*NBINCL*NPHADJ*NDP0S, " iterations"
	CALL CPU_TIME(TSTART)
c     Iterate over grid, compute
	NITER = 0
      DO IDXI=1, NINCL
      DO IDXB=1, NBINCL
      DO IDXPH=1, NPHADJ
      DO IDXDP=1, NDP0S
		ITER = NITER+1
c		WRITE(6, *) "[FLDCURV] ITERATION ", NITER  
c		WRITE(6, *) "	I=", INCLS(IDXI)
c		WRITE(6, *) "	B=", BINCLS(IDXB)
c		WRITE(6, *) "	DP0=", DP0S(IDXDP)
c		WRITE(6, *) "	PHADJ=", PHADJS(IDXPH)

c		Compute the Bl over all the phases
		DO I=1,NPH
c			Phase corrected by grid parameter - FLDCURV computes
c			according to a fixed rotational phase, but we want the
c			phase to be a free parameter. Therefore, we apply an
c			offset to the phases that FLDCURV sees and then just
c			ignore them otherwise
			PHASE=INPPH(I)-PHADJS(IDXPH)
c			UNFLD written by J.D. Landstreet. Computes Bl and Bs given
c			input parameters I,B, Rotational phase, dipolar field
c			strength, decentering, quadropole strength, octopole
c			strength, limb darkening
			CALL UNFLD(INCLS(IDXI),BINCLS(IDXB),PHASE,DP0S(IDXDP),
     &		A,BQ,BOCT,EPS,BL,BS)
			R = BL/BS
c			Current model is stored in CRVMDL
			CRVMDL(I) = BL
		ENDDO
c		Chi Squared is used to measure goodness of fit
c		Best fit parameters are stored in BESTI, BESTB,
c		BESTPH, BESTDP, and BCHISQ.
		CALL CCHISQ(CRVMDL, GRBLS, SIGBLS, CHISQ, NPH)
		IF ((BCHISQ .GT. CHISQ) .OR. (BCHISQ .EQ. 0)) THEN
			WRITE(6, *) "[FLDCURV] Found improved fit..."
			WRITE(6, *) "	Previous best chisq: ", BCHISQ
			WRITE(6, *) "	New best chisq: ", CHISQ
			
			BCHISQ = CHISQ
			BESTI = INCLS(IDXI)
			BESTB = BINCLS(IDXB)
			BESTPH = PHADJS(IDXPH)
			BESTDP = DP0S(IDXDP)
			WRITE(6, *) "	New best pars:", BESTI, 
     &		BESTB, BESTPH, BESTDP
		END IF 

      ENDDO
      ENDDO
      ENDDO
      ENDDO
	CALL CPU_TIME(TEND)

	WRITE(6, *) "[FLDCURV] Grid computation complete in ", TEND-TSTART
	WRITE(6, *) "[FLDCURV] BEST FIT:"
	WRITE(6, *) "	I=", BESTI
	WRITE(6, *) "	B=", BESTB
	WRITE(6, *) "	DP0=", BESTDP
	WRITE(6, *) "	PH_ADJ=", BESTPH

	WRITE(6, *)"[FLDCURV] Writing best fit model to fldcurv.out"

c	Recompute best fit model, then write it to fldcurv.out
	DO I=1, NPH
		CALL UNFLD(BESTI,BESTB,INPPH(I) - BESTPH,BESTDP,
     &		A,BQ,BOCT,EPS,BL,BS)
		WRITE(2, *) INPPH(I), BL
	ENDDO
	END PROGRAM

