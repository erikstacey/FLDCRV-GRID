      PROGRAM FLDCURV_GRID
c	Written by Erik Stacey, 30 Jan 2023. Based on FLDCRV program, originally
c	written by J.D. Landstreet, 1986. Utilizes UNFLD and MAGFLD subroutines
c	of original program.

c	A program for computing grids of synthetic longitudinal field curves
c	and comparing them against observed values to determine the best
c	fit physical parameters for a magnetic star.
c	As presently configured, this program holds fixed the following:
c	Limb Darkening: 0.3
c	Quadropolar component strength: 0
c	Octopolar component strength: 0
c	Decentering: 0

c	Requires following file in execution directory:
c	"bls.dat": input observed bls. Format HJD, Phase, Bl, SigmaBl, Null, SigmaNull
c	Outputs:
c	"fldcurv.out": Best fit synthetic longitudinal field curve, format: phase, Bl
c	"fldcurv.log": Diagnostic messages
c	"fldcurv_bestpars.out": Best pars as identified through the grid search

      INTEGER :: NINCL, NBINCL, NPHADJ, NDP0S
      REAL :: CHISQ, BCHISQ
	REAL :: BESTI, BESTB, BESTPH, BESTDP
c     These parameters define the grid. UL=upper limit, LL=lowerlimit
c     N* = number of values
	REAL :: LLI, ULI, LLB, ULB, LLPH, ULPH, LLDP0, ULDP0
      PARAMETER(NINCL=100, NBINCL=100, NPHADJ=10, NDP0S=40)
      PARAMETER(LLI=70, LLB=55, LLPH=-0.1, LLDP0=1000)
      PARAMETER(ULI=90, ULB=80, ULPH=0.1, ULDP0=3000)
      REAL, DIMENSION(1000) :: INPPH, CRVMDL, HJDS, GRBLS, SIGBLS, GRNUL
	REAL, DIMENSION(1000) :: SIGNUL
      REAL, DIMENSION(NINCL) :: INCLS
      REAL, DIMENSION(NBINCL) :: BINCLS
      REAL, DIMENSION(NPHADJ) :: PHADJS
      REAL, DIMENSION(NDP0S) :: DP0S
c     NPH is the number of phases to compute over, not be confused with the phase adjustments
c	LSTPCT stores the last reported percentage value for progress, as an integer from 1-100
c	if the measured percentage doesn't exceed this, it doesn't get reported
c	As these are determined through integer division, only integers between 1-100 can be reported
      INTEGER NPH, NITER, TITER, LSTPCTad
c	Timing
	REAL :: TSTART, TEND, CTIME
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
      OPEN(3,FILE='./bls.dat',STATUS='OLD')
	OPEN(11, FILE='./fldcurv.log', STATUS="unknown")
	OPEN(12, FILE='./fldcurv_bestpars.out', STATUS="unknown")
      

c	Set fixed pars
	BQ=0
	BOCT=0
	A=0
	EPS=0.3

c	Set Best chisq to 0, if bchisq is 0 always accept challenging chisq in loop
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

	TITER = NINCL*NBINCL*NPHADJ*NDP0S
c	Record the boundaries and steps to the log file
	WRITE(11, *) "Incl LB, UB, N: ", LLI, ULI, NINCL
	WRITE(11, *) "Beta Incl LB, UB, N:", LLB, ULB, NBINCL
	WRITE(11, *) "Phase corr LB, UB, N:", LLPH, ULPH, NPHADJ
	WRITE(11, *) "Dipolar FS LB, UB, N:", LLDP0, ULDP0, NDP0S

	WRITE(6, *) "[FLDCURV] Computing for ", NPH, " phases over "
	WRITE(6, *) TITER, " iterations"
	WRITE(6, *) "[FLDCURV] Estimated runtime: ",
     & FLOAT(TITER)*FLOAT(NPH)*0.000015, " s"

     	WRITE(11, *) "[FLDCURV] Computing for ", NPH, " phases over "
     	WRITE(11, *) TITER, " iterations"
     	WRITE(11, *) "[FLDCURV] Estimated runtime: ",
     &FLOAT(TITER)*FLOAT(NPH)*0.000015, " s"



	CALL CPU_TIME(TSTART)
c     Iterate over grid, compute
	NITER = 0
	LSTPCT = 0


      DO IDXI=1, NINCL
      DO IDXB=1, NBINCL
      DO IDXPH=1, NPHADJ
      DO IDXDP=1, NDP0S
		
c		Check for progress
c		overhead is about 2.5% in computation time
		if (NITER*100 / TITER .gt. LSTPCT) THEN
			LSTPCT = NITER*100 / TITER
			CALL CPU_TIME(CTIME)
			WRITE(6, *) "	Completed ", LSTPCT, "% of iterations (T=",
     &		CTIME-TSTART, " s)"
     			WRITE(11, *) "	Completed ", LSTPCT, "% of iterations (T=",
     &		CTIME-TSTART, " s)"
		END IF
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
			WRITE(11, *) "[FLDCURV] Found improved fit..."
			WRITE(11, *) "	Previous best chisq: ", BCHISQ
			WRITE(11, *) "	New best chisq: ", CHISQ
			
			BCHISQ = CHISQ
			BESTI = INCLS(IDXI)
			BESTB = BINCLS(IDXB)
			BESTPH = PHADJS(IDXPH)
			BESTDP = DP0S(IDXDP)
			WRITE(11, *) "	New best pars:", BESTI, 
     &		BESTB, BESTPH, BESTDP
		END IF
		NITER = NITER+1

      ENDDO
      ENDDO
      ENDDO
      ENDDO
	CALL CPU_TIME(TEND)

	WRITE(6, *) "[FLDCURV] Grid computation complete in ", TEND-TSTART, " s"
	WRITE(11, *) "[FLDCURV] Grid computation complete in ", TEND-TSTART, " s"
	WRITE(11, *) "[FLDCURV] BEST FIT:"
	WRITE(11, *) "	I=", BESTI
	WRITE(11, *) "	B=", BESTB
	WRITE(11, *) "	DP0=", BESTDP
	WRITE(11, *) "	PH_ADJ=", BESTPH

	WRITE(12, *) "I: ", BESTI
	WRITE(12, *) "B: ", BESTB
	WRITE(12, *) "DP0: ", BESTDP
	WRITE(12, *) "PHADJ: ", BESTPH
	WRITE(12, *) "CHISQ: ", BCHISQ

	WRITE(6, *)"[FLDCURV] Writing best fit model to fldcurv.out"

c	Recompute best fit model, then write it to fldcurv.out
	DO I=1, NPH
		CALL UNFLD(BESTI,BESTB,INPPH(I) - BESTPH,BESTDP,
     &		A,BQ,BOCT,EPS,BL,BS)
		WRITE(2, *) INPPH(I), BL
	ENDDO
	END PROGRAM

