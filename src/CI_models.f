Cd
C=================================================================
C
      Subroutine ModImpose(Eta)
C
C Impose parameter relations in given model


      Implicit None
      include 'CI_models.inc'
      include 'steering.inc'
C       Real*8 par(8)
C       REAL, DIMENSION(4,6) :: Eta

C       Real :: Eta(4,2)
C
C  Models are:  
C  =========== 
C
C  One parameter (as used by ZEUS - David Williams talk at DIS'98)
C  ---------------------------------------------------------------
C     
C    Model  Label    u_LL  u_LR  u_RL  u_RR    d_LL  d_LR  d_RL  d_RR
C
C    101     VV       +     +     +     +       +     +     +     +
C
C    102     AA       +     -     -     +       +     -     -     +
C
C    103     VA       +     -     +     -       +     -     +     -
C
C    104     X1       +     -                   +     -
C
C    105     X2       +           +             +           +
C
C    106     X3       +                 +       +                 +
C
C    107     X4             +     +                   +     +
C
C    108     X5             +           +             +           +
C
C    109     X6                   +     -                   +     -
C
C    110     U1       +     -              
C
C    111     U2       +           +        
C
C    112     U3       +                 +  
C
C    113     U4             +     +        
C
C    114     U5             +           +  
C
C    115     U6                   +     -  
C
C  Extra Dimentions
C  ================
C
C   201 == 16     ED            (special coupling structure)
C
C  LQ scenarios
C  ============
C
C Scalar
C ------
C
C    301     S_o     1/2               1/2
C
C    302     S_o_L   1/2 
C
C    303     S_o_R                     1/2
C
C    304    ~S_o                                                 1/2
C
C    305     S_1/2        -1/2  -1/2                       -1/2
C
C    306     S_1/2_L      -1/2
C
C    307     S_1/2_R            -1/2                       -1/2
C
C    308    ~S_1/2                                   -1/2
C
C    309     S_1     1/2                        1
C
C
C Vector
C ------
C
C    311     V_o                               -1                -1
C
C    312     V_o_L                             -1
C
C    313     V_o_R                                               -1
C
C    314    ~V_o                       -1
C
C    315     V_1/2                1                   1     1
C
C    316     V_1/2_L                                  1
C
C    317     V_1/2_R              1                         1
C
C    318    ~V_1/2          1
C
C    319     V_1     -2                        -1
C 
C
C
C  Global (multiparameter fits)
C  ----------------------------
C
C  LW 25.01.16 -> not implementing yet.
C
C     24     -   4 parameter (assume u_?? = d_?? )
C
C     26     -   6 parameter (assume u_?L = d_?L )
C
C     27     -   7 parameter  (assume u_RL = d_RL)
C
C     28     -   8 parameter fit
C
C
C Test fits (one Etha at the time)
C
C  LW 25.01.16 -> not implementing yet.
C
C    -n      -   only parameter n fitted, n=1,8
C
C   -1n      -   only parameter n fitted, but assuming u_??=d_??, n=1,4
C
C
C======================================================================
C
C
C
CCC      Integer Ipar

CC        print *,'ModImpose:  CIindex = ', CIindex, 
CC     $ '  CIvarval = ', CIvarval 

C
C  VV
C
      If(CIindex.EQ.101)then
        par(1) =  CIvarval
        par(2) =  CIvarval 
        par(3) =  CIvarval
        par(4) =  CIvarval
        par(5) =  CIvarval
        par(6) =  CIvarval
        par(7) =  CIvarval
        par(8) =  CIvarval
C
C  AA
C
      ElseIf(CIindex.EQ.102)then
        par(1) =  CIvarval
        par(2) = -CIvarval
        par(3) = -CIvarval
        par(4) =  CIvarval
        par(5) =  CIvarval
        par(6) = -CIvarval
        par(7) = -CIvarval
        par(8) =  CIvarval
C
C  VA
C
      ElseIf(CIindex.EQ.103)then
        par(1) =  CIvarval
        par(2) = -CIvarval
        par(3) =  CIvarval
        par(4) = -CIvarval
        par(5) =  CIvarval
        par(6) = -CIvarval
        par(7) =  CIvarval
        par(8) = -CIvarval
C
C  X1
C
      ElseIf(CIindex.EQ.104)then
        par(1) =  CIvarval
        par(2) = -CIvarval
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval
        par(6) = -CIvarval
        par(7) =  0
        par(8) =  0
C
C  X2
C
      ElseIf(CIindex.EQ.105)then
        par(1) =  CIvarval
        par(2) =  0
        par(3) =  CIvarval
        par(4) =  0
        par(5) =  CIvarval
        par(6) =  0
        par(7) =  CIvarval
        par(8) =  0
C
C  X3
C
      ElseIf(CIindex.EQ.106)then
        par(1) =  CIvarval
        par(2) =  0
        par(3) =  0
        par(4) =  CIvarval
        par(5) =  CIvarval
        par(6) =  0
        par(7) =  0
        par(8) =  CIvarval
C
C  X4
C
      ElseIf(CIindex.EQ.107)then
        par(1) =  0
        par(2) =  CIvarval
        par(3) =  CIvarval
        par(4) =  0
        par(5) =  0
        par(6) =  CIvarval
        par(7) =  CIvarval
        par(8) =  0
C
C  X5
C
      ElseIf(CIindex.EQ.108)then
        par(1) =  0
        par(2) =  CIvarval
        par(3) =  0
        par(4) =  CIvarval
        par(5) =  0
        par(6) =  CIvarval
        par(7) =  0
        par(8) =  CIvarval
C
C  X6
C
      ElseIf(CIindex.EQ.109)then
        par(1) =  0
        par(2) =  0
        par(3) =  CIvarval
        par(4) = -CIvarval
        par(5) =  0
        par(6) =  0
        par(7) =  CIvarval
        par(8) = -CIvarval
C
C
C  U1
C
      ElseIf(CIindex.EQ.110)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval
        par(6) = -CIvarval
        par(7) =  0
        par(8) =  0
C
C  U2
C
      ElseIf(CIindex.EQ.111)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval
        par(6) =  0
        par(7) =  CIvarval
        par(8) =  0
C
C  U3
C
      ElseIf(CIindex.EQ.112)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval
        par(6) =  0
        par(7) =  0
        par(8) =  CIvarval
C
C  U4
C
      ElseIf(CIindex.EQ.113)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  CIvarval
        par(7) =  CIvarval
        par(8) =  0
C
C  U5
C
      ElseIf(CIindex.EQ.114)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  CIvarval
        par(7) =  0
        par(8) =  CIvarval
C
C  U6
C
      ElseIf(CIindex.EQ.115)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  CIvarval
        par(8) = -CIvarval
C  LL    
      ElseIf(CIindex.EQ.116)then
        par(1) =  CIvarval
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval
        par(6) =  0
        par(7) =  0
        par(8) =  0
C  LR    
      ElseIf(CIindex.EQ.117)then
        par(1) =  0
        par(2) =  CIvarval
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  CIvarval
        par(7) =  0
        par(8) =  0
C  RL    
      ElseIf(CIindex.EQ.118)then
        par(1) =  0
        par(2) =  0
        par(3) =  CIvarval
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  CIvarval
        par(8) =  0

C  RR 
      ElseIf(CIindex.EQ.118)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  CIvarval
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  CIvarval

C
C
C  6 parameter fit
C
CCC      ElseIf(Model.EQ.26)then
CCC        Par(5)=  par(1)
CCC        Par(7)=  par(3)
C
C  7 parameter fit
C
CCC      ElseIf(Model.EQ.27)then
CCC        Par(7)=  par(3)
CCC
C
C  Extra Dimensions
C
      ElseIf(CIindex.EQ.201)then
        par(1) =  CIvarval
        par(2) =  CIvarval
        par(3) =  CIvarval
        par(4) =  CIvarval
        par(5) =  CIvarval
        par(6) =  CIvarval
        par(7) =  CIvarval
        par(8) =  CIvarval
C
C  S_o
C
      ElseIf(CIindex.EQ.301)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval * 0.5
        par(6) =  0
        par(7) =  0
        par(8) =  CIvarval * 0.5
C
C  S_o_L
C
      ElseIf(CIindex.EQ.302)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval * 0.5
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  S_o_R
C
      ElseIf(CIindex.EQ.303)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  CIvarval * 0.5
C
C  ~S_o
C
      ElseIf(CIindex.EQ.304)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  CIvarval * 0.5
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  S_1/2
C
      ElseIf(CIindex.EQ.305)then
        par(1) =  0
        par(2) =  0
        par(3) = -CIvarval * 0.5
        par(4) =  0
        par(5) =  0
        par(6) = -CIvarval * 0.5
        par(7) = -CIvarval * 0.5
        par(8) =  0
C
C  S_1/2_L
C
      ElseIf(CIindex.EQ.306)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) = -CIvarval * 0.5
        par(7) =  0
        par(8) =  0
C
C  S_1/2_R
C
      ElseIf(CIindex.EQ.307)then
        par(1) =  0
        par(2) =  0
        par(3) = -CIvarval * 0.5
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) = -CIvarval * 0.5
        par(8) =  0
C
C  ~S_1/2
C
      ElseIf(CIindex.EQ.308)then
        par(1) =  0
        par(2) = -CIvarval * 0.5
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  S_1
C
      ElseIf(CIindex.EQ.309)then
        par(1) =  CIvarval
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  CIvarval * 0.5
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  V_o
C
      ElseIf(CIindex.EQ.311)then
        par(1) = -CIvarval
        par(2) =  0
        par(3) =  0
        par(4) = -CIvarval
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  V_o_L
C
      ElseIf(CIindex.EQ.312)then
        par(1) = -CIvarval
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  V_o_R
C
      ElseIf(CIindex.EQ.313)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) = -CIvarval
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  ~V_o
C
      ElseIf(CIindex.EQ.314)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) = -CIvarval

C
C  V_1/2
C
      ElseIf(CIindex.EQ.315)then
        par(1) =  0
        par(2) =  CIvarval
        par(3) =  CIvarval
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  CIvarval
        par(8) =  0
C
C  V_1/2_L
C
      ElseIf(CIindex.EQ.316)then
        par(1) =  0
        par(2) =  CIvarval
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  0
        par(8) =  0
C
C  V_1/2_R
C
      ElseIf(CIindex.EQ.317)then
        par(1) =  0
        par(2) =  0
        par(3) =  CIvarval
        par(4) =  0
        par(5) =  0
        par(6) =  0
        par(7) =  CIvarval
        par(8) =  0
C
C  ~V_1/2
C
      ElseIf(CIindex.EQ.318)then
        par(1) =  0
        par(2) =  0
        par(3) =  0
        par(4) =  0
        par(5) =  0
        par(6) =  CIvarval
        par(7) =  0
        par(8) =  0
C
C  V_1
C
      ElseIf(CIindex.EQ.319)then
        par(1) = -CIvarval
        par(2) =  0 
        par(3) =  0
        par(4) =  0
        par(5) = -CIvarval * 2.0
        par(6) =  0
        par(7) =  0
        par(8) =  0

      EndIf
C
C In most models u sector is the same as d sector
C
CCC      If(Model.LT.-10 .or. 
CCC     +  (Model.GE.1 .AND. Model.LE.9) .OR.
CCC     +   Model.EQ.24                           )then
CCC
CCC         Do Ipar=1,4
CCC            Par(Ipar+4)=par(Ipar)
CCC            EndDo
CCC
CCC         Endif
CCC------------HERE WE FORM ETA FROM PAR

      Eta = reshape([
     $ par(1), par(2), par(3), par(4), par(5), par(6), par(7), par(8),     
     $ par(1), par(2), par(3), par(4), par(5), par(6), par(7), par(8),
     $ par(1), par(2), par(3), par(4), par(5), par(6), par(7), par(8)  
     $], [4,6])
      write (190,*) 'CIvarval (in ModImpose)=',CIvarval
        
      Return
      End
C
C=================================================================
C
       Subroutine KQCD(x,Ucor,Dcor)
C
C No K_QCD for electron-photon
C
       DOUBLE PRECISION Ucor,Dcor,x
       Ucor=1.
       Dcor=1.

       Return
       End
C
C====================================================================
C
C As ContDen routine is linked instead of David Williams "contact2"
C package, we still have to include what is missing:
C
C     ==========================================
C     ContAlph
C     ==========================================
C
C     Purpose:
C       Calculate running alpha-em
C
C     Arguments:
C       Q2  - Scale (GeV^2)
C
C     Return value:
C       Alpha at Q2
C
C     Notes:
C       * Code stolen directly out of Jetset 7.3
C*********************************************************************
C*                   Copyright Torbjorn Sjostrand                   **
C*                                                                  **
C*                    CERN/TH, CH-1211 Geneva 23                    **
C*                BITNET/EARN address TORSJO@CERNVM                 **
C*********************************************************************
C
      DOUBLE PRECISION Function ContAlph( Q2 )
      IMPLICIT NONE
      DOUBLE PRECISION Q2
C
      DOUBLE PRECISION AEMPI, RPIGG
C
      DOUBLE PRECISION PARU_101, PARU_001
      Data PARU_001/3.14159265/,
     +     PARU_101/0.00729735/

C...Purpose: to calculate the running alpha_electromagnetic.

C...Calculate real part of photon vacuum polarization.
C...For leptons simplify by using asymptotic (Q^2 >> m^2) expressions.
C...For hadrons use parametrization of H. Burkhardt et al.
C...See R. Kleiss et al, CERN 89-08, vol. 3, pp. 129-131.
      AEMPI=PARU_101/(3.*PARU_001)
      IF(Q2.LT.2E-6) THEN
        RPIGG=0.
      ELSEIF(Q2.LT.0.09) THEN
        RPIGG=AEMPI*(13.4916+LOG(Q2))+0.00835*LOG(1.+Q2)
      ELSEIF(Q2.LT.9.) THEN
        RPIGG=AEMPI*(16.3200+2.*LOG(Q2))+0.00238*LOG(1.+3.927*Q2)
      ELSEIF(Q2.LT.1E4) THEN
        RPIGG=AEMPI*(13.4955+3.*LOG(Q2))+0.00165+0.00299*LOG(1.+Q2)
      ELSE
        RPIGG=AEMPI*(13.4955+3.*LOG(Q2))+0.00221+0.00293*LOG(1.+Q2)
      ENDIF
C...Calculate running alpha_em.
      ContAlph = PARU_101/(1.-RPIGG)
      RETURN
      END
C
C====================================================================
C

C
C===================================================================
C
C     ==========================================
C     ContTamp
C     ==========================================
C     Purpose:
C       Calculates qe->qe amplitudes
C       with and without a contact term.
C
C     Arguments:
C       S        - (in) Center of mass energy squared (GeV^2)
C       T        - (in) Four-momentum transfer (-Q^2)
C       Eta(4,6) - (in) Size and sign of contact term coefficients:
C       Zmass0   - (in) Mass of Z boson. Set to -1 for default.
C       Alpha    - (in) Alpha for photon. 0=fixed (true born
C       cross-section),
C                       -1=running, >0 that value used. -1 recommended.
C
C       Ampl(4,6)- (out) calculated amplitude squares (M(s,t)^2)
C       Status   - (out) Error status.
C
      SUBROUTINE ContTamp( S, T, Eta, Zmass0, Alpha, Ampl, Status )
C
      IMPLICIT NONE
      DOUBLE PRECISION S, T, Zmass0, Alpha
      DOUBLE PRECISION  Ampl(4,6)
C      Real*8 Eta(4,6)
      Integer Status
C
C Model parameters - LQ mass
C
CCCC      include 'cimodel.inc'
      include 'steering.inc'
      include 'couplings.inc'
      include 'CI_models.inc'
C
C     --- First: a few constants
C
C      Real*8 Pi, Qcharge(6), Qiso(6), Zmass, ZmassDef, Sin2T,
C     +     AlphaB, Zgamma
      DOUBLE PRECISION Qcharge(6), Qiso(6), Zmass

C
CCC      Data Pi/3.14159265/
C
C     Quark charges, times 3:
      Data Qcharge/-1, 2,-1, 2,-1,2/
C
C     Quark strong isospin, times 2:
      Data Qiso/-1, 1, -1, 1, -1, 1/ 
C
C     Default mass of the Z:
CCC      Data ZmassDef/91.1882/
C      Data Zgamma/2.490/
C
C     Weinberg angle (z-scheme, including oblique radiative corrections,
C     which is strickly not correct, since we're not at Q^2 = Mz^2).
CCC      Data Sin2T/0.232/
C
C     Alpha (bare)
CCC      Data AlphaB/7.29735E-3/
C
C     --- and derived constants
C
      Logical Initialized
      Data Initialized/.FALSE./
      DOUBLE PRECISION Cos2T, gLq(6), gRq(6), gLe, gRe
      Save Cos2T, gLq, gRq, gLe, gRe
C
C     --- Program variables
C
      Integer Ic, Iq
      DOUBLE PRECISION    e2cor, g1, g2, gg, g3,g3R(6),g3I(6),Wcorr(6)
      DOUBLE PRECISION    ReDen,ImDen, Den, ReAmp,ImAmp
C
C  Function calculating running alpha_EM
C
      DOUBLE PRECISION    ContAlph,Ucor,Dcor
      DOUBLE PRECISION U, LQEX,LQefwid(6)
      DOUBLE PRECISION GrCorr(4)
      
C LW: additional parameters
      DOUBLE PRECISION LQFN, LQSign, LQSP

      double precision hf_get_alphas
      double precision hf_q2
      
C      Integer Mvar

C LW: set additional parameters
      If(CIindex.GE.301)Then

        If(CIindex.LE.304 .OR. CIindex.EQ.309 .OR.
     +   (CIindex.GE.315 .AND. CIindex.LE.318) )then
           LQFN=2
           LQSign=+1
        Else
           LQFN=0
           LQSign=-1
        EndIf
        
        IF(CIindex.GE.311)then
           LQSP=1
        Else
           LQSP=0
        EndIf
      
      Else
        LQsign=0
      EndIf  

C      if(CIindex.LE.106)then
C         Mvar=1
C      elseif(CIindex.LE.108)then
C         Mvar=2
C      elseif(CIindex.EQ.109)then
C         Mvar=3
C      elseif(CIindex.LE.112)then
C         Mvar=5
C      elseif(CIindex.LE.114)then
C         Mvar=6
C      elseif(CIindex.EQ.115)then
C         Mvar=7
C      elseif(CIindex.EQ.201)then
C         Mvar=9
C
C LQ scenarios
C
C      elseif(CIindex.EQ.309 .OR. CIindex.EQ.311 .OR. 
C      If(CIindex.EQ.309 .OR. CIindex.EQ.311 .OR. 
C     +       CIindex.EQ.312 .OR. CIindex.EQ.319 )then
C         Mvar=1
C      elseif(CIindex.EQ.308 .OR. CIindex.EQ.315 .OR.
C     +       CIindex.EQ.316 )then
C         Mvar=2
C      elseif(CIindex.EQ.305 .OR. CIindex.EQ.307 .OR. 
C     +       CIindex.EQ.317 )then
C         Mvar=3
C      elseif(CIindex.EQ.304 .OR. CIindex.EQ.313 )then
C         Mvar=4
C      elseif(CIindex.EQ.301 .OR. CIindex.EQ.302 )then
C         Mvar=5
C      elseif(CIindex.EQ.306 .OR. CIindex.EQ.318 )then
C         Mvar=6
C      elseif(CIindex.EQ.303 .OR. CIindex.EQ.314 )then
C         Mvar=8
C      EndIf
      
C LW: undefined (temporaly) LQ parameters
C      LQmass = 0
C      LQwidth = 0
C      LQres = 0
C
      IF (.NOT.Initialized) THEN
C
C       --- Fill in derived constants
C
        Cos2T = 1 - sin2thw
CC        print *,'     ContTamp check, sin2thw = ', sin2thw,
CC     +          '  Cos2T = ', Cos2T
C
C       --- Electroweak couplings.
C
C  electron
C
        gLe = -0.5 + sin2thw
        gRe = sin2thw
C
C  quarks
C
        DO Iq = 1, 6
          Qcharge(Iq) = Qcharge(Iq)/3.0
          Qiso(Iq) = Qiso(Iq)/2.0

          gLq(Iq) = Qiso(Iq) - Qcharge(Iq)*sin2thw
          gRq(Iq) = - Qcharge(Iq)*sin2thw

        ENDDO
C
        Initialized = .TRUE.
      ENDIF
C
C----------------------------------------------------------
C
C     --- Check input values
C
      IF (T.GT.-1E-20) THEN
        Status = -1
        RETURN
        ENDIF
C
C Neglect quark and electron masses
C
      U = -S - T
C
C   Zo mass
C
      IF (Zmass0.GE.0) THEN
        Zmass = Zmass0
      ELSE
CCC        Zmass = ZmassDef
        Zmass = Mz
CC        print *,'     ContTamp check, Mz = ', Mz,
CC     +          '  Zmass = ', Zmass
      ENDIF
C
C Alpha EM
C
CC      print *,'     ContTamp check, Alpha = ', Alpha
      IF (Alpha.LT.0.0) THEN
C LW: running alphaem, copied from Filip's "mypdf.f"
        e2cor = 4.*Pi*ContAlph(ABS(T))
      ELSE IF (Alpha.EQ.0.0) THEN
        e2cor = 4.*pi*alphaem
C        print *,'     ContTamp check, alphaem = ', alphaem,
C     +   '  e2cor = ', e2cor
      ELSE
        e2cor = 4.*pi*Alpha
C        print *,'     ContTamp check, Alpha = ', Alpha,
C     +   '  e2cor = ', e2cor
      ENDIF
C
C LQ coupling
C ===========
C
      Do iq=1,6
        LQefwid(iq)=0.0
        Wcorr(iq)=1.0
      EndDo

      LQEX=0.0

      If(CIindex.GE.301 .AND. CILQmass.GT.0.0  .AND.
CCC     +   (LQwidth.GT.0. .OR. Eta(Mvar,1).NE.0.) )Then
     +   (CILQwidth.GT.0. .OR. CIvarval.NE.0.) )Then
C
C Calculate LQwidth if not available
C
         If(CILQwidth.LE.0.)then
C
C general formula for vector LQ
C
CCC            LQefwid(1)=LQsign*Eta(Mvar,1)*LQmass**3/(24.0*Pi)
            LQefwid(1)=LQsign*CIvarval*CILQmass**3/(24.0*pi)
C
C for scalar LQ: factor 3/2 from spin * factor 2 from eta-coupling
C relation
C
            If(LQSP.EQ.0)LQefwid(1)=LQefwid(1)*3.0
C
C Correct for BR (CC decay channel: models 39 and 49)
C         
            IF(Eta(1,1).NE.Eta(1,2))LQefwid(1)=LQefwid(1)*2.0
C
C Store -width for debug purposes :-)
C
            CILQwidth=-LQefwid(1)
C
            Do iq=2,6
              Lqefwid(iq)=Lqefwid(1)
              EndDo
C
C Take into account QCD correction
C
C LW TMP, missing func:
             Call KQCD(DBLE(-1.),Ucor,Dcor)
C
C             Dcor = 1
C             Ucor = 1
C
             Lqefwid(1)=Lqefwid(1)/Dcor
             Lqefwid(2)=Lqefwid(2)/Ucor
C
         Else
            Do iq=1,6
              LQefwid(iq)=CILQwidth
            EndDo
         EndIf
C
C Avoid to small widths
C
         Do iq=1,6
            If(LQefwid(Iq).LT.CILQres)then
              Wcorr(iq)=CILQres/LQefwid(Iq)
              LQefwid(Iq)=CILQres
            EndIf
         EndDo
 
C
C LQ production/exchange channel
C
         If(LQFN.EQ.2)Then
           LQEX=S
         ELSE
           LQEX=U
         EndIf
        
      EndIf
C
C  pure foton coupling
C
      g1= e2cor/T
C
C  terms coming from Z propagator - same alpha_em as for photon (?)
C
      g2= e2cor/(sin2thw*cos2T)/(T - Zmass*Zmass)
C
C      g2= 4.*Pi*AlphaB/(sin2T*cos2T)/(T - Zmass*Zmass)
C
C CI/LQ term
C ==========
C CI approximation:
C
      If(LQefwid(1).LE.0.0)Then

        DO Iq = 1, 6
          g3R(iq) = 1.0
          g3I(Iq) = 0.0
          EndDo
C
C u-channel exchange
C
      ElseIf( LQEX.LE.0.0)Then

        DO Iq = 1, 6
          g3R(Iq) = (CILQmass*CILQmass)/(CILQmass*CILQmass - LQEX)
          g3I(Iq) = 0.
          EndDo

C
C s-channel production
C
      Else

        DO Iq = 1, 6
          g3 = (CILQmass*CILQmass - LQEX)**2 +
     &    (LQEX*LQefwid(iq)/CILQmass)**2
          g3R(iq) = (CILQmass*CILQmass - LQEX)*(CILQmass*CILQmass)/g3
C
C Corection to interference term
C
C          g3R(Iq) = g3R(Iq)/sqrt(Wcorr(Iq))
C          g3R(iq) = g3R(Iq)/Wcorr(Iq)    
C
C ===> No correction! (as found out in tests)
C
          If(Wcorr(Iq).le.1.0)then
            g3I(iq) = (LQEX*LQefwid(iq)*CILQmass)/g3
          Else
C
C No correction to interference term:
C
            g3i(iq) = 
     +        (Wcorr(iq)-1.0)*(CILQmass*CILQmass - LQEX)**2 +
     +         Wcorr(iq)*(LQEX*LQefwid(iq)/CILQmass)**2
C
C Gamma^(-0.5) correction to interference term
C
C            g3i(iq) = 
C     +        (Wcorr(iq)-1/Wcorr(iq))*(LQmass*LQmass - LQEX)**2 +
C     +         Wcorr(iq)*(LQEX*LQefwid(iq)/LQmass)**2
C
C Ganna^(-1) correction to interference term
C
C            g3i(iq) = 
C     +        (Wcorr(iq)-1/Wcorr(iq)**2)*(LQmass*LQmass - LQEX)**2 +
C     +         Wcorr(iq)*(LQEX*LQefwid(iq)/LQmass)**2
C
            g3i(iq)=sqrt(g3i(iq))*CILQmass*CILQmass/g3

            EndIf
          EndDo

        EndIf
C
C Kinemaitc corrections for model 16 (Extra dimentions)
C

       If(CIindex.EQ.201)Then
          GrCorr(1)= -3.1415926/2. * (4.*u+t)/1000000.
          GrCorr(2)= -3.1415926/2. * (4.*u+3.*t)/1000000.
          GrCorr(3)=GrCorr(2)
          GrCorr(4)=GrCorr(1)

       Else
          Do ic=1,4
            GrCorr(ic)=1.
          EndDo
       EndIf
C
C Final calculation
C
       Do iq=1,6
         Do Ic=1,4
            
           If(Ic.EQ.1)then
              gg = gLe * gLq(iq)
           elseIf(Ic.EQ.2)then
              gg = gLe * gRq(iq)
           elseIf(Ic.EQ.3)then
              gg = gRe * gLq(iq)
           else
              gg = gRe * gRq(iq)
           EndIf
C
           ReAmp = -g1 * Qcharge(Iq) + g2 * gg  + 
     +                                 g3R(iq)*GrCorr(ic)*Eta(ic,iq)
           ImAmp = g3I(iq)*Eta(ic,iq)
C
           Ampl(ic,iq) = ReAmp*ReAmp + ImAmp*ImAmp
C
         EndDo

      EndDo
C
      Status = 0
C
      RETURN
      END
C
C===================================================================
C

C
C===================================================================
C
C     ==========================================
C     ContNC
C     ==========================================
C     Purpose:
C       Calculates for a particular x, Q2, and S the tree-level NC DIS
C       differential cross-section with and without an addition contact
C       term. Can be used to reweigh ordinary NC DIS Monte Carlo and
C       investigate the effects of addition contact terms.
C
C     Arguments:
C       X        - (in) Kinematic variable X (0.0 < X <= 1.0)
C       Q2       - (in) Kinematic variable Q^2 (GeV^2)
C       S        - (in) Center of mass energy squared (GeV^2)
C       Eta(4,6) - (in) Size and sign of contact term coefficients:
C                    Eta(1,i) = left-left coupling 
C                    Eta(2,i) = left-right   "
C                    Eta(3,i) = right-left   "
C                    Eta(4,i) = rigth-right  "
C                    i = 1  for down quark
C                        2  for up 
C                        3  for strange
C                        4  for charm
C                        5  for bottom
C                        6  for top
C       Electron - (in) TRUE if electron, otherwise positron
C       Polar    - (in) Electron/Positron polarization
C       XQfract  - (in) X times Fraction of quarks in hadron:
C                    XQfract(i,1) = X * fraction down 
C                    XQfract(i,2) = X * fraction up
C                    XQfract(i,3) = X * fraction strange
C                    XQfract(i,4) = X * fraction charm
C                    XQfract(i,5) = X * fraction bottom
C                    XQfract(i,6) = X * fraction top (if you really want
C                    it!)
C                  i=1 is quark, i=2 is anti-quark.
C                  Reasonable values should have SUM(Qfract) < 1/X,
C                  although this is not explicitly checked.
C       Zmass    - (in) Mass of Z boson. Set to -1 for default.
C       Alpha    - (in) Alpha for photon. 0=fixed (true born
C       cross-section),
C                       -1=running, >0 that value used. -1 recommended.
C                       See routine ContAlph.
C       DisCross - (out) NC Dis tree-level differential cross-sect
C       (pb/GeV^2).
C       ConCross - (out) NC + Contact term DIS tree-level differential 
C                  cross-section (pb/GeV^2).
C       Status   - (out) Error status. If zero, there is no problem.
C                  IF < 0, DisCross and ConCross are undefined.
C                  Errors:
C                      1 = X is unphysical ( < Q2/S )
C                     -1 = Divergent input argument(s)
C                     -2 = ABS(Polar) > 1.0
C
C     Convention: 
C        The contact term is calculated by introducing
C        addition terms in the standard model lagrangian:
C
C        L = Lsm + Eta1(eL-g-eL)(qL-g-qL) + Eta2(eL-g-eL)(qR-g-qR)
C                + Eta3(eR-g-eR)(qL-g-qL) + Eta4(eR-g-eR)(qR-g-qR)
C
C
      SUBROUTINE DContNC( X, Q2, S, Eta, Electron, Polar, Zmass,
     +                    Alpha, XQfract, DisCross, ConCross, Status )
      IMPLICIT NONE
      DOUBLE PRECISION X, Q2, S, XQfract(2,7), Polar, Zmass, Alpha
      DOUBLE PRECISION DisCross, ConCross
C      Real*8 Eta(4,6)
      Integer Status
      Logical Electron
C
C Model parameters 
C
CCC      include 'cimodel.inc'
      include 'steering.inc'
      include 'couplings.inc'
      include 'CI_models.inc'
C
CCC      Real*8 Pi, GeVNb
CCC      Data Pi/3.14159265/
C
C     Conversion from GeV^-2 to nbarn
      DOUBLE PRECISION GeVNb
      Data GeVNb/0.38937966E9/
C
      DOUBLE PRECISION Eta0(4,6)
      Data Eta0/24*0./
      DOUBLE PRECISION Ampl0(4,6),Ampl(4,6),AmplB(4,6),EDgg
      DOUBLE PRECISION y, T1, T2, yfac(4), yfacB(4), PL, PR
      DOUBLE PRECISION T,U,SS
CCCC      Real Ucor,Dcor
      Integer Iq,Ilr
C
C     --- Check input values
C
      IF ((X.LT.0.OR.X.GT.1.0).OR.(S*X.LT.1E-10) ) THEN
        Status = -1
        print *,'   DContNC error, X: ', X, '   S: ', S
        RETURN
      ELSE IF (ABS(Polar).GT.1.0) THEN
        Status = -2
        RETURN
      ENDIF
C
C     --- Invoke ContTamp
C
      SS = X*S
      T = -Q2
      U =  Q2 - SS 
C
C LW: don`t see the reason for this call
CC      Call KQCD(X,Ucor,Dcor)
C
      Call ContTamp( SS, T, Eta0, Zmass, Alpha, Ampl0, Status )

      If (Electron) Then
        Call ContTamp( SS, T, Eta , Zmass, Alpha, Ampl , Status )
        Call ContTamp(  U, T, Eta , Zmass, Alpha, AmplB, Status )
      else
        Call ContTamp(  U, T, Eta , Zmass, Alpha, Ampl , Status )
        Call ContTamp( SS, T, Eta , Zmass, Alpha, AmplB, Status )
      EndIf
C
      IF (Status.NE.0) RETURN
C
C     --- Build cross-sections
C
      y = Q2/(S*X)

      IF(Electron)then
          T1 = 1
          T2 = (1-y)**2 
          PL = (1.0 - Polar)
          PR = (1.0 + Polar)
      else
          T1 = (1-y)**2 
          T2 = 1
          PL = (1.0 + Polar)
          PR = (1.0 - Polar)
      endif

      yfac(1)=T1*PL
      yfac(2)=T2*PL
      yfac(3)=T2*PR
      yfac(4)=T1*PR

      yfacB(1)=T2*PL
      yfacB(2)=T1*PL
      yfacB(3)=T1*PR
      yfacB(4)=T2*PR

      DisCross = 0.
      ConCross = 0.

      Do iq=1,6

        Do ilr=1,4

          DisCross = DisCross + XQfract(1,iq)*Ampl0(ilr,iq)*yfac(ilr)
     +                        + XQfract(2,iq)*Ampl0(ilr,iq)*yfacB(ilr)

          ConCross = ConCross + XQfract(1,iq)*Ampl(ilr,iq) *yfac(ilr)
     +                        + XQfract(2,iq)*AmplB(ilr,iq)*yfacB(ilr)

          EndDo

        EndDo
C
C Special addition for Large Extra Dimentions (Model=16)
C                 - scattering on gluons !
C It is assumed that all Eta(i) are equal for this model ( = 1/M_s^4)
C
      If(CIindex.eq.201 .and. XQfract(1,7).GT.0)then
         EDgg= -U*(U*U+SS*SS)/SS
         EDgg= 2*pi*pi* Eta(1,1)*Eta(1,1)/1.0D+12 * EDgg * XQfract(1,7)
         ConCross = ConCross + EDgg
         EndIf
C
      DisCross = GeVnb/(16.*pi*X)  * DisCross
      ConCross = GeVnb/(16.*pi*X)  * ConCross 
C
      IF (S*X.LT.Q2) THEN
        Status = 1
      ELSE
        Status = 0
      ENDIF
C
      RETURN
      END
C
C=========================================================================
C
C     ==========================================
C     ContCC
C     ==========================================
C     Purpose:
C       Calculates for a particular x, Q2, and S the tree-level CC DIS
C       differential cross-section with and without an addition contact
C       term. To be used to weigh an ordinary CC DIS Monte Carlo and
C       thus
C       convert it into a contact term Monte Carlo.
C
C     ------------------------------------------------------------
C     Rewritten by A.F.Zarnecki - consider only effects induced by 
C     NC CI couplings, when assuming SU(2) x U(1).
C
C     !!! Now including also proper LQ mass term !!!
C
C     ------------------------------------------------------------
C
C     Arguments:
C       X        - (in) Kinematic variable X (Q2/S <= X <= 1.0)
C       Q2       - (in) Kinematic variable Q^2 (GeV^2, > 0)
C       S        - (in) Center of mass energy squared (GeV^2)
C       Eta(3)   - (in) Size and sign of contact term coefficients
C                       induced by imposing SU(2)xU(1)
C                       Eta(1) = Eta^ed_LL - Eta^eu_LL
C                       Eta(2) = Eta^es_LL - Eta^ec_LL
C                       Eta(3) = Eta^eb_LL - Eta^et_LL
C
C       Electron - (in) TRUE if electron, otherwise positron
C       Polar    - (in) Electron/Positron polarization
C       XQfract  - (in) X times Fraction of quarks in hadron:
C                    XQfract(i,1) = X * fraction down 
C                    XQfract(i,2) = X * fraction up
C                    XQfract(i,3) = X * fraction strange
C                    XQfract(i,4) = X * fraction charm
C                    XQfract(i,5) = X * fraction bottom
C                    XQfract(i,6) = X * fraction top (if you really want
C                    it!)
C                  i=1 is quark, i=2 is anti-quark.
C                  Reasonable values should have SUM(Qfract) < 1/X,
C                  although this is not explicitly checked.
C       DisCross - (out) NC Dis tree-level differential cross-sect
C       (pb/GeV^2).
C       ConCross - (out) NC + Contact term DIS tree-level differential 
C                  cross-section (pb/GeV^2).
C       Status   - (out) Error status. If zero, there is no problem.
C                  IF < 0, DisCross and ConCross are undefined.
C                  Errors:
C                      1 = X is unphysical ( < Q2/S )
C                     -1 = Divergent input argument(s)
C
C     Convention: 
C        The contact term is calculated by introducing
C        additional terms in the standard model lagrangian:
C
C        eta * (e*gam_u*(1+gam_5)*n)(u*gam_u*(1+gam_5)*d)   
C
C        Here the contact terms do NOT couple KM mixed quark species.
C 
      SUBROUTINE DContCC( X, Q2, S, Eta3, Electron, 
     +             Polar, XQfract, DisCross, ConCross, Status )
      IMPLICIT NONE
      DOUBLE PRECISION X, Q2, S, XQfract(2,7), Polar
      DOUBLE PRECISION  Eta3(3), DisCross, ConCross
      Integer Status
      Logical Electron
C
C     --- First: a few constants
C

      include 'steering.inc'
      include 'couplings.inc'
      include 'CI_models.inc'

CC      Real*8 Pi, GeVNb, Wmass, Sin2T, KM(3,3)
      DOUBLE PRECISION GeVNb, KM(3,3)
C
C      Data Pi/3.14159265/
C
C     Conversion from GeV^-2 to nbarn
      Data GeVNb/0.38937966E9/
C
C     Mass of the W:
C      Data Wmass/80.41/
C
C     Weinberg angle (z-scheme, including oblique radiative
C     corrections):
C      Data Sin2T/0.232/
C
C     Kobayashi-Maskawa matrix (PDG 94, middle of range)
C      Data KM/ 0.975, 0.220, 0.003,
C     +         0.220, 0.975, 0.039,
C     +         0.009, 0.038, 0.999  /

C      Data KM/ Vud, Vus, Vub,
C     +         Vcd, Vcs, Vcb,
C     +         Vtd, Vts, Vtb  /

C
C     --- and derived constants
C
      Logical Initialized
      Data Initialized/.FALSE./
C
      DOUBLE PRECISION KM2(3,3)
      Save KM2
C
C     --- Program variables
C
      Integer I1, I2, Iup, Idn, Iu1, Iu2
      DOUBLE PRECISION  SumUK, SumDK, SumU1, SumD1,
     +        MQ2, Hratio, y, T1, Alpha, Tu
C
      DOUBLE PRECISION g3,CIfacSR,CIfacSI,CIfacX,CIfacU,CIfacD
      DOUBLE PRECISION U,Shat,LQefwid,Wcorr
      
C LW: additional parameters
      DOUBLE PRECISION LQFN
C
C Model parameters - LQ mass
C
CCC      include 'cimodel.inc'
      
C
C LW: set additional parameters
      If(CIindex.GE.301)Then
        If(CIindex.LE.304 .OR. CIindex.EQ.309 .OR.
     +   (CIindex.GE.315 .AND. CIindex.LE.318) )then
           LQFN=2
        Else
           LQFN=0
        EndIf
      EndIf

      KM = reshape((/ Vud, Vus, Vub,
     +                Vcd, Vcs, Vcb,
     +                Vtd, Vts, Vtb /), shape(KM))

CCCCC

C LW: undefined (temporaly) LQ parameters
C      LQmass = 0
C      LQwidth = 0
C      LQres = 0
      
CCCCC

      IF (.NOT.Initialized) THEN
C
C       --- Fill in derived constants
C
        DO I1 = 1, 3
          DO I2 = 1, 3
            KM2(I1,I2) = KM(I1,I2)*KM(I1,I2)
          ENDDO
        ENDDO
C
        Initialized = .TRUE.
      ENDIF
C
C     --- Check input values
C
      IF ((X.LT.0.OR.X.GT.1.0).OR.(S*X.LT.Q2).OR.Q2.LE.0.) THEN
        Status = -1
        RETURN

      ELSE IF (ABS(Polar).GT.1.0) THEN
        Status = -2
        RETURN

      ELSE
        Status = 0
      ENDIF
C
      DisCross=0.
      ConCross=0.
C
C     --- Calculate common factors
C
      Alpha = 1.0/128.0 
C
      MQ2 = sin2thw*(Mw*Mw + Q2)
      Hratio = MQ2/(2*pi*Alpha)
C
C LQ coupling
C ===========
C
C  !!! LQ width has to be set from outside !!!
C
      Shat=X*S
      U = -Shat + Q2

      If(CIindex.GE.301 .AND. CILQmass.GT.0.0 .AND. Eta3(1).NE.0.0 )Then

         If(CILQwidth.LE.0.0)then
            status=-2
            return            
         endif
C
C Avoid to small widths
C
         If(CILQwidth.LT.CILQres)then
            Wcorr=CILQres/CILQwidth
            LQefwid=CILQres
         else
            Wcorr=1.0
            LQefwid=CILQwidth
         EndIf

C
C u-channel exchange
C
         CIfacX = (CILQmass*CILQmass)/(CILQmass*CILQmass - U)
C
C s-channel production
C
        g3 = (CILQmass*CILQmass - Shat)**2 + (Shat*Lqefwid/CILQmass)**2
        CIfacSR = (CILQmass*CILQmass - Shat)
     &  *(CILQmass*CILQmass)/g3/sqrt(Wcorr)

        If(Wcorr.le.1.0)then
           CIfacSI = (Shat*Lqefwid*CILQmass)/g3
        Else
           CIfacSI = (Wcorr-1.0/Wcorr)*(CILQmass*CILQmass - Shat)**2 +
     +                           Wcorr*(Shat*Lqefwid/CILQmass)**2
           CIfacSI = (sqrt(CIfacSI)*CILQmass*CILQmass)/g3
C           CIfacSI = sqrt(CIfacSI)
C           CIfacSI = CIfacSI*CILQmass
C           CIfacSI = CIfacSI*CILQmass
C           CIfacSI = CIfacSI/g3
        EndIf

      EndIf
C
C     --- Build effective quark densities for electrons (positrons):
C
C            SumUK  - Sum of KM^2 x Up (anti-up) density
C            SumDK  - Sum of KM^2 x anti-Down (down) density
C            SumU1  - SumUK including CI contribution
C            SumD1  - SumDK including CI contribution
C
      IF (Electron) THEN
        Iu1 = 1
        Iu2 = 2
      ELSE
        Iu1 = 2
        Iu2 = 1
      ENDIF
C
      SumUK = 0.0
      SumDK = 0.0
      SumU1 = 0.0
      SumD1 = 0.0
C
C Loop over final state quark family
C
      DO I1 = 1, 3
C
C Loop over initial state quark family
C
        DO I2 = 1, 3
C
          Iup = I2*2
          Idn = Iup - 1
C
C CI/LQ contribution - only within one family !
C
          IF (I1.NE.I2 .OR. Eta3(1).EQ.0.0) THEN
            CIfacU = 1.0
            CIfacD = 1.0

          ELSEIF(CIindex.LT.301 .OR. CILQmass.LE.0.0)Then
            CIfacU = (1.0-Hratio*Eta3(I1)/KM(I1,I2))**2
            CIfacD = CIfacU

          ELSEIF(LQFN.EQ.0)Then
            CIfacU =  (1.0-Hratio*CifacX*Eta3(I1)/KM(I1,I2))**2
            CIfacD =  (1.0-Hratio*CifacSR*Eta3(I1)/KM(I1,I2))**2
     +                  + (Hratio*CifacSI*Eta3(I1)/KM(I1,I2))**2

          ELSE
            CIfacU =  (1.0-Hratio*CifacSR*Eta3(I1)/KM(I1,I2))**2
     +                  + (Hratio*CifacSI*Eta3(I1)/KM(I1,I2))**2
            CIfacD =  (1.0-Hratio*CifacX*Eta3(I1)/KM(I1,I2))**2
            ENDIF
C
C Sum of (weighted) up-type quark densities
C
          SumUK = SumUK + XQfract(Iu1,Iup)*KM2(I1,I2)
          SumU1 = SumU1 + XQfract(Iu1,Iup)*KM2(I1,I2)*CIfacU
C
C Sum of (weighted) down-type quark densities
C Exclude production of top !!!
C
          If(I1.NE.3)Then
             SumDK = SumDK + XQfract(Iu2,Idn)*KM2(I2,I1)
             SumD1 = SumD1 + XQfract(Iu2,Idn)*KM2(I2,I1)*CIfacD
             Endif
C
        ENDDO
      ENDDO
C
C     --- Build left (right) handed electron (positron) cross-sections
C
      y = Q2/(S*X)
      T1 = GeVNb*(pi*Alpha*Alpha)/(4.0*X*MQ2*MQ2)
      Tu = (1.0-y)
C
      DisCross = T1*(SumUK + Tu*Tu*SumDK)
C
      ConCross = T1*(SumU1 + Tu*Tu*SumD1)
C

C LW 30.01.16: adding lepton polarization

      IF(Electron)then
          DisCross = (1.0 - Polar) * DisCross
          ConCross = (1.0 - Polar) * ConCross
      else
          DisCross = (1.0 + Polar) * DisCross
          ConCross = (1.0 + Polar) * ConCross
      endif

      RETURN
      END
C

      Subroutine CIsavePdfLO
      IMPLICIT NONE
      
      Integer isdx
      include 'steering.inc'
      include 'CI_models.inc'
      open(112,file="CIpdfLO_out.txt",action="write",status="replace")
      
      Do isdx=1,5000
        write(112,*) CIXQfractLO(isdx,1,1), CIXQfractLO(isdx,2,1), 
     +               CIXQfractLO(isdx,1,2), CIXQfractLO(isdx,2,2), 
     +               CIXQfractLO(isdx,1,3), CIXQfractLO(isdx,2,3), 
     +               CIXQfractLO(isdx,1,4), CIXQfractLO(isdx,2,4), 
     +               CIXQfractLO(isdx,1,5), CIXQfractLO(isdx,2,5), 
     +               CIXQfractLO(isdx,1,6), CIXQfractLO(isdx,2,6), 
     +               CIXQfractLO(isdx,1,7), CIXQfractLO(isdx,2,7)
      EndDo
      
      close (112)
      
      Return
      End



      Subroutine CIreadPdfLO
      IMPLICIT NONE

      Integer irdx
      include 'steering.inc'
      include 'CI_models.inc'

      open(114,file="CIpdfLO_in.txt",action="read")

      Do irdx=1,5000
        read(114,*) CIXQfractLO(irdx,1,1), CIXQfractLO(irdx,2,1),
     +              CIXQfractLO(irdx,1,2), CIXQfractLO(irdx,2,2),
     +              CIXQfractLO(irdx,1,3), CIXQfractLO(irdx,2,3),
     +              CIXQfractLO(irdx,1,4), CIXQfractLO(irdx,2,4),
     +              CIXQfractLO(irdx,1,5), CIXQfractLO(irdx,2,5),
     +              CIXQfractLO(irdx,1,6), CIXQfractLO(irdx,2,6),
     +              CIXQfractLO(irdx,1,7), CIXQfractLO(irdx,2,7)
      EndDo

      close (114)

      Return
      End


      Subroutine CIprintLOratio
      IMPLICIT NONE

      DOUBLE PRECISION EtaParNC(4,6), EtaVal
      Integer iqg, iqt, iq, ilr
      Integer iNx, iNq2, iNelectron, iNeta
      Logical isElectron
      DOUBLE PRECISION XQfract(2,7), CIX(6), CIQ2, CIS, CIPolar, CIAlpha
      DOUBLE PRECISION Q2min0, Q2max0, Q2max, Q2logstep
      DOUBLE PRECISION DisCross, ConCross
      Integer Status

      double precision dbPdf, hfX, hfQ2
      dimension dbPdf(-6:6)

      character*64 OutFile, ChEta
      Integer iFile

      Data CIX/0.1, 0.2, 0.3, 0.4, 0.5, 0.6/
c      Data OutFile/"./CI_output/CI_LOration.txt"/
      
      include 'steering.inc'
      include 'CI_models.inc'

      CIS = (318.0)**2

      Q2min0 = 900
      Q2max0 = 70000

      CIPolar = 0
      CIAlpha = -1.

      EtaVal = (4 * 3.14159265) / (5000**2)

      Do iNeta=1,2

        if (iNeta.eq.1) then
          ChEta='eta+'
          Do iqg=1,6
            Do ilr=1,4
              EtaParNC(ilr,iqg) = EtaVal
            EndDo
          EndDo
        else
          ChEta='eta-'
          Do iqg=1,6
            Do ilr=1,4
              EtaParNC(ilr,iqg) = -EtaVal
            EndDo
          EndDo
        endif

        Do iNelectron=1,2

          if (iNelectron.eq.1) then
            isElectron=.true.
          else
            isElectron=.false.
          endif

          Do iNx=1,6

            if (isElectron) then
              write(OutFile, '(A14,F3.1,A4,A4)') 
     +          'LOratioNCe-pX=',CIX(iNx),ChEta,'.dat'
            else
              write(OutFile, '(A14,F3.1,A4,A4)') 
     +          'LOratioNCe+pX=',CIX(iNx),ChEta,'.dat'
            endif
          
            iFile = iNelectron*100 + iNx
            open(iFile,file=OutFile,action="write",status="replace")

            Q2max = CIS*CIX(iNx)*0.99
            if (Q2max.gt.Q2max0) then
              Q2max = Q2max0
            endif

            Q2logstep = (log(Q2max) - log(Q2min0)) / 499

            Do iNq2=1,500

              CIQ2 = exp(log(Q2min0) + Q2logstep*(iNq2-1))

              hfX = CIX(iNx)
              hfQ2 = CIQ2
              Call hf_get_pdfs(hfX,hfQ2,dbPdf)
              Do iq=1,6
                XQfract(1,iq) = dbPdf(iq)
                XQfract(2,iq) = dbPdf(-iq)
              EndDo

              Call DContNC( CIX(iNx), CIQ2, CIS, EtaParNC, isElectron,
     +            CIPolar,DBLE(-1.), CIAlpha, XQfract, DisCross, ConCross, 
     +            Status )

              write(iFile,*) CIQ2, (ConCross/DisCross), Status,CIX(iNx),
     +              XQfract

            EndDo

            close (iFile)

          EndDo
        EndDo
      EndDo

      Return
      End


