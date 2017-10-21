C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine fill_CI_rules(doCI_in, CIDoSimpFit_in, CISimpFirStep_in)
      include "steering.inc"

      logical, intent(out) :: doCI_in, CIDoSimpFit_in
      character*(*), intent(out) :: CISimpFirStep_in

      doCI_in = doCI
      CIDoSimpFit_in = CIDoSimpFit
      CISimpFirStep_in = CISimpFitStep
      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~      
      subroutine calc_theo(g_dummy, parminuit, iflag)
      implicit none

      include 'fcn.inc'
      include 'endmini.inc'
      include 'for_debug.inc'
      include 'steering.inc'
      include 'pdfparam.inc'
      include 'alphas.inc'
      include 'couplings.inc'
      include 'ntot.inc'
      include 'datasets.inc'
      include 'systematics.inc'
      include 'theo.inc'
      include 'indata.inc'
      include 'thresholds.inc'
      include 'polarity.inc'
      include 'fractal.inc'

      double precision g_dummy(*), parminuit(*)
      integer iflag

      integer i, j, kflag, k, nwds, idataset
      double precision alphaszero, hf_get_alphas, epsi
      external LHAPDFsubr

C     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --

*    Store params in a common block:
      do i=1,MNE
         parminuitsave(i) = parminuit(i)
      end do

*    PDF parameterisation at the starting scale
      call PDF_Param_Iteration(parminuit,iflag)

      if(doCI) CIvarval = parminuit(idxCIval)


      do i=1,ntot
         THEO(i) = 0.d0
         THEO_MOD(i) = 0.d0
      end do

*    Extra constraints on input PDF due to momentum and quark 
*    counting sum rules:
      kflag=0
      if (Itheory.eq.0)  then 
         print '(''-------------------------------------------------'')'
         call SumRules(kflag)
         if (kflag.eq.1) then
            write(6,*) ' --- problem in SumRules, kflag = 1'
            call HF_errlog(12020516,
     +           'F: FCN - problem in SumRules, kflag = 1')
         end if
      end if

*    set alphas
      if(itheory.eq.0) then 
         call setalf(dble(alphas),Mz*Mz)
         alphaSzero= hf_get_alphas(1.0D0)
         call RT_SetAlphaS(alphaSzero)
         if(IPDFSET.eq.5)
     $      call PDFINP(LHAPDFsubr, IPDFSET, dble(0.001), epsi, nwds)
      end if 

*    Call evolution
      call Evolution
   
*    Initialise theory calculation per iteration
      call GetTheoryIteration

*    Calculate theory for datasets:
      do idataset=1,NDATASETS
         if(NDATAPOINTS(idataset).gt.0)
     $      call GetTheoryForDataset(idataset)
      end do

      end


C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_calc_derivatives
      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'fcn.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'
      include 'for_debug.inc'
      include "endmini.inc"
      include "CI.inc"

c    counters      
      integer :: i_par, i_dat, i, idx, II

      character :: ans

      integer :: idxQ2, idxX, idxY, idxS
      double precision :: X(NTOT),Y(NTOT),Q2(NTOT),S(ndatasets)

c    buffers
      double precision :: THEO_buffer(NTOT), RqTrue, Rqerr, Parbuffer

      integer :: needed_pars(mne), nd_prs_end

c    external functions and subroutines
      double precision :: getPar, getParerr
      external :: calc_theo

      integer, parameter :: MNI = 50

      double precision GRD(50),G2(50),GSTEP(50),g_dummy(200),DGRD(50)
      common/MN7DER/ GRD, G2, GSTEP, g_dummy, DGRD
      
      double precision :: parminuit(200), ALIM(200), BLIM(200)
      common/MN7EXT/ parminuit, ALIM, BLIM
      
      character*10 pname(mne)
      common/MN7NAM/ pname
      
      character*3 num

      double precision WERR(50), ern(50), erp(50), globc(50)
      common/MN7ERR/ erp, ern, WERR, globc

      integer nvarl(200), NIOFEX(200), neofix(50)
      common/MN7INX/ nvarl, NIOFEX, neofix

*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      RqTrue = parminuit(idxCIval)
      Rqerr  = WERR(NIOFEX(idxCIval))


      call check_CIvarval
      call derivs_to_zero
      call get_Q2xys

      nd_prs_end = 1
*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      THEO = THEO*0D0
      call calc_theo(g_dummy, parminuit, 4)
      m0 = THEO

      do i_par = 1, mne
         if(NIOFEX(i_par) .eq. 0) cycle
         if(WERR(NIOFEX(i_par)) .eq. 0D0) cycle

         if(i_par .ne. idxCIval) then
         needed_pars(nd_prs_end) = i_par
            nd_prs_end = nd_prs_end + 1
         end if

         print*, " ~ ~ ~ ~ ~ ~ ~ ~ ~ "//TRIM(pname(i_par))//
     $           " ~ ~ ~ ~ ~ ~ ~ ~ ~ "
         print*, "parerr:", WERR(NIOFEX(i_par))

         Parbuffer = parminuit(i_par)

         
         if(i_par .ne. idxCIval) then
            parminuit(i_par) = Parbuffer + 0.5D0*WERR(NIOFEX(i_par))
         else
            parminuit(i_par) = Parbuffer + WERR(NIOFEX(i_par))
         endif
         call calc_theo(g_dummy, parminuit, 4)
         
        THEO_buffer = THEO

         if(i_par .ne. idxCIval) then
            parminuit(i_par) = Parbuffer - 0.5D0*WERR(NIOFEX(i_par))
         else
            parminuit(i_par) = Parbuffer - WERR(NIOFEX(i_par))
         endif
         call calc_theo(g_dummy, parminuit, 4)
         
         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)
               if(i_par .ne. idxCIval) then
                  theta_0(idx, i_par) = (THEO_buffer(idx) - THEO(idx))/
     $                                        WERR(NIOFEX(i_par))
               else if(doCI) then
                  m1(idx) = (THEO_buffer(idx) - THEO(idx))/
     $                                (2D0*Rqerr)
                  m2(idx) = (THEO(idx) + THEO_buffer(idx) -2D0*m0(idx))/
     $                                   (2D0*(Rqerr**2))
               endif
            end do
         end do

         if(i_par .eq. idxCIval .or. .not. doCI) goto 663
c               . . . . . . . . . . . . . . . . . . . . . . 

         parminuit(idxCIval) = RqTrue + Rqerr
         
         parminuit(i_par) = Parbuffer + 0.5D0*WERR(NIOFEX(i_par))
         call calc_theo(g_dummy, parminuit, 4)

         THEO_buffer = THEO

         parminuit(i_par) = Parbuffer - 0.5D0*WERR(NIOFEX(i_par))
         call calc_theo(g_dummy, parminuit, 4)
         
         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)
               theta_p(idx, i_par) = (THEO_buffer(idx) - THEO(idx)) /
     $                                      WERR(NIOFEX(i_par))
            end do
         end do

c               . . . . . . . . . . . . . . . . . . . . . . 

         parminuit(idxCIval) = RqTrue - Rqerr

         parminuit(i_par) = Parbuffer + 0.5D0*WERR(NIOFEX(i_par))
         call calc_theo(g_dummy, parminuit, 4)

         THEO_buffer = THEO

         parminuit(i_par) = Parbuffer - 0.5D0*WERR(NIOFEX(i_par))
         call calc_theo(g_dummy, parminuit, 4)

         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)
               theta_m(idx, i_par) = (THEO_buffer(idx) - THEO(idx)) /
     $                                     WERR(NIOFEX(i_par))
            end do
         end do

 663  continue
  
         parminuit(i_par) = Parbuffer
         parminuit(idxCIval) = RqTrue
      end do
*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --

      if(doCI) then
         theta_1 = (theta_p - theta_m) / (2.D0*Rqerr)
         theta_2 = (theta_p + theta_m - theta_0*2.D0) / (2.D0*Rqerr**2)
      end if


      nd_prs_end = nd_prs_end - 1

      open(103, file = "CIDerivatives.txt", status = 'unknown')
*     - - - - - -
      write(103, '(A, I3)') "Npars ", nd_prs_end
*     - - - - - -
      write(103, '(''CIvar '', F10.7, F10.7)')
     $                    parminuit(idxCIval), WERR(NIOFEX(idxCIval))
*     - - - - - -
      do 6600 i = 1, nd_prs_end
 6600    write(103,'(A10)',advance='no') pname(needed_pars(i))
      write(103,*) ! new line
*     - - - - - -
      write(103, 668, advance='no') "'Data Set'", "'Q2'", "'x'", "'y'"
      do i = 1, 3
         select case (i)
            case(1)
               write(103, '(A25)', advance='no') "'m0'"
               num = "0"
            case(2)
               write(103, '(A25)', advance='no') "'m1'"
               num = "1"
            case(3)
               num = "2"
               write(103, '(A25)', advance='no') "'m2'"
         end select

         do 6604 i_par = 1, nd_prs_end
 6604       write(103, '(A25)', advance='no') "'theta_"//TRIM(num)//
     $                         "_"//TRIM(pname(needed_pars(i_par)))//"'"
      end do
      write(103,*) ! new line
*     - - - - - -
      do i_dat = 1, ndatasets
         do i = 1, NDATAPOINTS(i_dat)
            idx = DATASETIDX(i_dat, i)
            write(103,665,advance='no') TRIM(DATASETLABEL(i_dat))
     $                                , Q2(idx), X(idx), Y(idx), m0(idx)
            do 6601 i_par = 1, nd_prs_end
 6601          write(103,667,advance='no') theta_0(idx,needed_pars(i_par))
            write(103,666,advance='no') m1(idx)
            do 6602 i_par = 1, nd_prs_end
 6602          write(103,667,advance='no') theta_1(idx,needed_pars(i_par))
            write(103,666,advance='no') m2(idx)
            do 6603 i_par = 1, nd_prs_end
 6603          write(103,667,advance='no') theta_2(idx,needed_pars(i_par))
            write(103,*) ! new line
         end do
      end do
*            'Data Set' 'Q2'  'x'    'y'
 668  format(    A15,   A10,  A25,   A25) !    m0
 665  format(    A15,  F10.2,F25.20,F25.20, E25.15E3)
 666  format(E25.15E3)
 667  format(E25.15E3)

      close(103)

      return

*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      contains
      
      subroutine derivs_to_zero
      implicit none
      do i = 1, NTOT
         m0(i) = 0D0
         m1(i) = 0D0
         m2(i) = 0D0
         do i_par = 1, mne
            theta_0(i, i_par) = 0D0
            theta_1(i, i_par) = 0D0
            theta_2(i, i_par) = 0D0
         enddo
      enddo
      end subroutine derivs_to_zero

      subroutine check_CIvarval
      implicit none
      if(RqTrue .ne. 0d0) then
         print*,"Warning: CIvarval is not zero! Continue anyway? (y/n):"
 
 701     continue
         read *, ans
         
         if(ans .eq. 'n') then
            call hf_stop
         else if(ans .ne. 'y') then
            print *, "Please enter 'y' or 'n':"
            goto 701
         endif
      end if
      end subroutine check_CIvarval

      subroutine get_Q2xys
      implicit none
      integer GetBinIndex, GetInfoIndex
      do i_dat = 1, ndatasets
         idxQ2 = GetBinIndex(i_dat, 'Q2')
         idxX  = GetBinIndex(i_dat, 'x')
         idxY  = GetBinIndex(i_dat, 'y')
            
         if(idxY .eq. 0) then
            idxS = GetInfoIndex(i_dat, 'sqrt(S)')
            if(idxS .gt. 0) then
               S(i_dat) = (DATASETInfo(GetInfoIndex(i_dat, 'sqrt(S)')
     $                                                       ,i_dat))**2
            endif
         endif

         do i = 1, NDATAPOINTS(i_dat)
            idx = DATASETIDX(i_dat, i)

            X(idx)  = AbstractBins(idxX, idx)
            Q2(idx) = AbstractBins(idxQ2, idx)
            if(idxY .eq. 0) then
               Y(idx) = Q2(idx) / ( X(idx) * S(i_dat) )
            else
               Y(idx) = AbstractBins(idxY, idx)
            endif
         end do
      end do
      end subroutine get_Q2xys

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_parameters(read_pars, read_errors)
      include "d506dp.inc"
      include "d506cm.inc"

      logical, intent(in) :: read_pars, read_errors

      integer parnum, i
      character*10 parname
      real parval, parerr
      external CI_read_Rq, mnpout

      integer i_int
      character*32 mn_pname
      double precision lolim, uplim

   
      if(read_pars) then
         do i = 1, mne
            U(i) = 0D0
            if(i .gt. 50) cycle
            WERR(i) = 0D0
         end do
      end if


      open(103, file = "CISimpFitData.txt", status = 'old', err = 2631)
      do
         read(103, *, end = 2661) parnum, parname, parval, parerr
         if(TRIM(parname) .eq. "CI_Rq") cycle

         if(read_pars) U(parnum) = parval
         if(read_errors.and.parerr.ne.0D0.and.NIOFEX(parnumm).ne.0) 
     $      WERR(NIOFEX(parnum))=parerr

      end do

 2661 continue
      close(103)

      call CI_read_CIvar

      return


 2631 print *, "Error: file 'CISimpFitData.txt' not found."
      close(103)
      call hf_stop

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_missing_parerrors
      include "d506dp.inc"
      include "d506cm.inc"

      integer parnum
      character*10 parname
      double precision parval, parerr
      external CI_try_read_CIvarerr


      open(103, file = "CISimpFitData.txt", status = 'old', err = 2632)
      do
         read(103, *, end = 2662) parnum, parname, parval, parerr
         if(NIOFEX(parnum)       .eq.  0      ) cycle
         if(WERR(NIOFEX(parnum)) .ne.  0d0    ) cycle
         if(TRIM(parname)        .eq.  "CI_Rq") cycle

         WERR(NIOFEX(parnum)) = parerr
      end do

      goto 2666

 2632 print *, "CI: file 'CISimpFitData.txt' not found."
      goto 2666
 2662 print *, "CI: error while reading 'CISimpFitData.txt'."

   
 2666 close(103)

      call CI_try_read_CIvarerr

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_CIvar
      implicit none

      include "steering.inc"

      integer parnum
      character*10 parname
      double precision parval
      double precision parerr

      double precision WERR(50), ern(50), erp(50), globc(50)
      common/MN7ERR/ erp, ern, WERR, globc

      integer nvarl(200), NIOFEX(200), neofix(50)
      common/MN7INX/ nvarl, NIOFEX, neofix

      double precision :: parminuit(200), ALIM(200), BLIM(200)
      common/MN7EXT/ parminuit, ALIM, BLIM

      external setParValue, setParError

      if(.not. doCI) return

      open(104, file = "CIval_in.txt", status = 'old', err = 3012)
      read(104, *) parnum, parname, parval, parerr
      parminuit(idxCIval) = parval
      WERR(NIOFEX(idxCIval)) =  parerr
      close(104)
      return

 3012 print *, "Error: file 'CIval_in.txt' not found."
      close(104)
      call hf_stop

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_try_read_CIvarerr
      implicit none

      include "steering.inc"

      integer parnum
      character*10 parname
      double precision parval, parerr
      
      double precision WERR(50), ern(50), erp(50), globc(50)
      common/MN7ERR/ erp, ern, WERR, globc

      integer nvarl(200), NIOFEX(200), neofix(50)
      common/MN7INX/ nvarl, NIOFEX, neofix

      double precision :: parminuit(200), ALIM(200), BLIM(200)
      common/MN7EXT/ parminuit, ALIM, BLIM


      open(103, file = "CIval_in.txt", status = 'old', err = 3023)
      read(103, *) parnum, parname, parval, parerr
      WERR(NIOFEX(idxCIval)) = parerr
      goto 3063

 3023  continue
      print *, "CI: file 'CIval_in.txt' not found."
 3063  continue
      close(103)

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine simpfcn(g_dummy, chi2out, parminuit, iflag)
      implicit none

      include 'ntot.inc'
      include 'steering.inc'
      include 'datasets.inc'
      include 'indata.inc'
      include 'theo.inc'
      include 'fcn.inc'
      include 'couplings.inc'
      include 'qcdnumhelper.inc'
      include 'for_debug.inc'
      include "endmini.inc"
      include "CI.inc"

      double precision, intent(in)  :: g_dummy(*), parminuit(*)
      integer         , intent(in)  :: iflag
      double precision, intent(out) :: chi2out

c    counters      
      integer i, i_par, i_dat, idx

c    external functions and subroutines
      external CI_read_missing_parerrors
      double precision chi2data_theory

      double precision WERR(50), ern(50), erp(50), globc(50)
      common/MN7ERR/ erp, ern, WERR, globc

      integer nvarl(200), NIOFEX(200), neofix(50)
      common/MN7INX/ nvarl, NIOFEX, neofix
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      print '('' --- simpfcn...'')'

      THEO = THEO*0.
      do i = 1, NTOT
         THEO(i) = m0(i)
      end do

      do i_par = 1, mne
         if(NIOFEX(i_par) .eq. 0) cycle

         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)

               if(i_par .eq. idxCIval) then
                  THEO(idx) = THEO(idx) +
     $                        m1(idx)*parminuit(idxCIval) +
     $                        m2(idx)*parminuit(idxCIval)**2
               else
                  THEO(idx) = THEO(idx) + 
     $               (
     $                  theta_0(idx,i_par) +
     $                  theta_1(idx,i_par)*parminuit(idxCIval) +
     $                  theta_2(idx,i_par)*parminuit(idxCIval)**2
     $               ) * (parminuit(i_par) - p0(i_par))
               end if
            end do
         end do
      end do
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --

c    get Chi2
      chi2out = chi2data_theory(iflag)
      print *, "chi2out = ", chi2out
      print *, ' '
      print *, ' '

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_derivatives
      implicit none
   
      include "ntot.inc"
      include "endmini.inc"
      include "datasets.inc"
      include "CI.inc"
      character*20 dsname, trash

      character*10 pnames(mne)
      double precision Q2, X, Y
      integer II, i_dat, idx, i, i_par, npars, prespars(mne)

      external hf_stop

      print '('' - - - - - - - - - - reading derivatives... '//
     $         ' - - - - - - - - - - '')'

      call derivs_to_zero
*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      open(103, file = "CIDerivatives.txt", status = 'old', err = 1301)

      read(103, *) trash, npars

c    second line in trash
      read(103,*) trash

      read(103,*,end=1304,err=1304) pnames(1:npars)

c    fourth line in trash
      read(103,*) trash

      do 1305 i = 1, npars
 1305    prespars(i) = name_to_idx(pnames(i))
      
      do i_dat = 1, ndatasets
         do i = 1, NDATAPOINTS(i_dat)
            idx = DATASETIDX(i_dat, i)
            read(103,1303,end=1302,err=1304,advance='no') dsname, Q2, X
     $                                                  , Y
            
            read(103,'(E25.15E3)',end=1302,err=1304,advance='no')m0(idx)
            
            do 1306 i_par = 1, npars
 1306          read(103,'(E25.15E3)',end=1302,err=1304,advance='no')
     $            theta_0(idx, prespars(i_par))

            read(103,'(E25.15E3)',end=1302,err=1304,advance='no')m1(idx)
            
            do 1307 i_par = 1, npars
 1307          read(103,'(E25.15E3)',end=1302,err=1304,advance='no')
     $            theta_1(idx, prespars(i_par))
            
            read(103,'(E25.15E3)',end=1302,err=1304,advance='no')m2(idx)
            
            do 1308 i_par = 1, npars - 1
 1308          read(103,'(E25.15E3)',end=1302,err=1304,advance='no') 
     $            theta_2(idx, prespars(i_par))
            read(103,'(E25.15E3)',end=1302,err=1304)
     $         theta_2(idx, prespars(npars))
            
         end do
      end do

      close(103)
      return
*           dsname  Q2      X        Y     m0
 1303 format(A15, F10.2, F25.20, F25.20)

 1302 continue
      print *, "Error: unexpected end of file ""CIDerivatives.txt"""
      call hf_stop

 1301 continue
      print *, "Error: can't read file ""CIDerivatives.txt"""
      call hf_stop

 1304 continue
      print *, "Error while reading ""CIDerivatives.txt"""
      do i = 1, Npars
         print *, "theta_0: ", theta_0(1,prespars(i))
      end do
      call hf_stop


      contains
      
      integer function name_to_idx(pname)
      implicit none

      character*10, intent(in) :: pname
      integer :: i
      character*10 allnames(mne)
      common/MN7NAM/ allnames
      
      external hf_stop

      do i = 1, mne
         if(TRIM(pname) .eq. TRIM(allnames(i))) then
            name_to_idx = i
            return
         end if
      end do
      
      print *, "Error while reading ""CIDerivatives.txt"""
      call hf_stop
      end function name_to_idx
      
      subroutine derivs_to_zero
      implicit none
      integer :: i, i_par

      do i = 1, NTOT
         m0(i) = 0D0
         m1(i) = 0D0
         m2(i) = 0D0
         do i_par = 1, mne
            theta_0(i, i_par) = 0D0
            theta_1(i, i_par) = 0D0
            theta_2(i, i_par) = 0D0
         enddo
      enddo
      end subroutine derivs_to_zero

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_p0
      include "d506dp.inc"
      include "d506cm.inc"
      include "CI.inc"

      integer i
      double precision U_save(MNE)
      external CI_read_parameters

      U_save = U
      do i = 1, MNE
         U(i) = 0.
         p0(i) = 0.
      enddo

      call CI_read_parameters(.true.,.false.)
      p0 = U
      U = U_save

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_simp_fit_data
      implicit none
      external CI_read_derivatives, CI_read_p0
      call CI_read_derivatives
      call CI_read_p0
      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_scan_parerrors
      include "d506dp.inc"
      include "d506cm.inc"
      integer i
      double precision err
      logical fail


      fail = .true.

      do i = 1, 200
         if(NIOFEX(i) .eq. 0) cycle

         err = WERR(NIOFEX(i))
         if(err .ne. 0D0) then
            print *, "par #", i, ", error:", err
            fail = .false.
         end if
      end do

      if(fail) print*, "...fail"

      end

