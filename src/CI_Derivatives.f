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
      double precision function getPar(i_par)
      include "d506dp.inc"
      include "d506cm.inc"
      
      integer i_par

      getPar = U(i_par)
      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision function getParerr(i_par)
      include "d506dp.inc"
      include "d506cm.inc"

      integer i_par

      getParerr = WERR(NIOFEX(i_par))
      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine setParValue(i_par, val)
      include "d506dp.inc"
      include "d506cm.inc"
      
      integer i_par
      double precision val

      U(i_par) = val
      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine setParError(i_par, err)
      include "d506dp.inc"
      include "d506cm.inc"
      
      integer i_par
      double precision err

      WERR(NIOFEX(i_par)) = err
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
      double precision :: dm1dp_p, dm1dp_m

c    external functions and subroutines
      double precision :: getPar, getParerr
      external :: calc_theo

      integer, parameter :: MNI = 50
      double precision GRD(MNI),G2(MNI),GSTEP(MNI),g_dummy(MNE),DGRD(MNI)
      common/MN7DER/ GRD, G2, GSTEP, g_dummy, DGRD
      double precision :: parminuit(MNE), ALIM(MNE), BLIM(MNE)
      common/MN7EXT/ parminuit, ALIM, BLIM


      RqTrue = getPar(idxCIval)
      Rqerr  = getParerr(idxCIval)


      call check_CIvarval
      call derivs_to_zero
      call get_Q2xys


      if(.not.doCI) goto 664 ! skip coutings
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      do i_par = 1, mne
         if(getParerr(i_par) .eq. 0.D0) cycle

         Parbuffer = getPar(i_par)

c                          dm0/dp, dm2/dp, m1, m2
         
         if(i_par .ne. idxCIval) then
            parminuit(i_par) = Parbuffer + 0.5*getParerr(i_par)
         else
            parminuit(i_par) = Parbuffer + getParerr(i_par)
         endif
         call calc_theo(g_dummy, parminuit, 4)
         
         do i = 1, NTOT
            THEO_buffer(i) = THEO(i)
         end do

         if(i_par .ne. idxCIval) then
            parminuit(i_par) = Parbuffer - 0.5*getParerr(i_par)
         else
            parminuit(i_par) = Parbuffer - getParerr(i_par)
         endif
         call calc_theo(g_dummy, parminuit, 4)
         
         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)
               if(i_par .ne. idxCIval) then
                  dm0dp(idx,i_par) = (THEO(idx) - THEO_buffer(idx))/
     $                                      getParerr(i_par)
                  if(.not. doCI) cycle
     
                  dm2dp(idx, i_par) = 2.*(THEO_buffer(idx) - THEO(idx))
               else if(doCI) then
                  m1(idx) = (THEO(idx) - THEO_buffer(idx))/
     $                          (2.*getParerr(i_par))
                  m2(idx) = (THEO(idx) + THEO_buffer(idx) - 2.*m0(idx))/
     $                                   (2.*Rqerr**2)
               endif
            end do
         end do

         if(.not. doCI) cycle
         if(i_par .eq. idxCIval) goto 663
c               ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
c                                dm1/dp, dm2/dp

c          /*  dm1/dp = (dm1dp_p - dm1dp_m) / dp,                *
c           *  where dm1dp_p = (m(p+, Rq+) - m(p+, Rq-)) / dRq,  *
c           *        dm1dp_m = (m(p-, Rq+) - m(p-, Rq-)) / dRq.  */

c                              dm1dp_p

         parminuit(i_par) = Parbuffer + 0.5*getParerr(i_par)
         
         parminuit(idxCIval) = RqTrue + Rqerr
         call calc_theo(g_dummy, parminuit, 4)

         do i = 1, NTOT
            THEO_buffer(i) = THEO(i)
         end do

         parminuit(idxCIval) = RqTrue - Rqerr
         call calc_theo(g_dummy, parminuit, 4)
         
         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)
               dm1dp_p = (THEO(idx) - THEO_buffer(idx))/(2.*Rqerr)
     
               dm2dp(idx, i_par) = dm2dp(idx, i_par) + THEO(idx)
     $                           + THEO_buffer(idx)
            end do
         end do

c               . . . . . . . . . . . . . . . . . . . . . . 
c                            dm1dp_m

         parminuit(i_par) = Parbuffer - getParerr(i_par)

         parminuit(idxCIval) = RqTrue + Rqerr
         call calc_theo(g_dummy, parminuit, 4)

         do i = 1, NTOT
            THEO_buffer(i) = THEO(i)
         end do

         parminuit(idxCIval) = RqTrue - Rqerr
         call calc_theo(g_dummy, parminuit, 4)

         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)
               dm1dp_m = (THEO(idx) - THEO_buffer(idx))/(2.*Rqerr)
               dm1dp(idx, i_par) = (dm1dp_p - dm1dp_m)/
     $                               getParerr(i_par)

               dm2dp(idx,i_par)=dm2dp(idx,i_par)-(THEO(idx)+THEO_buffer(idx))/
     $                              (2. * Rqerr**2 * getParerr(i_par))
            end do
         end do

 663  continue
  
         parminuit(i_par) = Parbuffer
         parminuit(idxCIval) = RqTrue
      end do
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --

 664  continue


      open(103, file = "CIDerivatives.txt", status = 'unknown')
      write(103, '(A, I3)') "Npars ", mne
      write(103, '(''CIvar '', F10.7, F10.7)')
     $                    getPar(idxCIval), getParerr(idxCIval)
      do i_dat = 1, ndatasets
         do i = 1, NDATAPOINTS(i_dat)
            idx = DATASETIDX(i_dat, i)
            write(103,665,advance='no') TRIM(DATASETLABEL(i_dat))
     $                           , Q2(idx), X(idx), Y(idx), m0(idx) 
            write(103,667,advance='no') dm0dp(idx, : )
            write(103,666,advance='no') m1(idx)
            write(103,667,advance='no') dm1dp(idx, : )
            write(103,666,advance='no') m2(idx)
            write(103,667,advance='no') dm2dp(idx, : )
            write(103,*)
         end do
      end do
 665  format(A15, F10.2, F25.20, F25.20, E25.15E3)
 666  format(E25.15E3)
 667  format(200E25.15E3)
      close(103)

      return

      
      contains
      
      subroutine derivs_to_zero
      implicit none
      do i = 1, NTOT
         m0(i) = THEO(i)

         m1(i) = 0D0
         m2(i) = 0D0
         do i_par = 1, mne
            dm0dp(i, i_par) = 0D0
            dm1dp(i, i_par) = 0D0
            dm2dp(i, i_par) = 0D0
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
               S(i_dat)=(DATASETInfo(GetInfoIndex(i_dat,'sqrt(S)'),i_dat))**2
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

      integer parnum
      character*10 parname
      double precision parval, parerr
      external CI_read_Rq


      open(103, file = "CISimpFitData.txt", status = 'old', err = 263)
      do
         read(103, *, end = 266) parnum, parname, parval, parerr
         if(TRIM(parname) .eq. "CI_Rq") cycle

         if(read_pars)     U(parnum) = parval
         if(read_errors)   WERR(NIOFEX(parnum)) = parerr
      end do

 266  continue
      close(103)

      call CI_read_Rq

      return


 263  print *, "Error: file 'CISimpFitData.txt' not found."
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
      external CI_try_read_Rq


      open(103, file = "CISimpFitData.txt", status = 'old', err = 263)
      do
         read(103, *, end = 266) parnum, parname, parval, parerr
         if(WERR(NIOFEX(parnum)) .ne. 0d0) cycle
         if(TRIM(parname) .eq. "CI_Rq" .and.WERR(NIOFEX(parnum)).eq.0d0)
     $      call CI_try_read_Rq

         WERR(NIOFEX(parnum)) = parerr
      end do

 263  print *, "CI: file 'CISimpFitData.txt' not found."
 266  continue
      close(103)

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_Rq
      implicit none

      include "steering.inc"

      integer parnum
      character*10 parname
      double precision parval, parerr


      if(.not.doCI) then
         call setParValue(idxCIval, 0D0)
         call setParError(idxCIval, 0D0)
         return
      endif


      open(103, file = "CIval_in.txt", status = 'old', err = 301)
      read(103, *) parnum, parname, parval, parerr
      call setParValue(idxCIval, parval)
      call setParError(idxCIval, parerr)
      goto 306

 301  continue
      print *, "Error: file 'CIval_in.txt' not found."
      close(103)
      call hf_stop

 306  continue
      close(103)

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_try_read_Rq
      implicit none

      include "steering.inc"

      integer parnum
      character*10 parname
      double precision parval, parerr


      open(103, file = "CIval_in.txt", status = 'old', err = 301)
      read(103, *) parnum, parname, parval, parerr
      call setParError(idxCIval, parerr)
      goto 306

 301  continue
      print *, "CI: file 'CIval_in.txt' not found."
 306  continue
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
      double precision getPar, getParerr
      external CI_read_missing_parerrors
      double precision chi2data_theory

c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      print '('' --- simpfcn...'')'


      THEO = m0

      do i_par = 1, mne
         do i_dat = 1, ndatasets
            do i = 1, NDATAPOINTS(i_dat)
               idx = DATASETIDX(i_dat, i)

               if(i_par .eq. idxCIval) then
                  THEO(idx) = THEO(idx) +
     $                            m1(idx)*getPar(idxCIval) +
     $                            m2(idx)*getPar(idxCIval)**2
               else
                  THEO(idx) = THEO(idx) + 
     $               (
     $                  dm0dp(idx,i_par) +
     $                  dm1dp(idx,i_par)*getPar(idxCIval) +
     $                  dm2dp(idx,i_par)*getPar(idxCIval)**2
     $               ) * (getPar(i_par) - p0(i_par))
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
      character*15 name, trash
      double precision Q2, X, Y
      integer II, i_dat, idx, i

      external hf_stop

      print '('' - - - - - - - - - - reading derivatives... '//
     $         ' - - - - - - - - - - '')'
*     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - --
      open(103, file = "CIDerivatives.txt", status = 'old', err = 1301)

c    first two lines in trash
      read(103,*) (trash, i=1,2)
      
      do i_dat = 1, ndatasets
         do i = 1, NDATAPOINTS(i_dat)
            idx = DATASETIDX(i_dat, i)

            read(103,1303,end=1302,err=1304) name,   Q2,   X,   Y
     $                                     , m0(idx), dm0dp(idx,:)
     $                                     , m1(idx), dm1dp(idx,:)
     $                                     , m2(idx), dm2dp(idx,:)
         end do
      end do

      close(103)
      return
*            Name  Q2      X        Y
 1303 format(A15, F10.2, F25.20, F25.20,
     $  E25.15E3,200E25.15E3,E25.15E3,200E25.15E3,E25.15E3,200E25.15E3)
*          m0       dm0dp       m1        dm1dp      m2       dm2dp

 1302 continue
      print *, "Error: unexpected end of file ""CIDerivatives.txt"""
      call hf_stop

 1301 continue
      print *, "Error: can't read file ""CIDerivatives.txt"""
      call hf_stop

 1304 continue
      print *, "Error while reading ""CIDerivatives.txt"""
      call hf_stop

      end

C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine CI_read_p0
      include "d506dp.inc"
      include "d506cm.inc"
      include "CI.inc"

      integer i, Ulen
      parameter(Ulen = size(U))
      double precision U_save(Ulen)
      external CI_read_parameters

      U_save = U
      do i = 1, Ulen
         U(i) = 0
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
         err = WERR(NIOFEX(i))
         if(err .ne. 0) then
            print *, "par #", i, ", error:", err
            fail = .false.
         end if
      end do

      if(fail) print*, "...fail"

      end

