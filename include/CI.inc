      logical CI_simpfcn_mode
      common/CIrules/ CI_simpfcn_mode


      logical CI_extra_reading_done
      common/CIstatics/ CI_extra_reading_done


      double precision m0(2500), dm0dp(2500, 200)
     $               , m1(2500), dm1dp(2500, 200)
     $               , m2(2500), dm2dp(2500, 200)
      
      double precision p0(200)
      
      common/CISimpFitData/ m0, dm0dp, m1, dm1dp, m2, dm2dp, p0
