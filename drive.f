      PROGRAM NEKTON
C
C       VERSION 2.7
C
C------------------------------------------------------------------------
C
C       3-D Isoparametric Legendre Spectral Element Solver.
C
C
C
C       N a v i e r     S t o k e s     S o l v e r 
C
C
C       LAPLACIAN FORMULATION
C
C       The program solves the equations 
C
C       (1) ro ( dv/dt + (v*grad)v ) = mu*div(grad(v)) - grad(p) + f 
C       (2) div(v) = 0,
C
C       subject to periodic, Dirichlet or Neumann boundary 
C       conmonditions for the velocity v.
C
C       STRESS FORMULATION
C
C       (1) ro*( dv/dt + (v*grad)v ) = mu*( div grad v + grad div v )
C                                      - grad p + f 
C       (2) div(v) = 0,
C
C       subject to velocity and traction boundary conditions.
C
C
C       In both formulations, the velocity v is solved on a 
C       Gauss-Legendre Lobatto spectral element mesh, while the
C       pressure is solved on a Gauss-Legendre mesh (staggered mesh).
C
C
C
C       P a s s i v e     S c a l a r     S o l v e r
C
C
C       The program solves the equation
C
C       (3)  rhocp ( dT/dt + v*grad T ) = k div(grad T) + q 
C
C       subject to Dirichlet and Neumann boundary conditions
C       for the passive scalar T.
C
C--------------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'DEALIAS'
      include 'DOMAIN'
      include 'ZPER'
c
      include 'OPCTR'
      include 'CTIMER'

C     Declare scratch arrays
C     NOTE: no initial declaration needed. Linker will take 
c           care about the size of the CBs
c
c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)
  
      REAL e, oe
      integer WDS
      REAL*8 t0,tp

      call iniproc !  processor initialization 
      if (nid.eq.0) write(6,*) 'Number of Processors ::',np

      TIME0  = dnekclock()
      etimes = dnekclock()
      ISTEP  = 0
      tpp    = 0.0

      call opcount(1)

      call initdim
C
C     Data initialization
C
      call initdat
      call files
      t0 = dnekclock()

      call readat  ! Read processor map, followed by data.
      if (nid.eq.0) write(6,*) 'readat time ::',dnekclock()-t0,
     &                         ' seconds'

      call setvar  ! initialize some variables
      call echopar ! echo back the parameter stack

c     Check for zero steps

      instep=1
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0


C     Geometry initialization

      igeom = 2

      call connect
      call genwz
      call usrdat
      call gengeom (igeom)
      call usrdat2
      call geom_reset(1)    ! recompute Jacobians, etc.

      if (nid.eq.0) write(6,*) 'NELGV/NELGT/NX1:',nelgv,nelgt,nx1
      if (nid.eq.0) write(6,*) nid,' call vrdsms',instep

      call vrdsmsh  ! verify mesh topology


C     Field initialization
      if (nid.eq.0) write(6,*) nid,' call setlog',instep
      call setlog

      call bcmask
C
C     Need eigenvalues to set tolerances in presolve (SETICS)
C
      if (nid.eq.0) write(6,*) 'this is ifflow:',ifflow,nsteps
      if (fintim.ne.0.0.or.nsteps.ne.0) call geneig (igeom)
      call vrdsmsh
C
C     Solver initialization  (NOTE:  Uses "SOLN" space as scratch...)

      if (ifflow.and.nsteps.gt.0) then

         if (ifsplit) then
            call set_up_h1_crs
         else

            call estrat

            if (nid.eq.0) 
     $         write(6,*) nid,'estrat',iesolv,nsteps,fintim,solver_type

            if (fintim.ne.0.0.or.nsteps.ne.0) then
               if (iftran.and.solver_type.eq.'itr') then
                  call set_overlap
               elseif (solver_type.eq.'fdm'.or.solver_type.eq.'pdm')then
                  call gfdm_init
               elseif (solver_type.eq.'25D') then
                  call g25d_init
               endif
            endif

         endif
      endif


      call usrdat3

C     The properties are set if PRESOLVE is used in SETICS,
C     otherwise they are set in the beginning of the time stepping loop

      call setics

      CALL SETPROP

      CALL USERCHK
      CALL COMMENT
      CALL SSTEST (ISSS) 

      CALL TIME00
      CALL opcount(2)
      CALL dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow

      DO 1000 ISTEP=1,NSTEPS

         IF (IFTRAN) CALL SETTIME
         if (ifmhd ) call cfl_check

         CALL SETSOLV
         CALL COMMENT

         if (ifsplit) then

            if (ifheat)      call heat     (0)

                             call setprop
                             call qthermal

            if (ifflow)      call fluid    (0)

         else

           call setprop
           do igeom=1,2

              if (ifgeom) then
                 call gengeom (igeom)
                 call geneig  (igeom)
              endif

              if (ifmhd) then
                             call induct   (igeom)
                 if (ifheat) call heat     (igeom)
              else
                 if (ifheat)      call heat     (igeom)
                 if (ifflow)      call fluid    (igeom)
                 if (ifmvbd)      call meshv    (igeom)
              endif

              if (ifpert) then
                 if (ifflow) call fluidp   (igeom)
                 if (ifheat) call heatp    (igeom)
              endif

           enddo

         endif

         if (.not.ifmhd) then  ! (filter in induct.f for ifmhd)
             if (param(103).gt.0) alpha_filt=param(103)
             if (param(103).gt.0) call q_filter(alpha_filt)
         endif

         call prepost (.false.,'his')
         call userchk

         if (lastep .eq. 1) goto 1001
 1000 CONTINUE
 1001 CONTINUE
c
      call opcount(3)
      call timeout
C
C----------------------------------------------------------------------
C     Time stepping loop end
C----------------------------------------------------------------------
C
c     call prepost (.true.,'his')
c
      if (instep.eq.0) then
         lastep=1
         t0 = dnekclock()
         call prepost (.true.,'his')
         tpp = tpp + (dnekclock()-t0)
         nsteps=0
         call userchk
      endif
C
      if (nid.eq.0) then
         write(6,*) 'prepost time ::',tpp,' seconds'
      endif
C
      CALL COMMENT
      CALL DIAGNOS
      call crs_stats(xxth)
      call exitt
      END
C
      subroutine initdim
C-------------------------------------------------------------------
C
C     Transfer array dimensions to common
C
C-------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
C
      NX1=LX1
      NY1=LY1
      NZ1=LZ1
C
      NX2=LX2
      NY2=LY2
      NZ2=LZ2
C
      NX3=LX3
      NY3=LY3
      NZ3=LZ3
C
      NELT=LELT
      NELV=LELV
      NDIM=LDIM
C
      RETURN
      END
C
      subroutine initdat
C--------------------------------------------------------------------
C
C     Initialize and set default values.
C
C--------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      COMMON /DOIT/ IFDOIT
      LOGICAL       IFDOIT
C
C     Set default logicals
C
      IFFLOW  = .FALSE.
      IFMVBD  = .FALSE.
      IFHEAT  = .TRUE.
      IFSPLIT = .FALSE.
      IFDOIT  = .FALSE.
      ifxxt   = .false.
      IFCVODE = .false.

      if (lx1.eq.lx2) ifsplit=.true.


C     Turn off (on) diagnostics for communication
C
      IFGPRNT= .FALSE.
C
      CALL RZERO (PARAM,200)
C
C     The initialization of CBC is done in READAT
C
C      LCBC = 3*6*LELT*(LDIMT1+1)
C      CALL BLANK(CBC,LCBC)
C
      CALL BLANK(CCURVE ,8*LELT)
      NEL8 = 8*LELT
      CALL RZERO(XC,NEL8)
      CALL RZERO(YC,NEL8)
      CALL RZERO(ZC,NEL8)
C
      NTOT=NX1*NY1*NZ1*LELT
      CALL RZERO(ABX1,NTOT)
      CALL RZERO(ABX2,NTOT)
      CALL RZERO(ABY1,NTOT)
      CALL RZERO(ABY2,NTOT)
      CALL RZERO(ABZ1,NTOT)
      CALL RZERO(ABZ2,NTOT)
      CALL RZERO(VGRADT1,NTOT)
      CALL RZERO(VGRADT2,NTOT)

      RETURN
      END
C
      subroutine comment
C---------------------------------------------------------------------
C
C     No need to comment !!
C
C---------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP'
      LOGICAL  IFCOUR
      SAVE     IFCOUR
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      REAL*8 ETIME0,ETIME1,ETIME2
      SAVE   ETIME0,ETIME1,ETIME2
      DATA   ETIME0,ETIME1,ETIME2 /0.0, 0.0, 0.0/
      REAL*8 DNEKCLOCK
C
C     Only node zero makes comments.
      IF (NID.NE.0) RETURN
C
C
      IF (ETIME0.EQ.0.0) ETIME0=DNEKCLOCK()
      ETIME1=ETIME2
      ETIME2=DNEKCLOCK()
C
      IF (ISTEP.EQ.0) THEN
         IFCOUR  = .FALSE.
         DO 10 IFIELD=1,NFIELD
            IF (IFADVC(IFIELD)) IFCOUR = .TRUE.
 10      CONTINUE
         IF (IFWCNO) IFCOUR = .TRUE.
         WRITE (6,*) ' '
         WRITE (6,*) 'Initialization successfully completed'
         IF (TIME.NE.0.0) WRITE (6,*) 'Initial time is:',TIME
         WRITE (6,*) ' '
         WRITE (6,*) 'START OF SIMULATION'
         WRITE (6,*) ' '
      ELSEIF (ISTEP.GT.0 .AND. LASTEP.EQ.0 .AND. IFTRAN) THEN
         ETIME=ETIME2-ETIME1
         TTIME=ETIME2-ETIME0
         IF (     IFCOUR) 
     $      WRITE (6,100) ISTEP,TIME,DT,COURNO/10,TTIME,ETIME
         IF (.NOT.IFCOUR) WRITE (6,101) ISTEP,TIME,DT
      ELSEIF (LASTEP.EQ.1) THEN
         WRITE (6,*) ' '
         WRITE (6,*) 'Simulation successfully completed'
      ENDIF
 100  FORMAT('Step',I6,', t=',1pE14.7,', DT=',1pE14.7
     $,', C=',F7.3,2(1pE11.4))
 101  FORMAT('Step',I6,', time=',1pE12.5,', DT=',1pE11.3)
C      call flush_io()
      RETURN
      END
C
      subroutine exit2
C     This is here because calling Sun-4 Fortran's EXIT causes a core dump (!)
      call exitt
      END
C
      subroutine setvar
C------------------------------------------------------------------------
C
C     Initialize variables
C
C------------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'DEALIAS'
      include 'TSTEP'
C
C     Enforce splitting/Uzawa according to the way the code was compiled
C
c
      nxd = lxd
      nyd = lyd
      nzd = lzd
C
C     Geometry on Mesh 3 or 1?
C
      IFGMSH3 = .TRUE.
      IF ( IFSTRS )           IFGMSH3 = .FALSE.
      IF (.NOT.IFFLOW)        IFGMSH3 = .FALSE.
      IF ( IFSPLIT )          IFGMSH3 = .FALSE.
C
      NFIELD = 1
      IF (IFHEAT) THEN
         NFIELD = 2 + NPSCAL
         NFLDTM = 1 + NPSCAL
      ENDIF
c
      nfldt = nfield
      if (ifmhd) then
         nfldt  = nfield + 1
         nfldtm = nfldtm + 1
      endif
c
      IF (IFMODEL) CALL SETTMC
      IF (IFMODEL.AND.IFKEPS) THEN
         NPSCAL = 1
         NFLDTM = NPSCAL + 1
         IF (LDIMT.LT.NFLDTM) THEN
            WRITE (6,*) 'k-e turbulence model activated'
            WRITE (6,*) 'Insufficient number of field arrays'
            WRITE (6,*) 'Rerun through PRE or change SIZE file'
            call exitt
         ENDIF
         NFIELD = NFIELD + 2
         CALL SETTURB
      ENDIF
      MFIELD = 1
      IF (IFMVBD) MFIELD = 0
C
      DO 100 IFIELD=MFIELD,nfldt
         IF (IFTMSH(IFIELD)) THEN
             NELFLD(IFIELD) = NELT
         ELSE
             NELFLD(IFIELD) = NELV
         ENDIF
 100  CONTINUE
C
      NMXH   = 1000 !  1000
      NMXP   = 2000 !  2000
      NMXE   = 100 !  1000
      NMXNL  = 10  !  100
C
      PARAM(86) = 0
C
      BETAG  = PARAM(3)
      GTHETA = PARAM(4)
      DT     = abs(PARAM(12))
      DTINIT = DT
      FINTIM = PARAM(10)
      NSTEPS = PARAM(11)
      IOCOMM = PARAM(13)
      TIMEIO = PARAM(14)
      IOSTEP = PARAM(15)
      LASTEP = 0
      TOLPDF = abs(PARAM(21))
      TOLHDF = abs(PARAM(22))
      TOLREL = abs(PARAM(24))
      TOLABS = abs(PARAM(25))
      CTARG  = PARAM(26)
      NBDINP = PARAM(27)
      NABMSH = PARAM(28)

      if(abs(PARAM(16)).eq.2) IFCVODE = .true.
      if(abs(PARAM(16)).eq.3) IFEXPL = .true.

C
C     Check accuracy requested.
C
      IF (TOLREL.LE.0.) TOLREL = 0.01
C
C     Relaxed pressure iteration; maximum decrease in the residual.
C
      PRELAX = 0.1*TOLREL
      IF (.NOT.IFTRAN .AND. .NOT.IFNAV) PRELAX = 1.E-5
C
C     Tolerance for nonlinear iteration
C
      TOLNL  = 1.E-4
C
C     Fintim overrides nsteps
C
      IF (FINTIM.NE.0.) NSTEPS = 1000000000
      IF (.NOT.IFTRAN ) NSTEPS = 1
C
C     Print interval defaults to 1
C
      IF (IOCOMM.EQ.0)  IOCOMM = nsteps+1

C
C     Set logical for Boussinesq approx (natural convection)
C
      IFNATC = .FALSE.
      IF (BETAG.GT.0.) IFNATC=.TRUE.
      IF(IFLOMACH) IFNATC = .FALSE.
C
C     Set default for mesh integration scheme
C
      IF (NABMSH.LE.0 .OR. NABMSH.GT.3) THEN
         NABMSH    = NBDINP
         PARAM(28) = (NABMSH)
      ENDIF
C
C     Set default for mixing length factor
C
      TLFAC = 0.14
      IF (PARAM(49) .LE. 0.0) PARAM(49) = TLFAC
C
C     Courant number only applicable if convection in ANY field.
C
      IADV  = 0
      IFLD1 = 1
      IF (.NOT.IFFLOW) IFLD1 = 2
      DO 200 IFIELD=IFLD1,nfldt
         IF (IFADVC(IFIELD)) IADV = 1
 200  CONTINUE
C
C     If characteristics, need number of sub-timesteps (DT/DS).
C     Current sub-timeintegration scheme: RK4.
C     If not characteristics, i.e. standard semi-implicit scheme,
C     check user-defined Courant number.
C
      IF (IADV.EQ.1) CALL SETCHAR
C
C     Initialize order of time-stepping scheme (BD)
C     Initialize time step array.
C
      NBD    = 0
      CALL RZERO (DTLAG,10)
C
C     Useful constants
C
      one = 1.
      PI  = 4.*ATAN(one)
C
      RETURN
      END
C
      subroutine echopar
C
C     Echo the nonzero parameters from the readfile to the logfile
C
      include 'SIZE'
      include 'INPUT'
      CHARACTER*80 STRING
      CHARACTER*1  STRING1(80)
      EQUIVALENCE (STRING,STRING1)
C
      IF (nid.ne.0) RETURN
C
      OPEN (UNIT=9,FILE=REAFLE,STATUS='OLD')
      REWIND(UNIT=9)
C
C
      READ(9,*,ERR=400)
      READ(9,*,ERR=400) VNEKTON
      NKTONV=VNEKTON
      VNEKMIN=2.5
      IF(VNEKTON.LT.VNEKMIN)THEN
         PRINT*,' Error: This NEKTON Solver Requires a .rea file'
         PRINT*,' from prenek version ',VNEKMIN,' or higher'
         PRINT*,' Please run the session through the preprocessor'
         PRINT*,' to bring the .rea file up to date.'
         call exitt
      ENDIF
      READ(9,*,ERR=400) NDIM
c     error check
      IF(NDIM.NE.LDIM)THEN
         WRITE(6,10) LDIM,NDIM
   10       FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $              /,2X,'       for spatial dimension equal to',I2,'.'
     $              /,2X,'       The data file has dimension',I2,'.')
         CALL exitt
      ENDIF
C
      CALL BLANK(STRING,80)
      CALL CHCOPY(STRING,REAFLE,80)
      Ls=LTRUNC(STRING,80)
      READ(9,*,ERR=400) NPARAM
      WRITE(6,82) NPARAM,(STRING1(j),j=1,Ls)
C
      DO 20 I=1,NPARAM
         CALL BLANK(STRING,80)
         READ(9,80,ERR=400) STRING
         Ls=LTRUNC(STRING,80)
         IF (PARAM(i).ne.0.0) WRITE(6,81) I,(STRING1(j),j=1,Ls)
   20 CONTINUE
   80 FORMAT(A80) 
   81 FORMAT(I4,3X,80A1)
   82 FORMAT(I4,3X,'Parameters from file:',80A1)
      CLOSE (UNIT=9)

      if(param(2).ne.param(8).and.nid.eq.0) then
         write(6,*) 'Note VISCOS not equal to CONDUCT!'
         write(6,*) 'Note VISCOS  =',PARAM(2)
         write(6,*) 'Note CONDUCT =',PARAM(8)
      endif
c
      return
C
C     Error handling:
C
  400 CONTINUE
      WRITE(6,401)
  401 FORMAT(2X,'ERROR READING PARAMETER DATA'
     $    ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
      CALL exitt
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING LOGICAL DATA'
     $    ,/,2X,'ABORTING IN ROUTINE ECHOPAR.')
      CALL exitt
C
      RETURN
      END
C
      subroutine gengeom (igeom)
C----------------------------------------------------------------------
C
C     Generate geometry data
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'GEOM'
      include 'WZ'
C
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C
      IF (IGEOM.EQ.1) THEN
         RETURN
      ELSEIF (IGEOM.EQ.2) THEN
         CALL LAGMASS
         IF (ISTEP.EQ.0) CALL GENCOOR (XM3,YM3,ZM3)
         IF (ISTEP.GE.1) CALL UPDCOOR
         CALL GEOM1 (XM3,YM3,ZM3)
         CALL GEOM2
         CALL UPDMSYS (1)
         CALL VOLUME
         CALL SETINVM
         CALL SETDEF
         CALL SFASTAX
         IF (ISTEP.GE.1) CALL EINIT
      ELSEIF (IGEOM.EQ.3) THEN
c
c        Take direct stiffness avg of mesh
c
         ifieldo = ifield
         CALL GENCOOR (XM3,YM3,ZM3)
         if (ifheat) then
            ifield = 2
            CALL dssum(xm3,nx3,ny3,nz3)
            call col2 (xm3,tmult,ntot3)
            CALL dssum(ym3,nx3,ny3,nz3)
            call col2 (ym3,tmult,ntot3)
            if (if3d) then
               CALL dssum(xm3,nx3,ny3,nz3)
               call col2 (xm3,tmult,ntot3)
            endif
         else
            ifield = 1
            CALL dssum(xm3,nx3,ny3,nz3)
            call col2 (xm3,vmult,ntot3)
            CALL dssum(ym3,nx3,ny3,nz3)
            call col2 (ym3,vmult,ntot3)
            if (if3d) then
               CALL dssum(xm3,nx3,ny3,nz3)
               call col2 (xm3,vmult,ntot3)
            endif
         endif
         CALL GEOM1 (XM3,YM3,ZM3)
         CALL GEOM2
         CALL UPDMSYS (1)
         CALL VOLUME
         CALL SETINVM
         CALL SETDEF
         CALL SFASTAX
         ifield = ifieldo
      ENDIF
C
      RETURN
      END
C
      subroutine files
C----------------------------------------------------------------------
C
C     Defines machine specific input and output file names.
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
C
      CHARACTER*132 NAME
      CHARACTER*1   SESS1(132),PATH1(132),NAM1(132)
      EQUIVALENCE  (SESSION,SESS1)
      EQUIVALENCE  (PATH,PATH1)
      EQUIVALENCE  (NAME,NAM1)
      CHARACTER*1  DMP(4),FLD(4),REA(4),HIS(4),SCH(4) ,ORE(4), NRE(4)
      CHARACTER*1  RE2(4)
      CHARACTER*4  DMP4  ,FLD4  ,REA4  ,HIS4  ,SCH4   ,ORE4  , NRE4
      CHARACTER*4  RE24  
      EQUIVALENCE (DMP,DMP4), (FLD,FLD4), (REA,REA4), (HIS,HIS4)
     $          , (SCH,SCH4), (ORE,ORE4), (NRE,NRE4)
     $          , (RE2,RE24)
      DATA DMP4,FLD4,REA4 /'.dmp','.fld','.rea'/
      DATA HIS4,SCH4      /'.his','.sch'/
      DATA ORE4,NRE4      /'.ore','.nre'/
      DATA RE24           /'.re2'       /
      CHARACTER*78  STRING
C
C     Find out the session name:
C
      CALL BLANK(SESSION,132)
      CALL BLANK(PATH   ,132)
      OPEN (UNIT=8,FILE='SESSION.NAME',STATUS='OLD')
      READ(8,10) SESSION
      READ(8,10) PATH
      CLOSE(UNIT=8)
   10 FORMAT(A132)


      CALL BLANK(REAFLE,132)
      CALL BLANK(RE2FLE,132)
      CALL BLANK(FLDFLE,132)
      CALL BLANK(HISFLE,132)
      CALL BLANK(SCHFLE,132)
      CALL BLANK(DMPFLE,132)
      CALL BLANK(OREFLE,132)
      CALL BLANK(NREFLE,132)
      CALL BLANK(NAME  ,132)
C
C     Construct file names containing full path to host:
C
      LS=LTRUNC(SESSION,132)
      LPP=LTRUNC(PATH,132)
      LSP=LS+LPP
c
      call chcopy(nam1(    1),path1,lpp)
      call chcopy(nam1(lpp+1),sess1,ls )
      l1 = lpp+ls+1
      ln = lpp+ls+4
c
c
c .rea file
      call chcopy(nam1  (l1),rea , 4)
      call chcopy(reafle    ,nam1,ln)
c      write(6,*) 'reafile:',reafle
c
c .re2 file
      call chcopy(nam1  (l1),re2 , 4)
      call chcopy(re2fle    ,nam1,ln)
c
c .fld file
      call chcopy(nam1  (l1),fld , 4)
      call chcopy(fldfle    ,nam1,ln)
c
c .his file
      call chcopy(nam1  (l1),his , 4)
      call chcopy(hisfle    ,nam1,ln)
c
c .sch file
      call chcopy(nam1  (l1),sch , 4)
      call chcopy(schfle    ,nam1,ln)
c
c
c .dmp file
      call chcopy(nam1  (l1),dmp , 4)
      call chcopy(dmpfle    ,nam1,ln)
c
c .ore file
      call chcopy(nam1  (l1),ore , 4)
      call chcopy(orefle    ,nam1,ln)
c
c .nre file
      call chcopy(nam1  (l1),nre , 4)
      call chcopy(nrefle    ,nam1,ln)
c
C     Write the name of the .rea file to the logfile.
C
      IF (NID.EQ.0) THEN
         CALL CHCOPY(STRING,REAFLE,78)
         WRITE(6,1000) STRING
         WRITE(6,1001) 
 1000    FORMAT(//,2X,'Beginning session:',/,2X,A78)
 1001    FORMAT(/,' ')
      ENDIF
C
      RETURN
      END
C
      subroutine settime
C----------------------------------------------------------------------
C
C     Store old time steps and compute new time step, time and timef.
C     Set time-dependent coefficients in time-stepping schemes.
C
C----------------------------------------------------------------------
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      COMMON  /CPRINT/ IFPRINT
      LOGICAL          IFPRINT
      SAVE
C
      irst = param(46)
C
C     Set time step.
C
      DO 10 ILAG=10,2,-1
         DTLAG(ILAG) = DTLAG(ILAG-1)
 10   CONTINUE
      CALL SETDT
      DTLAG(1) = DT
      IF (ISTEP.EQ.1 .and. irst.le.0) DTLAG(2) = DT
C
C     Set time.
C
      TIMEF    = TIME
      TIME     = TIME+DT
C
C     Set coefficients in AB/BD-schemes.
C
      CALL SETORDBD
      if (irst.gt.0) nbd = nbdinp
      CALL RZERO (BD,10)
      CALL SETBD (BD,DTLAG,NBD)
      NAB = 3
      IF (ISTEP.LE.2 .and. irst.le.0) NAB = ISTEP
      CALL RZERO   (AB,10)
      CALL SETABBD (AB,DTLAG,NAB,NBD)
      IF (IFMVBD) THEN
         NBDMSH = 1
         NABMSH = PARAM(28)
         IF (NABMSH.GT.ISTEP .and. irst.le.0) NABMSH = ISTEP
         IF (IFSURT)          NABMSH = NBD
         CALL RZERO   (ABMSH,10)
         CALL SETABBD (ABMSH,DTLAG,NABMSH,NBDMSH)
      ENDIF
C
C     Set logical for printout to screen/log-file
C
      IFPRINT = .FALSE.
      IF (IOCOMM.GT.0.AND.MOD(ISTEP,IOCOMM).EQ.0) IFPRINT=.TRUE.
      IF (ISTEP.eq.1  .or. ISTEP.eq.0           ) IFPRINT=.TRUE.
C
      RETURN
      END
C
C
      subroutine geneig (igeom)
C-----------------------------------------------------------------------
C
C     Compute eigenvalues. 
C     Used for automatic setting of tolerances and to find critical
C     time step for explicit mode. 
C     Currently eigenvalues are computed only for the velocity mesh.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'EIGEN'
      include 'INPUT'
      include 'TSTEP'
C
      IF (IGEOM.EQ.1) RETURN
C
C     Decide which eigenvalues to be computed.
C
      IF (IFFLOW) THEN
C
         IFAA  = .FALSE.
         IFAE  = .FALSE.
         IFAS  = .FALSE.
         IFAST = .FALSE.
         IFGA  = .TRUE.
         IFGE  = .FALSE.
         IFGS  = .FALSE.
         IFGST = .FALSE.
C
C        For now, only compute eigenvalues during initialization.
C        For deforming geometries the eigenvalues should be 
C        computed every time step (based on old eigenvectors => more memory)
C
         IMESH  = 1
         IFIELD = 1
         TOLEV  = 1.E-3
         TOLHE  = TOLHDF
         TOLHR  = TOLHDF
         TOLHS  = TOLHDF
         TOLPS  = TOLPDF
         CALL EIGENV
         CALL ESTEIG
C
      ELSEIF (IFHEAT.AND..NOT.IFFLOW) THEN
C
         CALL ESTEIG
C
      ENDIF
C
      RETURN
      END
C-----------------------------------------------------------------------
      subroutine fluid (igeom)
C
C     Driver for solving the incompressible Navier-Stokes equations.
C
C     Current version:
C     (1) Velocity/stress formulation.
C     (2) Constant/variable properties.
C     (3) Implicit/explicit time stepping.
C     (4) Automatic setting of tolerances .
C     (5) Lagrangian/"Eulerian"(operator splitting) modes
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'

      ifield = 1
      imesh  = 1
      call unorm
      call settolv

      if (ifsplit) then

c        PLAN 4: TOMBO SPLITTING
c                - Time-dependent Navier-Stokes calculation (Re>>1).
c                - Same approximation spaces for pressure and velocity.
c                - Weakly compressible (div u .ne. 0).

         call plan4
         call twalluz (igeom) ! Turbulence model
         call chkptol         ! check pressure tolerance
         call vol_flow        ! check for fixed flow rate

      elseif (iftran) then

         if (param(181).eq.0) then   !  Same as PLAN 1 w/o nested iteration
            call plan3 (igeom)       !  Std. NEKTON time stepper  !
         else
            call plan1 (igeom)
         endif

         if (ifmodel)    call twalluz (igeom) ! Turbulence model
         if (igeom.eq.2) call chkptol         ! check pressure tolerance
         if (igeom.eq.2) call vol_flow        ! check for fixed flow rate

      else   !  steady Stokes, non-split

c             - Steady/Unsteady Stokes/Navier-Stokes calculation.
c             - Consistent approximation spaces for velocity and pressure.
c             - Explicit treatment of the convection term. 
c             - Velocity/stress formulation.

         call plan1 (igeom) ! The NEKTON "Classic".

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine heat (igeom)
C
C     Driver for temperature or passive scalar.
C
C     Current version:
C     (1) Varaiable properties.
C     (2) Implicit time stepping.
C     (3) User specified tolerance for the Helmholtz solver
C         (not based on eigenvalues).
C     (4) A passive scalar can be defined on either the 
C         temperatur or the velocity mesh.
C     (5) A passive scalar has its own multiplicity (B.C.).  
C
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'TURBO'

      if (ifcvode) then

         call cdscal_cvode

      elseif (ifsplit) then

         do igeo=1,2
         do ifield=2,nfield
            intype        = -1
            if (.not.iftmsh(ifield)) imesh = 1
            if (     iftmsh(ifield)) imesh = 2
            call unorm
            call settolt
            call cdscal (igeo)
         enddo
         enddo

      else  ! PN-PN-2

         do ifield=2,nfield
            intype        = -1
            if (.not.iftmsh(ifield)) imesh = 1
            if (     iftmsh(ifield)) imesh = 2
            call unorm
            call settolt
            call cdscal (igeom)
         enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine meshv (igeom)

C     Driver for mesh velocity used in conjunction with moving geometry.
C
C-----------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
C
      IF (IGEOM.EQ.1) RETURN
C
      IFIELD = 0
      NEL    = NELFLD(IFIELD)
      IMESH  = 1
      IF (IFTMSH(IFIELD)) IMESH = 2
C
      CALL UPDMSYS (0)
      CALL MVBDRY  (NEL)
      CALL ELASOLV (NEL)
C
      RETURN
      END
      subroutine rescont (ind)
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
c
      if (np.gt.1) return
      irst = param(46)
      iwrf = 1
      if (irst.ne.0) then
         iwrf = mod(istep,iabs(irst))
         if (lastep.eq.1) iwrf = 0
      endif
c
      if (ind.eq.1 .and. irst.gt.0) call rstartc (ind)
      if (ind.eq.0 .and. iwrf.eq.0) call rstartc (ind)
c
      return     
      end
      subroutine rstartc (ind)
c
      include 'SIZE'
      include 'TOTAL'
      common /SCRSF/ xm3(lx3,ly3,lz3,lelv)
     $             , ym3(lx3,ly3,lz3,lelv)
     $             , zm3(lx3,ly3,lz3,lelv)
c
      integer icall1,icall2
      save    icall1,icall2
      data    icall1 /0/
      data    icall2 /0/
c
      if (np.gt.1) return
      ntov1=nx1*ny1*nz1*nelv
      ntov2=nx2*ny2*nz2*nelv
      ntot1=nx1*ny1*nz1*nelt
      ntfc1=nx1*nz1*6*nelv 
      ntow1=lx1m*ly1m*lz1m*nelfld(0)
      ntoe1=lx1m*ly1m*lz1m*nelv
      ntotf=ntot1*ldimt
      nlag =lorder-1
c
      if (ind.eq.1) then
c
          iru=22
          if (icall1.eq.0) then
             icall1=1
             open(unit=22,file=orefle,status='OLD')
          endif
c
          rewind iru
          read(iru,1100,end=9000) time,dt,courno,(dtlag(i),i=1,10)
          read(iru,1100,end=9000) eigaa, eigas, eigast, eigae,
     $                            eigga, eiggs, eiggst, eigge
c
          write (6,*) '  '
          write (6,*) 'READ RESTART FILE, TIME =',time
c
          iread = 0
          if (ifmvbd) then
             read (iru,1100,end=9000) (xm1(i,1,1,1),i=1,ntot1)
             read (iru,1100,end=9000) (ym1(i,1,1,1),i=1,ntot1)
             if (ndim.eq.3)
     $       read (iru,1100,end=9000) (zm1(i,1,1,1),i=1,ntot1)
             read (iru,1100,end=9000) (wx(i,1,1,1) ,i=1,ntow1)
             read (iru,1100,end=9000) (wy(i,1,1,1) ,i=1,ntow1)
             if (ndim.eq.3)
     $       read (iru,1100,end=9000) (wz(i,1,1,1) ,i=1,ntow1)
            if (nlag.ge.1) then
             read (iru,1100,end=9000) (wxlag(i,1,1,1,1) ,i=1,ntow1*nlag)
             read (iru,1100,end=9000) (wylag(i,1,1,1,1) ,i=1,ntow1*nlag)
             if (ndim.eq.3)
     $       read (iru,1100,end=9000) (wzlag(i,1,1,1,1) ,i=1,ntow1*nlag)
            endif
          endif
c
          iread = 1
          if (ifflow) then
             read (iru,1100,end=9000) (vx(i,1,1,1) ,i=1,ntov1)
             read (iru,1100,end=9000) (vy(i,1,1,1) ,i=1,ntov1)
             if (ndim.eq.3)
     $       read (iru,1100,end=9000) (vz(i,1,1,1) ,i=1,ntov1)
             read (iru,1100,end=9000) (pr(i,1,1,1) ,i=1,ntov2)
             read (iru,1100,end=9000) (abx2(i,1,1,1),i=1,ntov1)
             read (iru,1100,end=9000) (aby2(i,1,1,1),i=1,ntov1)
             if (ndim.eq.3)
     $       read (iru,1100,end=9000) (abz2(i,1,1,1),i=1,ntov1)
             read (iru,1100,end=9000) (abx1(i,1,1,1),i=1,ntov1)
             read (iru,1100,end=9000) (aby1(i,1,1,1),i=1,ntov1)
             if (ndim.eq.3)
     $       read (iru,1100,end=9000) (abz1(i,1,1,1),i=1,ntov1)
            if (nlag.ge.1) then
             read (iru,1100,end=9000) (vxlag (i,1,1,1,1),i=1,ntov1*nlag)
             read (iru,1100,end=9000) (vylag (i,1,1,1,1),i=1,ntov1*nlag)
             if (ndim.eq.3)
     $       read (iru,1100,end=9000) (vzlag (i,1,1,1,1),i=1,ntov1*nlag)
             read (iru,1100,end=9000) (bm1lag(i,1,1,1,1),i=1,ntot1*nlag)
            endif
          endif
c
          iread = 2
          if (ifheat) then
             read (iru,1100,end=9000) (t(i,1,1,1,1),i=1,ntotf)
             read (iru,1100,end=9000) (vgradt1(i,1,1,1,1),i=1,ntotf)
             read (iru,1100,end=9000) (vgradt2(i,1,1,1,1),i=1,ntotf)
            if (nlag.ge.1) then
             read (iru,1100,end=9000) (tlag(i,1,1,1,1,1),i=1,ntotf*nlag)
            endif
          endif
c
          iread = 3
          if (ifmodel .and. .not.ifkeps) then
             read (iru,1100,end=9000) tlmax,tlimul
             read (iru,1100,end=9000) (turbl(i,1,1,1),i=1,ntov1)
          endif
          if (ifcwuz) then
             read (iru,1100,end=9000) (zwall (i,1,1,1),i=1,ntfc1)
             read (iru,1100,end=9000) (uwall(i,1,1,1),i=1,ntfc1)
          endif 
c
          if (ifgeom) then
             call geom1 (xm3,ym3,zm3)
             call geom2
             call updmsys (1)
             call volume
             call setinvm
          endif
c
      elseif (ind.eq.0) then
c
          iwu=23
          if (icall2.eq.0) then
             icall2=1
             open(unit=23,file=nrefle,status='NEW')
          endif
c
          rewind iwu
          write(iwu,1100) time,dt,courno,(dtlag(i),i=1,10)
          write(iwu,1100) eigaa, eigas, eigast, eigae,
     $                    eigga, eiggs, eiggst, eigge
c
          write (6,*) '  '
          write (6,*) 'WRITE RESTART FILE, TIME =',time
c
          if (ifmvbd) then
             write (iwu,1100) (xm1(i,1,1,1),i=1,ntot1)
             write (iwu,1100) (ym1(i,1,1,1),i=1,ntot1)
             if (ndim.eq.3)
     $       write (iwu,1100) (zm1(i,1,1,1),i=1,ntot1)
             write (iwu,1100) (wx(i,1,1,1) ,i=1,ntow1)
             write (iwu,1100) (wy(i,1,1,1) ,i=1,ntow1)
             if (ndim.eq.3)
     $       write (iwu,1100) (wz(i,1,1,1) ,i=1,ntow1)
            if (nlag.ge.1) then
             write (iwu,1100) (wxlag(i,1,1,1,1) ,i=1,ntow1*nlag)
             write (iwu,1100) (wylag(i,1,1,1,1) ,i=1,ntow1*nlag)
             if (ndim.eq.3)
     $       write (iwu,1100) (wzlag(i,1,1,1,1) ,i=1,ntow1*nlag)
            endif
          endif
c
          if (ifflow) then
             write (iwu,1100) (vx(i,1,1,1) ,i=1,ntov1)
             write (iwu,1100) (vy(i,1,1,1) ,i=1,ntov1)
             if (ndim.eq.3)
     $       write (iwu,1100) (vz(i,1,1,1) ,i=1,ntov1)
             write (iwu,1100) (pr(i,1,1,1) ,i=1,ntov2)
             write (iwu,1100) (abx2(i,1,1,1),i=1,ntov1)
             write (iwu,1100) (aby2(i,1,1,1),i=1,ntov1)
             if (ndim.eq.3)
     $       write (iwu,1100) (abz2(i,1,1,1),i=1,ntov1)
             write (iwu,1100) (abx1(i,1,1,1),i=1,ntov1)
             write (iwu,1100) (aby1(i,1,1,1),i=1,ntov1)
             if (ndim.eq.3)
     $       write (iwu,1100) (abz1(i,1,1,1),i=1,ntov1)
            if (nlag.ge.1) then
             write (iwu,1100) (vxlag (i,1,1,1,1),i=1,ntov1*nlag)
             write (iwu,1100) (vylag (i,1,1,1,1),i=1,ntov1*nlag)
             if (ndim.eq.3)
     $       write (iwu,1100) (vzlag (i,1,1,1,1),i=1,ntov1*nlag)
             write (iwu,1100) (bm1lag(i,1,1,1,1),i=1,ntot1*nlag)
            endif
          endif
c
          if (ifheat) then
             write (iwu,1100) (t(i,1,1,1,1),i=1,ntotf)
             write (iwu,1100) (vgradt1(i,1,1,1,1),i=1,ntotf)
             write (iwu,1100) (vgradt2(i,1,1,1,1),i=1,ntotf)
            if (nlag.ge.1) then
             write (iwu,1100) (tlag(i,1,1,1,1,1),i=1,ntotf*nlag)
            endif
          endif
c
          if (ifmodel .and. .not.ifkeps) then
             write (iwu,1100) tlmax,tlimul
             write (iwu,1100) (turbl(i,1,1,1),i=1,ntov1)
          endif
          if (ifcwuz) then
             write (iwu,1100) (zwall(i,1,1,1),i=1,ntfc1)
             write (iwu,1100) (uwall(i,1,1,1),i=1,ntfc1)
          endif 
c
      endif
c
      return
c
 1100 format ((5e16.8))
 9000 continue
c     
      write ( 6,*)  ' RECORD OUT-OF-ORDER DURING READING OF RESTART' 
      write ( 6,*)  ' FILE -- iread =',iread
c
      call exitt
      end
      subroutine time00
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
C
      nmxmf=0
      nmxms=0
      ndsum=0
      nvdss=0
      nsett=0
      ncdtp=0
      npres=0
      nmltd=0
      ngsum=0
      nprep=0
      ndsnd=0
      ndadd=0
      nhmhz=0
      naxhm=0
      ngop =0
      nusbc=0
      ncopy=0
      ninvc=0
      ninv3=0
      nsolv=0
      nslvb=0
      nddsl=0
      ncrsl=0
      ndott=0
      nbsol=0
c
      tmxmf=0.0
      tmxms=0.0
      tdsum=0.0
      tvdss=0.0
      tvdss=0.0
      tdsmn=9.9e9
      tdsmx=0.0
      tsett=0.0
      tcdtp=0.0
      tpres=0.0
      teslv=0.0
      tmltd=0.0
      tgsum=0.0
      tgsmn=9.9e9
      tgsmx=0.0
      tprep=0.0
      tdsnd=0.0
      tdadd=0.0
      thmhz=0.0
      taxhm=0.0
      tgop =0.0
      tusbc=0.0
      tcopy=0.0
      tinvc=0.0
      tinv3=0.0
      tsolv=0.0
      tslvb=0.0
      tddsl=0.0
      tcrsl=0.0
      tdott=0.0
      tbsol=0.0
      tbso2=0.0
      etims0= dnekclock()
C
      return
      end
C
      subroutine timeout
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      dimension vdsum(lp),vgsum(lp),vbsol(lp),vusbc(lp)
      dimension vvdss(lp),vdadd(lp),work(lp),vgop (lp)
      dimension vdsmn(lp),vdsmx(lp),vgsmn(lp),vgsmx(lp)
C
      real min_dsum, max_dsum, avg_dsum
      real min_vdss, max_vdss, avg_vdss
      real min_gop,  max_gop,  avg_gop
      real min_crsl, max_crsl, avg_crsl
c
      real dhc, dwork
C
      tstop=dnekclock()
      ttotal=tstop-etimes
      tttstp=tstop-etims0
c
c
      min_vdss = tvdss
      call gop(min_vdss,wwork,'m  ',1)
      max_vdss = tvdss
      call gop(max_vdss,wwork,'M  ',1)
      avg_vdss = tvdss
      call gop(avg_vdss,wwork,'+  ',1)
      avg_vdss = avg_vdss/np
c
      min_dsum = tdsum
      call gop(min_dsum,wwork,'m  ',1)
      max_dsum = tdsum
      call gop(max_dsum,wwork,'M  ',1)
      avg_dsum = tdsum
      call gop(avg_dsum,wwork,'+  ',1)
      avg_dsum = avg_dsum/np
c
      min_gop = tgop
      call gop(min_gop,wwork,'m  ',1)
      max_gop = tgop
      call gop(max_gop,wwork,'M  ',1)
      avg_gop = tgop
      call gop(avg_gop,wwork,'+  ',1)
      avg_gop = avg_gop/np
c
      min_crsl = tcrsl
      call gop(min_crsl,wwork,'m  ',1)
      max_crsl = tcrsl
      call gop(max_crsl,wwork,'M  ',1)
      avg_crsl = tcrsl
      call gop(avg_crsl,wwork,'+  ',1)
      avg_crsl = avg_crsl/np
c
      tttstp = tttstp + 1e-7
      if (nid.eq.0) then
         write(6,*) 'total time',ttotal,tttstp
         ttotal=tttstp
         pcopy=tcopy/tttstp
         write(6,*) 'copy time',ncopy,tcopy,pcopy
         pmxmf=tmxmf/tttstp
         write(6,*) 'mxmf time',nmxmf,tmxmf,pmxmf
         pmxms=tmxms/tttstp
         write(6,*) 'mxms time',nmxms,tmxms,pmxms
         pinv3=tinv3/tttstp
         write(6,*) 'inv3 time',ninv3,tinv3,pinv3
         pinvc=tinvc/tttstp
         write(6,*) 'invc time',ninvc,tinvc,pinvc
         pmltd=tmltd/tttstp
         write(6,*) 'mltd time',nmltd,tmltd,pmltd
         pcdtp=tcdtp/tttstp
         write(6,*) 'cdtp time',ncdtp,tcdtp,pcdtp
         peslv=teslv/tttstp 
         write(6,*) 'eslv time',neslv,teslv,peslv
         ppres=tpres/tttstp
         write(6,*) 'pres time',npres,tpres,ppres
         phmhz=thmhz/tttstp
         write(6,*) 'hmhz time',nhmhz,thmhz,phmhz
         pusbc=tusbc/tttstp
         write(6,*) 'usbc time',nusbc,tusbc,pusbc
         paxhm=taxhm/tttstp
         write(6,*) 'axhm time',naxhm,taxhm,paxhm
c
         pgop =tgop /tttstp
         write(6,*) 'gop  time',ngop ,tgop ,pgop 
         write(6,*) 'gop  min ',min_gop 
         write(6,*) 'gop  max ',max_gop 
         write(6,*) 'gop  avg ',avg_gop 
c
         pvdss=tvdss/tttstp
         write(6,*) 'vdss time',nvdss,tvdss,pvdss
         write(6,*) 'vdss min ',min_vdss
         write(6,*) 'vdss max ',max_vdss
         write(6,*) 'vdss avg ',avg_vdss
c
         pdsum=tdsum/tttstp
         write(6,*) 'dsum time',ndsum,tdsum,pdsum
         write(6,*) 'dsum min ',min_dsum
         write(6,*) 'dsum max ',max_dsum
         write(6,*) 'dsum avg ',avg_dsum
c
         pgsum=tgsum/tttstp
         write(6,*) 'gsum time',ngsum,tgsum,pgsum
         pdsnd=tdsnd/tttstp
         write(6,*) 'dsnd time',ndsnd,tdsnd,pdsnd
         pdadd=tdadd/tttstp
         write(6,*) 'dadd time',ndadd,tdadd,pdadd
         pdsmx=tdsmx/tttstp
         write(6,*) 'dsmx time',ndsmx,tdsmx,pdsmx
         pdsmn=tdsmn/tttstp
         write(6,*) 'dsmn time',ndsmn,tdsmn,pdsmn
         pgsmx=tgsmx/tttstp
         write(6,*) 'gsmx time',ngsmx,tgsmx,pgsmx
         pgsmn=tgsmn/tttstp
         write(6,*) 'gsmn time',ngsmn,tgsmn,pgsmn
         pslvb=tslvb/tttstp
         write(6,*) 'slvb time',nslvb,tslvb,pslvb
         pddsl=tddsl/tttstp
         write(6,*) 'ddsl time',nddsl,tddsl,pddsl
c
         pcrsl=tcrsl/tttstp
         write(6,*) 'crsl time',ncrsl,tcrsl,pcrsl
         write(6,*) 'crsl min ',min_crsl
         write(6,*) 'crsl max ',max_crsl
         write(6,*) 'crsl avg ',avg_crsl
c
         psolv=tsolv/tttstp
         write(6,*) 'solv time',nsolv,tsolv,psolv
         psett=tsett/tttstp
         write(6,*) 'sett time',nsett,tsett,psett
         pprep=tprep/tttstp
         write(6,*) 'prep time',nprep,tprep,pprep
         pbsol=tbsol/tttstp
         write(6,*) 'bsol time',nbsol,tbsol,pbsol
         pbso2=tbso2/tttstp
         write(6,*) 'bso2 time',nbso2,tbso2,pbso2
      endif
      if (np.gt.0) then
c        call rzero(vbsol,np)
c        vbsol(node)=tbsol
c        call gop(vbsol,work,'+  ',np)
c
         call rzero(vusbc,np)
         vusbc(node)=tusbc
         call gop(vusbc,work,'+  ',np)
         write(6,*) nid,' nusbc',nusbc,pusbc
c
         call rzero(vvdss,np)
         vvdss(node)=tdsnd
         call gop(vvdss,work,'+  ',np)
c
         call rzero(vdadd,np)
         vdadd(node)=tdadd
         call gop(vdadd,work,'+  ',np)
c
         call rzero(vgsum,np)
         vgsum(node)=tgsum
         call gop(vgsum,work,'+  ',np)
c
         call rzero(vvdss,np)
         vvdss(node)=tvdss
         call gop(vvdss,work,'+  ',np)
c
         call rzero(vdsum,np)
         vdsum(node)=tdsum
         call gop(vdsum,work,'+  ',np)
c
         call rzero(vgop ,np)
         vgop (node)=tgop 
         call gop(vgop ,work,'+  ',np)
c
         call rzero(vdsmx,np)
         vdsmx(node)=tdsmx
         call gop(vdsmx,work,'+  ',np)
c
         call rzero(vdsmn,np)
         vdsmn(node)=tdsmn
         call gop(vdsmn,work,'+  ',np)
c
         call rzero(vgsmx,np)
         vgsmx(node)=tgsmx
         call gop(vgsmx,work,'+  ',np)
c
         call rzero(vgsmn,np)
         vgsmn(node)=tgsmn
         call gop(vgsmn,work,'+  ',np)
c
         ndsum = max(ndsum,1)
         nvdss = max(nvdss,1)
         if (nid.eq.0) then
c
            write(6,202) np,nelgv,tttstp
            write(6,203) ndsum,nvdss,nbsol
            write(6,*) 'qqq ip tdsum tdsnd tdadd tgsum tgop tusbs'
            do 100 ip=1,np
               write(6,204) ip,vdsum(ip),vvdss(ip),vdadd(ip)
     $                         ,vgsum(ip),vgop(ip),vusbc(ip)
  100       continue
            write(6,*) 'qqq ip dsavg tdsmn tdsmx tgsmn tgsmx'
            do 200 ip=1,np
               dsavg = vdsum(ip)/(ndsum)
               write(6,204) ip,dsavg,vdsmn(ip),vdsmx(ip)
     $                        ,vgsmn(ip),vgsmx(ip)
  200       continue
  202       format('qqq  np,nel,tttstp:',2i8,f12.5)
  203       format('qqq  num procs',/,' dot,dsum,bsol:',3i8)
  204       format('qqq', i7,6f12.5)
C
C
            write(6,*) 'qqq ip tdsum tdsnd tdadd tgsum tgop'
            do 110 ip=1,np
               rdsum=vdsum(ip)/tttstp
               rdsnd=vvdss(ip)/tttstp
               rdadd=vdadd(ip)/tttstp
               rgsum=vgsum(ip)/tttstp
               rgop =vgop (ip)/tttstp
               rusbc=vusbc(ip)/tttstp
               write(6,204) ip,rdsum,rdsnd,rdadd,rgsum,rgop,rusbc
  110       continue
         endif
      endif
C
      return
      end

      subroutine opcount(ICALL)
C
      include 'SIZE'
      include 'OPCTR'
      character*6 sname(maxrts)
      integer     ind  (maxrts)
      integer     idum (maxrts)
C
      if (icall.eq.1) then
         nrout=0
      endif
      if (icall.eq.1.or.icall.eq.2) then
         dcount = 0.0
         do 100 i=1,maxrts
            ncall(i) = 0
            dct(i)   = 0.0
  100    continue
      endif
      if (icall.eq.3) then
C
C        Sort and print out diagnostics
C
         write(6,*) nid,' opcount',dcount
         dhc = dcount
         call gop(dhc,dwork,'+  ',1)
         if (nid.eq.0) then
            write(6,*) nid,' TOTAL OPCOUNT',dhc
         endif
C
         CALL DRCOPY(rct,dct,nrout)
         CALL SORT(rct,ind,nrout)
         CALL CHSWAPR(rname,6,ind,nrout,sname)
         call iswap(ncall,ind,nrout,idum)
C
         if (nid.eq.0) then
            do 200 i=1,nrout
               write(6,201) nid,rname(i),rct(i),ncall(i)
  200       continue
  201       format(2x,' opnode',i4,2x,a6,g18.7,i12)
         endif
      endif
      return
      end
C
      subroutine dofcnt
      include 'SIZE'
      include 'TOTAL'
      COMMON /CTMP0/ DUMMY0(LCTMP0)
      COMMON /CTMP1/ DUMMY1(LCTMP1)
      COMMON /SCRNS/ WORK(LCTMP1)
C
      ntot1=nx1*ny1*nz1*nelv
      ntot2=nx2*ny2*nz2*nelv
C
      if (ifflow) then
         call col3 (work,vmult,v1mask,ntot1)
      else
         call col3 (work,tmult,tmask,ntot1)
      endif
      vpts = glsum(work,ntot1) + .1
      nvtot=vpts
      work(1)=ntot2
      ppts = glsum(work,1) + .1
      nptot=ppts
C
      work(1)=0.0
      do 10 i=1,ntot1
         if (vmult(i,1,1,1).lt.0.5) work(1)=work(1)+vmult(i,1,1,1)
   10 continue
      epts = glsum(work,1) + .1
      netot=epts
      if (nid.eq.0) write(6,*) ' dofs:',nvtot,nptot,netot
      return
      end
c-----------------------------------------------------------------------
      subroutine add2col2(a,b,c,n)
      real a(1),b(1),c(1)
c
      do i=1,n
         a(i) = a(i) + b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine vol_flow
c
c
c     Adust flow volume at end of time step to keep flow rate fixed by
c     adding an appropriate multiple of the linear solution to the Stokes
c     problem arising from a unit forcing in the X-direction.  This assumes
c     that the flow rate in the X-direction is to be fixed (as opposed to Y-
c     or Z-) *and* that the periodic boundary conditions in the X-direction
c     occur at the extreme left and right ends of the mesh.
c
c     pff 6/28/98
c
      include 'SIZE'
      include 'TOTAL'
c
c     Swap the comments on these two lines if you don't want to fix the
c     flow rate for periodic-in-X (or Z) flow problems.
c
      parameter (kx1=lx1,ky1=ly1,kz1=lz1,kx2=lx2,ky2=ly2,kz2=lz2)
c     parameter (kx1=1,ky1=1,kz1=1,kx2=1,ky2=1,kz2=1)
c
      common /cvflow_a/ vxc(kx1,ky1,kz1,lelv)
     $                , vyc(kx1,ky1,kz1,lelv)
     $                , vzc(kx1,ky1,kz1,lelv)
     $                , prc(kx2,ky2,kz2,lelv)
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_i/ icvflow,iavflow
      common /cvflow_c/ chv(3)
      character*1 chv
c
      real bd_vflow,dt_vflow
      save bd_vflow,dt_vflow
      data bd_vflow,dt_vflow /-99.,-99./
c
c     Check list:
c
c     param (55) -- volume flow rate, if nonzero
c     forcing in X? or in Z?
c
c
c
      if (param(55).eq.0.) return
      if (kx1.eq.1) then
         write(6,*) 'ABORT. Recompile vol_flow with kx1=lx1, etc.'
         call exitt
      endif
c
      icvflow   = 1                                    ! Default flow dir. = X
      if (param(54).ne.0) icvflow = abs(param(54))
      iavflow   = 0                                    ! Determine flow rate from
      if (param(54).lt.0) iavflow = 1                  ! mean velocity
      flow_rate = param(55)
c
      chv(1) = 'X'
      chv(2) = 'Y'
      chv(3) = 'Z'
c
c     If either dt or the backwards difference coefficient change,
c     then recompute base flow solution corresponding to unit forcing:
c
      if (dt.ne.dt_vflow.or.bd(1).ne.bd_vflow.or.param(30).ne.0)
     $   call compute_vol_soln(vxc,vyc,vzc,prc)
      dt_vflow = dt
      bd_vflow = bd(1)
c
      ntot1 = nx1*ny1*nz1*nelv
      ntot2 = nx2*ny2*nz2*nelv
      if (icvflow.eq.1) current_flow=glsc2(vx,bm1,ntot1)/domain_length  ! for X
      if (icvflow.eq.2) current_flow=glsc2(vy,bm1,ntot1)/domain_length  ! for Y
      if (icvflow.eq.3) current_flow=glsc2(vz,bm1,ntot1)/domain_length  ! for Z
c
      if (iavflow.eq.1) then
         xsec = volvm1 / domain_length
         flow_rate = param(55)*xsec
      endif
c
      delta_flow = flow_rate-current_flow
c
c     Note, this scale factor corresponds to FFX, provided FFX has
c     not also been specified in userf.   If ffx is also specified
c     in userf then the true FFX is given by ffx_userf + scale.
c
      scale = delta_flow/base_flow
      scale_vf(icvflow) = scale
      if (nid.eq.0) write(6,1) istep
     $   ,time,scale,delta_flow,current_flow,flow_rate,chv(icvflow)
    1    format(i8,e14.7,1p4e13.5,' volflow',1x,a1)
c
      call add2s2(vx,vxc,scale,ntot1)
      call add2s2(vy,vyc,scale,ntot1)
      call add2s2(vz,vzc,scale,ntot1)
      call add2s2(pr,prc,scale,ntot2)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine compute_vol_soln(vxc,vyc,vzc,prc)
c
c     Compute the solution to the time-dependent Stokes problem
c     with unit forcing, and find associated flow rate.
c
c     pff 2/28/98
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
c
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_i/ icvflow,iavflow
      common /cvflow_c/ chv(3)
      character*1 chv
c
      integer icalld
      save    icalld
      data    icalld/0/
c
c
      ntot1 = nx1*ny1*nz1*nelv
      if (icalld.eq.0) then
         icalld=icalld+1
         xlmin = glmin(xm1,ntot1)
         xlmax = glmax(xm1,ntot1)
         ylmin = glmin(ym1,ntot1)          !  for Y!
         ylmax = glmax(ym1,ntot1)
         zlmin = glmin(zm1,ntot1)          !  for Z!
         zlmax = glmax(zm1,ntot1)
c
         if (icvflow.eq.1) domain_length = xlmax - xlmin
         if (icvflow.eq.2) domain_length = ylmax - ylmin
         if (icvflow.eq.3) domain_length = zlmax - zlmin
c
      endif
c
      if (ifsplit) then
         call plan2_vol(vxc,vyc,vzc,prc)
      else
         call plan3_vol(vxc,vyc,vzc,prc)
      endif
c
c     Compute base flow rate
c 
      if (icvflow.eq.1) base_flow = glsc2(vxc,bm1,ntot1)/domain_length
      if (icvflow.eq.2) base_flow = glsc2(vyc,bm1,ntot1)/domain_length
      if (icvflow.eq.3) base_flow = glsc2(vzc,bm1,ntot1)/domain_length
c
      if (nid.eq.0) write(6,1) 
     $   istep,base_flow,domain_length,flow_rate,chv(icvflow)
    1    format(i9,1p3e13.5,' basflow',1x,a1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine plan2_vol(vxc,vyc,vzc,prc)
c
c     Compute pressure and velocity using fractional step method.
c     (classical splitting scheme).
c
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
C
      COMMON /SCRNS/ RESV1 (LX1,LY1,LZ1,LELV)
     $ ,             RESV2 (LX1,LY1,LZ1,LELV)
     $ ,             RESV3 (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      COMMON /SCRVH/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)
c
      common /cvflow_i/ icvflow,iavflow
C
C
C     Compute pressure 
C
      ntot1  = nx1*ny1*nz1*nelv
c
      if (icvflow.eq.1) then
         call cdtp     (respr,v1mask,rxm2,sxm2,txm2,1)
      elseif (icvflow.eq.2) then
         call cdtp     (respr,v2mask,rxm2,sxm2,txm2,1)
      else
         call cdtp     (respr,v3mask,rxm2,sxm2,txm2,1)
      endif
c
      call ortho    (respr)
c
      call ctolspl  (tolspl,respr)
      call rone     (h1,ntot1)
      call rzero    (h2,ntot1)
c
      call hmholtz  ('pres',prc,respr,h1,h2,pmask,vmult,
     $                             imesh,tolspl,nmxh,1)
      call zaver1   (prc)
C
C     Compute velocity
C
      call opgrad   (resv1,resv2,resv3,prc)
      call opchsgn  (resv1,resv2,resv3)
      call add2col2 (resv1,bm1,v1mask,ntot1)
c
      intype = -1
      call sethlm   (h1,h2,intype)
      call ophinv   (vxc,vyc,vzc,resv1,resv2,resv3,h1,h2,tolhv,nmxh)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine plan3_vol(vxc,vyc,vzc,prc)
c
c     Compute pressure and velocity using fractional step method.
c     (PLAN3).
c
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
C
      COMMON /SCRNS/ rw1   (LX1,LY1,LZ1,LELV)
     $ ,             rw2   (LX1,LY1,LZ1,LELV)
     $ ,             rw3   (LX1,LY1,LZ1,LELV)
     $ ,             dv1   (LX1,LY1,LZ1,LELV)
     $ ,             dv2   (LX1,LY1,LZ1,LELV)
     $ ,             dv3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      COMMON /SCRVH/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)
      COMMON /SCRHI/ H2INV (LX1,LY1,LZ1,LELV)
      common /cvflow_i/ icvflow,iavflow
c
c
c     Compute velocity, 1st part 
c
      ntot1  = nx1*ny1*nz1*nelv
      ntot2  = nx2*ny2*nz2*nelv
      ifield = 1
c
      if (icvflow.eq.1) then
         call copy     (rw1,bm1,ntot1)
         call rzero    (rw2,ntot1)
         call rzero    (rw3,ntot1)
      elseif (icvflow.eq.2) then
         call rzero    (rw1,ntot1)
         call copy     (rw2,bm1,ntot1)
         call rzero    (rw3,ntot1)
      else
         call rzero    (rw1,ntot1)        ! Z-flow!
         call rzero    (rw2,ntot1)        ! Z-flow!
         call copy     (rw3,bm1,ntot1)    ! Z-flow!
      endif
      intype = -1
      call sethlm   (h1,h2,intype)
      call ophinv   (vxc,vyc,vzc,rw1,rw2,rw3,h1,h2,tolhv,nmxh)
      call ssnormd  (vxc,vyc,vzc)
c
c     Compute pressure  (from "incompr")
c
      intype = 1
      dtinv  = 1./dt
c
      call rzero   (h1,ntot1)
      call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
      call cmult   (h2,dtinv,ntot1)
      call invers2 (h2inv,h2,ntot1)
      call opdiv   (respr,vxc,vyc,vzc)
      call chsign  (respr,ntot2)
      call ortho   (respr)
c
c
c     Set istep=0 so that h1/h2 will be re-initialized in eprec
      i_tmp = istep
      istep = 0
      call esolver (respr,h1,h2,h2inv,intype)
      istep = i_tmp
c
      call opgradt (rw1,rw2,rw3,respr)
      call opbinv  (dv1,dv2,dv3,rw1,rw2,rw3,h2inv)
      call opadd2  (vxc,vyc,vzc,dv1,dv2,dv3)
c
      call cmult2  (prc,respr,bd(1),ntot2)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine a_dmp
c
      include 'SIZE'
      include 'TOTAL'
      COMMON /SCRNS/ w(LX1,LY1,LZ1,LELT)
      COMMON /SCRUZ/ v (LX1,LY1,LZ1,LELT)
     $             , h1(LX1,LY1,LZ1,LELT)
     $             , h2(LX1,LY1,LZ1,LELT)
c
      ntot = nx1*ny1*nz1*nelv
      call rone (h1,ntot)
      call rzero(h2,ntot)
      do i=1,ntot
         call rzero(v,ntot)
         v(i,1,1,1) = 1.
         call axhelm (w,v,h1,h2,1,1)
         call outrio (w,ntot,55)
      enddo
c     write(6,*) 'quit in a_dmp'
c     call exitt
      return
      end
c-----------------------------------------------------------------------
      subroutine outrio (v,n,io)
c
      real v(1)
c
      write(6,*) 'outrio:',n,io,v(1)
      write(io,6) (v(k),k=1,n)
    6 format(1pe19.11)
c
c     nr = min(12,n)
c     write(io,6) (v(k),k=1,nr)
c   6 format(1p12e11.3)
      return
      end
c-----------------------------------------------------------------------
      subroutine reset_prop
C------------------------------------------------------------------------
C
C     Set variable property arrays
C
C------------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
C
C     Caution: 2nd and 3rd strainrate invariants residing in scratch
C              common /SCREV/ are used in STNRINV and NEKASGN
C
      COMMON /SCREV/ SII (LX1,LY1,LZ1,LELT)
     $             , SIII(LX1,LY1,LZ1,LELT)
      COMMON /SCRUZ/ TA(LX1,LY1,LZ1,LELT)
C
      real    rstart
      save    rstart
      data    rstart  /1/
c
      rfinal   = 1./param(2) ! Target Re
c
      ntot  = nx1*ny1*nz1*nelv
      iramp = 200
      istpp = istep
c     istpp = istep+2033+1250
      if (istpp.ge.iramp) then
         vfinal=1./rfinal
         call cfill(vdiff,vfinal,ntot)
      else
         one = 1.
         pi2 = 2.*atan(one)
         sarg  = (pi2*istpp)/iramp
         sarg  = sin(sarg)
         rnew = rstart + (rfinal-rstart)*sarg
         vnew = 1./rnew
         call cfill(vdiff,vnew,ntot)
         if (nid.eq.0) write(6,*) istep,' New Re:',rnew,sarg,istpp
      endif
      return
      end
C-----------------------------------------------------------------------
