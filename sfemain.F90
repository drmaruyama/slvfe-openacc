! -*- F90 -*-
! ERmod - Eneregy Representation Module
! Copyright (C) 2000-2024 Nobuyuki Matubayasi
! Copyright (C) 2010-2024 Shun Sakuraba
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

module sysvars
   implicit none

   character(len=5) :: clcond = 'merge'
   character(len=3) :: uvread = 'yes',    slfslt = 'yes',   ljlrc = 'not'
   character(len=3) :: infchk = 'not',    meshread = 'not', cumuint = 'not'
   character(len=3) :: write_mesherror = 'cnd'
   character(len=4) :: zerosft = 'orig',  wgtfnform = 'harm'
   character(len=3) :: refmerge = 'yes',  extsln = 'lin'
   character(len=3) :: wgtf2smpl = 'yes', slncor = 'not'
   character(len=3) :: normalize = 'yes', showdst= 'not'
   character(len=3) :: wrtzrsft = 'not',  readwgtfl = 'yes'

   integer :: numprm = 0                           ! initialized to 0
   integer :: numprm_def_inf_yes = 11   ! default numprm at infchk = 'yes'
   integer :: numprm_def_inf_not = 5    ! default numprm at infchk = 'not'

   integer :: numsln = 0, numref = 0, numdiv = 0   ! initialized to 0
   integer :: maxsln, maxref, numrun, prmmax
   integer :: numslv, ermax

   real(kind=8) :: inptemp = 300.0        ! temperature in Kelvin, initialized
   real(kind=8) :: temp, kT, slfeng
   real(kind=8) :: avevolume = 0.0        ! average volume of system, initialized

   integer :: pickgr = 3
   integer :: msemin = 1, msemax = 5
   real(kind=8) :: mesherr = 0.1          ! allowed mesh error in kcal/mol

   integer :: extthres_soln = 1, extthres_refs = 1
   integer :: minthres_soln = 0, minthres_refs = 0
   real(kind=8), parameter :: zero = 0.0
!   real(kind=8) :: error = 1.0e-8, tiny = 1.0e-8
   real(kind=8) :: error = 1.0e-8, tiny = 1.0e-5
   integer :: ermax_limit = 15000
   integer :: large = 500000, itrmax = 100

   character(len=1024) :: solndirec = 'soln'
   character(len=1024) :: refsdirec = 'refs'
   character(len=1024) :: wgtslnfl  = 'weight_soln'
   character(len=1024) :: wgtreffl  = 'weight_refs'
   character(len=1024) :: slndnspf  = 'engsln'
   character(len=1024) :: slncorpf  = 'corsln4'
   character(len=1024) :: refdnspf  = 'engref'
   character(len=1024) :: refcorpf  = 'corref4'
   character(len=1024) :: aveuvfile = 'aveuv.tt'
   character(len=1024) :: engmeshfile = 'EngMesh'
   character(len=1024) :: cumuintfl = 'cumsfe'
   character(len=10), parameter :: numbers='0123456789'

   real(kind=8), dimension(:),     allocatable :: nummol
   integer, dimension(:),  allocatable :: rduvmax, rduvcore
   real(kind=8), dimension(:),     allocatable :: rdcrd, rddst, rddns
   real(kind=8), dimension(:,:),   allocatable :: rdslc, rdcor
   integer, dimension(:),  allocatable :: rdspec
   real(kind=8), dimension(:,:,:), allocatable :: chmpt
   real(kind=8), dimension(:),     allocatable :: aveuv
   real(kind=8), dimension(:,:),   allocatable :: uvene, blockuv
   integer, dimension(:),  allocatable :: svgrp, svinf
   real(kind=8), dimension(:),     allocatable :: wgtsln, wgtref

   logical :: force_calculation = .false., strict_ewald_parameters = .false.

   namelist /fevars/ clcond, numprm, numsln, numref, numdiv, &
      uvread, slfslt, infchk, meshread, zerosft, wgtfnform, &
      refmerge, extsln, extthres_soln, extthres_refs, &
      minthres_soln, minthres_refs, &
      wgtf2smpl, slncor, normalize, showdst, wrtzrsft, readwgtfl, &
      inptemp, pickgr, write_mesherror, msemin, msemax, mesherr, &
      ljlrc, avevolume, &
      solndirec, refsdirec, wgtslnfl, wgtreffl, &
      slndnspf, slncorpf, refdnspf, refcorpf, &
      aveuvfile, engmeshfile, cumuint, cumuintfl, &
      ermax_limit, large, itrmax, error, tiny, &
      force_calculation, strict_ewald_parameters

contains

   subroutine init_sysvars
      implicit none
      character(len=*), parameter :: parmfname = 'parameters_fe'
      character(len=3) :: file_suf
      character(len=1024) :: opnfile
      integer, parameter :: iounit = 191, sufmax = 99
      integer :: ioerr, count_suf, i, j, srcnt, count_soln, count_refs
      logical :: file_exist

      open(unit = iounit, file = parmfname, action = 'read', status = 'old', iostat = ioerr)
      if(ioerr /= 0) goto 99
      read(iounit, nml = fevars)
      close(iounit)
99    continue

      if(clcond == 'merge') then
         do srcnt = 1, 2
            do count_suf = 1, sufmax
               i = count_suf / 10
               j = mod(count_suf, 10)
               file_suf = '.' // numbers(i+1:i+1) // numbers(j+1:j+1)
               select case(srcnt)
                case(1)
                  opnfile = trim(solndirec) // '/' // trim(slndnspf) // file_suf
                case(2)
                  opnfile = trim(refsdirec) // '/' // trim(refdnspf) // file_suf
               end select
               inquire(file = opnfile, exist = file_exist)
               if( file_exist ) then
                  if(srcnt == 1) count_soln = count_suf
                  if(srcnt == 2) count_refs = count_suf
               else
                  exit
               endif
            enddo
         enddo

         open(unit = iounit, file = parmfname, action = 'read', status = 'old', iostat = ioerr)
         if(ioerr /= 0) goto 91
         read(iounit, nml = fevars)
         close(iounit)
91       continue

         if((numsln <= 0) .or. (numsln > count_soln)) numsln = count_soln
         if((numref <= 0) .or. (numref > count_refs)) numref = count_refs

         if((numdiv <= 0) .or. (numdiv >= numsln)) numdiv = numsln
         if(mod(numsln, numdiv) /= 0) then
            do i = numdiv + 1, numsln      ! find the larger and closest divisor
               if(mod(numsln, i) == 0) exit
            enddo
            numdiv = i
         endif
         if(refmerge == 'not') then        ! see subroutine datread for refmerge
            if(numdiv > numref) stop " With refmerge = 'not', numdiv needs to be not larger than numref"
            if(mod(numref, numdiv) /= 0) then
               write(6, "(A,i2,A,i2,A)") " Note: only ", numdiv * (numref / numdiv), &
               &" files out of ", numref, " engref and corref files prepared"
               numref = numdiv * (numref / numdiv)
            endif
         endif
         call check_params() ! check consistency in parameters_er
      endif

      if(numprm <= 0) then                 ! default setting
         if(infchk == 'yes') then
            numprm = numprm_def_inf_yes    ! default numprm at infchk = 'yes'
         else
            numprm = numprm_def_inf_not    ! default numprm at infchk = 'not'
         endif
      endif

      if(pickgr < msemin) stop " Incorrect setting: pickgr < msemin not allowed"
      if(pickgr > msemax) stop " Incorrect setting: pickgr > msemax not allowed"
      if(pickgr > numprm) stop " Incorrect setting: pickgr > numprm not allowed"

   end subroutine init_sysvars

   ! Check parameter consistency (issue #17)
   subroutine check_params
      use engmain, engmain_force_calculation=>force_calculation
      implicit none
      integer :: boxshp_s, estype_s, insposition_s, insstructure_s
      real(kind=8) :: lwreg_s, upreg_s, lwstr_s, upstr_s, inptemp_s
      integer :: ljformat_s, ljswitch_s, cmbrule_s
      real(kind=8) :: lwljcut_s, upljcut_s, elecut_s, screen_s, ewtoler_s
      integer :: splodr_s, cltype_s, ms1max_s, ms2max_s, ms3max_s
      logical :: inconsistent
      
      ! load both parameters_er and compare varibles
      lwreg = -1; upreg = -1; lwstr = -1; upstr = -1
      call init_params(solndirec)
      boxshp_s = boxshp; estype_s = estype; insposition_s = insposition; insstructure_s = insstructure
      lwreg_s = lwreg; upreg_s = upreg; lwstr_s = lwstr; upstr_s = upstr; inptemp_s = inptemp
      ljformat_s = ljformat; ljswitch_s = ljswitch; lwljcut_s = lwljcut; upljcut_s = upljcut; cmbrule_s = cmbrule
      elecut_s = elecut; screen_s = screen; ewtoler_s = ewtoler; splodr_s = splodr
      cltype_s = cltype; ms1max_s = ms1max; ms2max_s = ms2max; ms3max_s = ms3max

      call init_params(refsdirec)

      inconsistent = .false.

      if((boxshp /= boxshp_s) .or. (estype /= estype_s)) then
         print *, "box shape / system type inconsistent"
         inconsistent = .true.
      end if
      if((cltype /= cltype_s) .or. &
         (lwljcut /= lwljcut_s) .or. &
         (upljcut /= upljcut_s) .or. &
         (ljswitch /= ljswitch_s) .or. &
         (ljformat /= ljformat_s) .or. &
         (cmbrule /= cmbrule_s) .or. &
         (splodr /= splodr_s) .or. & ! splodr is technically ok, but different splodrs look *very* wrong...
         (elecut /= elecut_s)) then
         print *, "Nonbond calculation conditions (electrostatic, LJ) inconsistent"
         inconsistent = .true.
      end if
      if(((insposition /= insposition_s) .and. &
          ((insposition /= INSPOS_RANDOM) .and. (insposition /= INSPOS_NOCHANGE) .and. (insposition /= INSPOS_GAUSS))) .or. &
         ((insstructure /= insstructure_s) .and. (insstructure /= INSSTR_NOREJECT)) .or. &
         (lwreg /= lwreg_s) .or. &
         (upreg /= upreg_s) .or. &
         (lwstr /= lwstr_s) .or. &
         (upstr /= upstr_s)) then
         print *, "Insertion conditions inconsistent between solution and reference."
         print *, " Test particle insertion parameters are also used in the solution system to reject" // &
            " snapshots in the solution system. We thus need to use same parameters for:" // &
            "insposition, insstructure, lwreg, upreg, lwstr, and upstr"
         inconsistent = .true.
      end if

      if(inconsistent) then
         if(force_calculation) then
            print *, "Proceeding the calculation because force_calculation is set"
         else
            print *, "==>Aborting because solution and reference systems do not match"
            print *, "This typically mean you need to rerun the simulation with matching conditions"
            print *, 'If you dare to proceed, set "force_calculation=.true." in parameters_fe'
            stop "Stopping the calculation due to inconsistency"
         end if
      end if

      if((screen /= screen_s) .or. &
         (ms1max /= ms1max_s) .or. &
         (ms2max /= ms2max_s) .or. &
         (ms3max /= ms3max_s) .or. &
         (ewtoler /= ewtoler_s)) then
         print *, "Some of Ewald parameters are inconsistent"
         print *, "beta =", screen_s, "(soln) /", screen, "(refs)"
         print *, "grid =", ms1max_s, ",", ms2max_s, ",", ms3max_s, "(soln) /", &
            ms1max, ",", ms2max, ",", ms3max, "(refs)"
         print *, "Ewald tolerance =", ewtoler_s, "(soln)", ewtoler, "(refs)"
         print *, "This is mostly harmless and slvfe continues;"
         print *, " but keep in mind that the F.E. estimation may become unstable due to this difference"
         print *, " (if you cannot prevent this inconsistency, keep ewald tolerance low in the SIMULATION software)"
         if(strict_ewald_parameters) then
            stop "Aborting due to strict condition"
         end if
      end if
   end subroutine check_params
end module sysvars
