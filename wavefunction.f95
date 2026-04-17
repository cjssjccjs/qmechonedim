!***********************************************************************************************************************************
! wavefunction.f95
! currently set up for the ground state of molecular oxygen with the Morse potential
! last modified 1/26/2026 by cjs
!
! This source code computes the evolution of a one-dimensional quantum-mechanical wavefunction and can simultaneously display it to
! create a "movie" if told to. It can also calculate energy eigenvalues and position-space energy eigenfunctions for bound states.
! It's written in Fortran and generates graphics using PlPlot subroutines. When referencing it please cite the journal article it
! goes with: C. J. Sweeney, "A few old quantum mechanics problems... revisited... yet again," Journal of Numerical Simulations in
! Physics and Mathematics ??, ??-?? (2026). Use this article as guidance to understand how the program works and to troubleshoot 
! problems. The recommended compiler for the source code is GNU Fortran, which can be freely downloaded from the internet along 
! with PlPlot.
!
! Suggested usage instructions:
!    1. Choose a potential and initial wavefunction.
!    2. Set program parameters following the restrictions on them enumerated in the journal article just cited.
!    3. Decide on the graphics format:
!       A. Display evolution?
!          Yes--> qshowwvfn = true
!                 I.   Show "movie" instead of selected still frames?
!                      Yes--> qshowmovie = true
!                      No-->  qshowmovie = false
!                 II.  Show probability density instead of wavefunction?
!                      Yes--> qprob = true
!                      No-->  qprob = false
!                 III. Show initial configuration?
!                      Yes--> qshowinit = true
!                      No-->  qshowinit = false
!          No-->  qshowwvfn = false
!       B. Determine bound-state energy eigenvalues?
!          Yes--> qEIT = true
!                 Show them?
!                      Yes--> qshowEIT = true
!                      No-->  qshowEIT = false
!          No-->  qEIT = false
!       C. Determine bound-state position-space energy eigenfunctions?
!          Yes--> qEIF = true
!                 Show them?
!                      Yes--> qshowEIF = true
!                      No-->  qshowEIF = false
!          No-->  qEIF = true
!       D. Save eigenvalues, eigenfunctions and the latters' weights?
!          Yes--> qsave = true
!          No-->  qsave = false
!    4. Compile and run the program.
!
! Miscellany:
!    1. The usage instructions just described are hierarchical; eigenfunctions can't be calculated without first having eigen-
!       values and eigenvalues can't be calculated without first propagation the system's vavefunction.
!    2. This source code has been tested and to the best of our knowledge has no flaws. Nonetheless, we make no warranties, implied
!       or expressesd, that the source code herein is free of error. It should not be used for any purpose in which incorrect
!       results that it might produce could result in personal injury, loss of life, property damage or property loss. All use of
!       the source code, including modified versions of it, is at the user's own risk. We disclaim all liability for direct or
!       consequential damage resulting use of this source code or anything derived from it or its modifications.
!    3. The source code employs the pseudo-spectral method of Fleck, et al. [Cf., Feit M D, Fleck J A Jr and Steiger A 1982 "Solu-
!       tion of the Schroedinger Equation by a Spectral Method," Journal of Computational Physics 47 pp 412-33.].
!    4. The units used are "natural," where hbar = c = e = 1.
!    5. Most of the documentation for this source code appears before the program itself; you are reading the bulk of it right
!       now.
!    6. The source code can easily be set to any desired degree of numerical precision.
!    7. Some of the constructions are obsolete, yet perhaps clearer than the equivalent modern constructions. The obsolete ones
!       have been verified to run under the GNU compiler, but not others.
!    8. A color table can be user defined, or colors can be set according to this standard PlPlot table:
!       0                          black (default background)
!       1                          red (default foreground)
!       2                          yellow
!       3                          green
!       4                          aquamarine
!       5                          pink
!       6                          wheat
!       7                          gray
!       8                          brown
!       9                          blue
!       10                         blue violet
!       11                         cyan
!       12                         turquoise
!       13                         magenta
!       14                         salmon
!       15                         white
!    9. Upper and lower case variables have different meanings; this can easily be changed if too problematic.
!    10. Variables are case-sensitive; capital and lower-case letters have different meanings. Make sure the compiler can accomodate
!        this---or rename some variables.
!    11. Some variables have different meaning in different parts of the program.
!    12. The parameters and variables in this program and what they represent are:
!       Parameter or Variable      Meaning
!       a                          slope for graphics scaling; Morse exponential factor
!       autocorr                   array containing temporal autocorrelation of the wavefunction
!       b                          intercept for graphics scaling
!       D                          depth of Morse potential
!       dE                         energy increment
!       delta_prime                bound-state energy eigenvalue offset
!       delx                       position-space "spread" of wavefunction
!       dk                         wavenumber-space increment
!       ds                         generic increment
!       dt                         time increment between displayed steps of the wavefunction
!       dtpropagation              time increment between hidden steps of the wavefunction
!       dx                         position-space increment
!       E                          array of corrected bound-state energy eigenvalues
!       E_n                        analytic energy bound-state energy eigenvalue
!       energy                     array containing uncorrected bound-state energy eigenvalues
!       i                          positive square root of negative one
!       j                          generic counter
!       j1                         generic counter
!       j2                         generic counter
!       j3                         generic counter
!       j4                         generic counter
!       j5                         generic counter
!       j6                         generic counter
!       j7                         generic counter
!       jcolor                     color counter
!       jframe                     frame counter
!       jhidden                    counter for hidden propagations
!       jjcolor                    color array
!       jlocation                  jth bound-state energy eigenvalue counter
!       jminus                     -1
!       jpeak                      counter for locating bound-state energy eigenvalues
!       jplus                      1
!       jsign                      generic +/-1
!       jstep                      time step counter
!       k                          array of values of the wavenumbers treated at---not arranged in intuitive order!!!
!       k0                         initial average wavenumber of wavefunction
!       location                   array containing locations of bound-state energy eigenvalues
!       m                          power to which 2 will be raised to spell out size of array containing the wavefunction
!       mEIT                       array of eigenfunctions to be displayed
!       msteps                     power to which 2 will be raised to give the number displayed time steps of the wavefunction
!       mu                         mass (often reduced) of the system described by the wavefunction
!       n                          size m**2 of array containing the wavefunction
!       nEIT                       number of bound-state energy eigenvalues to be computed
!       nh                         number of horizontal frames high display is
!       nhalf                      n/2 for FFT
!       nhidden                    number of time steps between successive displayed ones
!       noise                      noise threshold for eigenvalue determination
!       nshowEIT                   number of eigenfunctions to display
!       nsteps                     number msteps**2 of displayed time steps of the wavefunction
!       nw                         number of horizontal frames wide display is
!       pi                         pi to way more decimal places than anyone would care about...
!       potential                  potential-energy function
!       psi                        array containing the wavefunction
!       psi_sum                    array containing wvefunction for eigenfunction computation
!       psi0                       array holding initial wavefunction
!       psin                       array holding nth bound-state eigenfunction
!       ptntl                      array containing calculated values of potential-energy function
!       qEIF                       .true. if eigenfunctions are to be calculated
!       qEIT                       .true. if eigenvalues are to be calculated
!       qnormalize                 .true. to normalize FFT
!       qprob                      .true. to display probability density instead of wavefunction
!       qshowEIF                   .true. to display eigenfunctions
!       qshowEIT                   .true. to display eigenvalues
!       qshowinit                  .true. to continuously display initial system configuration
!       qshowmovie                 .true. to display movie
!       qshowwvfn                  .true. to display wavefunction in any form, as movie or as frames from it
!       r                          internuclear separation
!       R                          variable for calculating bound-state energy eigenvalue correction
!       r_e                        bond length
!       renorm                     normalization prefactor for initial wavefunction; same for FFT
!       rr                         another variable for calculating bound-state energy eigenvalue correction
!       scale                      fraction of screen used for display
!       Smax                       maximum value of spectrum
!       spectrum                   array containg modulus of temporal autocorrelation of the wavefunction
!       sum                        temporary summation variable
!       t                          time
!       T                          array of values for the exponentiated kinetic-energy operator
!       temp                       intermediate result for FFT
!       time                       array of time steps values when wavefunction or probability density is displayed
!       theta                      phase angle for the FFT
!       u                          array of normalized bound-state energy eigenfunctions; intermediate result for FFT
!       umax                       maximum value of bound-state eigenfunctions
!       Vdouble                    array of values for exponentiated potential-energy operator for hidden propagations
!       Vsingle                    array containing values of exponentiated potential-enrgy operator for displayed propagations
!       vinf                       potential's value at infinity
!       vmax                       maximum potential-energy function value
!       vmin                       minimum potential-energy function value
!       vv                         array of displayed potential-energy function
!       w                          Hann weighting function; intermediate result for FFT
!       x                          array of values of the positions the wavefunction is treated at---is arranged in intuitive order
!       x0                         average position of intitial wavefunction if the latter is Guassian
!       xlo                        lower bound of range treated in position space
!       xrange                     range of position space that the wavefunction is treated on
!       ymin                       minimum value of y for bottom of graphics window
!       zz                         array of displayed probability density of wavefunction
!       zzim                       scaled imaginary part of psi
!       zzmax                      maximum value of wavefunction or probability density
!       zzre                       scaled real part of psi
!       zztop                      a rock band with an unusually skilled guitarist; also a scaling factor for the displayed wave-
!                                  function or probability density
!       zzz                        scaled initial probability density
!       zzzim                      imaginary part of scaled initial probability density
!       zzzre                      real part of scaled initial probability density


!***********************************************************************************************************************************
! This is the driver that calls the principal subroutines; the whole program starts and ends here.
       program                                    wavefunction
       implicit                                   none
       integer(kind = 8), parameter ::            m = 10, msteps = 18, n = 2**m, nEIT = 50, nsteps = 2**msteps
       real(kind = 8), dimension(nEIT) ::         E
       complex(kind = 8), dimension(nsteps) ::    autocorr
       complex(kind = 8), dimension(n) ::         psi

       call initialize(psi)
       call propagate(autocorr, psi)
       call EIT(autocorr, E)
       call EIF(E)
       stop

       end program                                wavefunction


!***********************************************************************************************************************************
! Here many of the parameters needed for computations are set up. This includes putting lots of static values in common blocks for
! future use with other subroutines.
       subroutine                                 initialize(psi)
       implicit                                   none
       integer(kind = 8) ::                       j
       integer(kind = 8), parameter ::            m = 10, n = 2**m, nhidden = 9
       real(kind = 8) ::                          dk, dtpropagation, dx, potential
       real(kind = 8), parameter ::               pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164063d0
       real(kind = 8), parameter ::               dt = 1.0d0, k0 = 3.0d0, mu = 5.0d-1*2.6566963d-26/9.1093837d-31, &
                                                  xrange = 1.2d1, x0 = 1.207d-10/5.29167d-11 + 3.0d-1, xlo = 0.0d0, &
                                                  delx = 1.0d-2*xrange
       real(kind = 8), dimension(n) ::            k, ptntl, x
       complex(kind = 8), parameter ::            i = (0.0d0, 1.0d0)
       complex(kind = 8), dimension(n) ::         psi, T, Vdouble, Vsingle
       common/evolution/                          T, Vdouble, Vsingle
       common/position/                           x
       common/ptl/                                ptntl
       common/xincr/                              dx

       dk            = 2.0d0*pi/xrange
       dtpropagation = dt/dble(nhidden + 1)
       dx            = xrange/dble(n)
       do j = 1, n/2
          k(j)         =  dble(j - 1)*dk
          k(n + 1 - j) = -dble(j)*dk
       end do
       do j = 1, n
          x(j) = dble(j - 1)*dx + xlo
       end do
       do j = 1, n
          ptntl(j) = potential(x(j))
       end do
       psi     = (1.0d0/(dsqrt(dsqrt(2.0d0*pi)*delx))) &
               * dexp(-2.5d-1*(x - x0)*(x - x0)/(delx*delx))*cdexp(i*k0*x)
       T       = cdexp(-5.0d-1*i*dtpropagation*k*k/mu)
       Vdouble = cdexp(-i*dtpropagation*ptntl)
       Vsingle = cdexp(-5.0d-1*i*dtpropagation*ptntl)
       return

       end subroutine                             initialize


!***********************************************************************************************************************************
! The potential energy function is set up here. Currently it's the Morse function for the ground state of molecular oxygen.
       function                                   potential(r)
       implicit                                   none
       real(kind = 8) ::                          potential, r
       real(kind = 8), parameter ::               a = 2.78d0*5.29167d-11/1.0d-10, D = 5.211d0*1.6021766d-19/4.35944d-18, &
                                                  r_e = 1.207d-10/5.29167d-11

       potential = D*dexp(-2.0d0*a*(r - r_e)) - 2.0d0*D*dexp(-a*(r - r_e)) + D
       return

       end function                               potential


!***********************************************************************************************************************************
! This subroutine propagates the wavefunction. If needed, then it also saves results required for determination of bound-state
! energy eigenvalues and eigenfunctions.
       subroutine                                 propagate(autocorr, psi)
       use                                        plplot
       implicit                                   none
       logical, parameter ::                      qEIT = .true., qnormalize = .true., qshowwvfn = .true.
       integer(kind = 8) ::                       j, jframe, jhidden, jstep
       integer(kind = 8), parameter ::            jminus = -1, jplus = 1, m = 10, msteps = 18, n = 2**m, nhidden = 9, &
                                                  nsteps = 2**msteps
       real(kind = 8) ::                          dx
       complex(kind = 8), dimension(nsteps) ::    autocorr
       complex(kind = 8), dimension(n) ::         psi, psi0, T, Vdouble, Vsingle
       common/evolution/                          T, Vdouble, Vsingle
       common/xincr/                              dx
       common/init/                               psi0

       if (qshowwvfn .eqv. .true.) then
          jframe = 1
          call start_graphics(psi)
       end if
       if (qEIT .eqv. .true.) then
          psi0 = psi
          open(unit = 1, file = 'wavefunction.dat', status = 'new', form = 'unformatted')
       end if
       do jstep = 1, nsteps
          if (qshowwvfn .eqv. .true.) then
             call plot_frame(jframe, jstep, psi)
          else
             write(0, 10) jstep, nsteps, nhidden
          end if
          if (qEIT .eqv. .true.) then
             write(1) psi
             autocorr(jstep) = (0.0d0, 0.0d0)
             do j = 1, n
                autocorr(jstep) = autocorr(jstep) + dconjg(psi0(j))*psi(j)
             end do
             autocorr(jstep) = autocorr(jstep)*dx
          end if
          psi = Vsingle*psi
          call FFT(dx, jminus, m, qnormalize, psi)
          psi = T*psi
          call FFT(dx, jplus, m, qnormalize, psi)
          do jhidden = 1, nhidden
             psi = Vdouble*psi
             call FFT(dx, jminus, m, qnormalize, psi)
             psi = T*psi
             call FFT(dx, jplus, m, qnormalize, psi)
          end do
          psi = Vsingle*psi
       end do
       if (qshowwvfn .eqv. .true.) then
          call plend()
       end if
       if (qEIT .eqv. .true.) then
          close(1)
       end if
       write(0, 20)
       write(0, 30)
       return

10     format(1x, 'At step', 1x, i6, 1x, 'of', 1x, i6, 1x, 'with', 1x, i3, 1x, 'hidden steps.')
20     format(1x, 'All done with wavefunction propagation...')
30     format(1x)
       end subroutine                             propagate


!***********************************************************************************************************************************
! This subroutine determines the energy eigenvalues of any bound states. When shifting is calculated there's a "mystery" plus sign
! that according to Fleck et al. should be a minus sign. I can't figure out the sign error, except that the plus sign always gives
! the right answer and the minus sign the wrong one...
       subroutine                                 EIT(autocorr, E)
       implicit                                   none
       logical, parameter ::                      qEIT = .true., qnormalize = .true., qshowEIT = .true.
       integer(kind = 8), parameter ::            jplus = 1, msteps = 18, nEIT = 50, nhidden = 9, nsteps = 2**msteps
       integer(kind = 8) ::                       j1, j2, jlocation, jpeak
       integer(kind = 8), dimension(nEIT) ::      location
       real(kind = 8) ::                          dE, delta_prime, E_n, Smax, R, rr, vmax
       real(kind = 8), parameter ::               pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164063d0
       real(kind = 8), parameter ::               D = 5.211d0*1.6021766d-19/4.35944d-18, dt = 1.0d0, noise = 1.0d-15
       real(kind = 8), dimension(nEIT) ::         E
       real(kind = 8), dimension(nsteps) ::       energy, spectrum
       complex(kind = 8), dimension(nsteps) ::    autocorr
       common/EITgraphics/                        energy, Smax, vmax

       if (qEIT .eqv. .false.) then
          return
       end if
       location = 0
       E        = 0.0d0
       dE       = 2.0d0*pi/(dble(nsteps)*dt)
       do j1 = 1, nsteps
          energy(j1)   = dE*dble(j1 - 1)
          autocorr(j1) = autocorr(j1)*(1.0d0 - dcos(2.0d0*pi*dble(j1 - 1)/dble(nsteps)))
       end do
       call FFT(dE, jplus, msteps, qnormalize, autocorr)
       vmax = D
       spectrum = cdabs(autocorr)
       Smax     = 0.0d0
       do j1 = 1, nsteps
          if (spectrum(j1) > Smax) then
             Smax = spectrum(j1)
          end if
       end do
       call plot_EIT(spectrum)
       Smax = Smax*noise
       jlocation = 1
       do j1 = 4, nsteps - 4
          jpeak = 0
          do j2 = j1 - 4, j1 + 4
             if (spectrum(j1) >= spectrum(j2)) then
                jpeak = jpeak + 1
             end if
          end do
          if((jpeak >= 9) .and. (spectrum(j1)) > Smax) then
             location(jlocation) = j1
             jlocation           = jlocation + 1
             if (jlocation > nEIT) then
                go to 10
             end if
          end if
       end do
10     if (qshowEIT .eqv. .true.) then
          write(0, 20)
       end if
       do j1 = 1, nEIT
          R  = spectrum(location(j1) + 1)/spectrum(location(j1) - 1)
          rr = (1.0d0 + R)/(1.0d0 - R)
          if (R < 1.0d0) then
             delta_prime = 5.0d-1*(-3.0d0*rr + dsqrt(9.0d0*rr*rr - 8.0d0))
          else if (R > 1.0d0) then
             delta_prime = 5.0d-1*(-3.0d0*rr - dsqrt(9.0d0*rr*rr - 8.0d0))
          else
             write(0, 30)
          end if
          E(j1) = energy(location(j1)) + 2.0d0*pi*delta_prime/(dble(nsteps*(nhidden + 1))*dt)
          if (qshowEIT .eqv. .true.) then
             write(0, 40) j1 - 1, E_n(j1 - 1), E(j1), -(D - E_n(j1 - 1)), -(D - E(j1))
          end if
       end do
       write(0, 50)
       write(0, 60)

       return
20     format(6x, 'n', 4x, 'Ecorrect', 10x, 'Ecalculated', 7x, 'Ecorrect', 10x, 'Ecalculated')
30     format(1x, 'There was an error in energy eigenvalue computation.')
40     format(4x, i3, 2x, 1pd16.8, 2x, 1pd16.8, 2x, 1pd16.8, 2x, 1pd16.8)
50     format(1x, 'All done with calculation of bound-state energy eigenvalues...')
60     format(1x)
       end subroutine                             EIT


!***********************************************************************************************************************************
! Position-space representation for bound-state energy eigenfunctions are determined here.
       subroutine                                 EIF(E)
       implicit                                   none
       logical ::                                 qEIF = .true., qEIT = .true., qshowEIF = .true.
       integer(kind = 8) ::                       j1, j2, j3
       integer(kind = 8), parameter ::            m = 10, msteps = 18, n = 2**m, nEIT = 50, nsteps = 2**msteps
       real(kind = 8) ::                          dx, sum, t, w
       real(kind = 8), parameter ::               pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164063d0
       real(kind = 8), parameter ::               dt = 1.0d0
       real(kind = 8), dimension(nEIT) ::         E
       real(kind = 8), dimension(n) ::            x
       real(kind = 8), dimension(n, nEIT) ::      u
       complex(kind = 8), parameter ::            i = (0.0d0, 1.0d0)
       complex(kind = 8), dimension(n) ::         psi, psi_sum
       common/position/                           x
       common/xincr/                              dx

       if (qEIT .eqv. .false.) then
          return
       else if (qEIF .eqv. .false.) then
          return
       end if
       do j1 = 1, nEIT
          write(0, 10) j1 - 1
          psi_sum = (0.0d0, 0.0d0)
          open(unit = 1, file = 'wavefunction.dat', status = 'old', form = 'unformatted')
          do j2 = 1, nsteps
             read(1) psi
             t = dble(j2 - 1)*dt
             w = 1.0d0 - dcos(2.0d0*pi*t/(dble((nsteps - 1))*dt))
             do j3 = 1, n
                psi_sum(j3) = psi_sum(j3) + cdexp(i*E(j1)*t)*w*psi(j3)
             end do
          end do
          close(1)
          do j2 = 1, n
             u(j2, j1) = dreal(psi_sum(j2))
          end do
          sum = 0.0d0
          do j2 = 1, n - 1
             sum = sum + 5.0d-1*(u(j2 + 1, j1)*u(j2 + 1, j1) + u(j2, j1)*u(j2, j1))
          end do
          sum = 1.0d0/dsqrt(dx*sum)
          do j2 = 1, n
             u(j2, j1) = sum*u(j2, j1)
          end do
       end do
       open(unit = 1, file = 'EITandEIF.dat', status = 'new', form = 'unformatted')
       write(1) x
       write(1) E
       write(1) u
       close(1)
       write(0, 20)
       write(0, 30)
       if (qshowEIF .eqv. .true.) then
          call plot_EIF(E, u)
       end if

       return
10     format(1x, 'Currently determining bound-state position-space energy eigenfunction for n =', 1x, i3, '...')
20     format(1x, 'All done calculating bound-state position-space energy eigenfunctions...')
30     format(1x)
       end subroutine                             EIF


!***********************************************************************************************************************************
! This Fast Fourier Transform (FFT) subroutine is based on the one by DeVries [DeVries, P. L., "A First Course in Computational
! Physics," John Wiley & Sons, Inc., New York, 1994), pp. 313-6], which is in turn based on the one by Cooley, Lewis, and Welch
! [Cooley, J W, Lewis P A W and Welch P D 1969 "The Fast Fourier Transform and its Applications," Institute of Electrical and
! Electronics Engineers (IEEE) Transactions on Education 12, pp. 27-34.<--Note: We *think* the similar reference DeVries put in his
! book isn't quite right...] We converted the sine and cosine pairs into the corresponding complex exponential functions; somehow
! deep in the "bowels" of the compiler this leads to machine instructions that substantially reduce "noise" and roundoff error.
! Also, this subroutine provides for normalization, though usually the cancellation that occurs when calculating the inverse FFT
! obviates it. If normalization isn't done, then the variable 'ds' still has to be given a value, but this value becomes irrelevant.
! The normalization is symmetric; both the transform integral and its inverse have the factor '1.0e0/sqrt(2.0e0)' premultiplying
! them. This is common in physics. 'jsign' is set to +/-1 according to the sign of 'i' in the kernel of the Fourier integral.
! ---cjs, last modified 12/8/2025
       subroutine                                 FFT(ds, jsign, m, qnormalize, psi)
       implicit                                   none
       logical ::                                 qnormalize
       integer(kind = 8) ::                       j1, j2, j3, j4, j5, j6, j7, jsign, m, n, nhalf
       real(kind = 8) ::                          ds, renorm, theta
       real(kind = 8), parameter ::               pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164063d0
       complex(kind = 8) ::                       u, w, temp
       complex(kind = 8), parameter ::            i = (0.0d0, 1.0d0)
       complex(kind = 8), dimension(2**m) ::      psi

       n     =  2**m
       nhalf =  n/2
       j2    =  1
       do j1 = 1, n - 1
          if (j1 < j2) then
             temp    = psi(j2)
             psi(j2) = psi(j1)
             psi(j1) = temp
          end if
          j3 = nhalf
10        if (j3 < j2) then
             j2 = j2 - j3
             j3 = j3/2
             go to 10
          end if
          j2 = j2 + j3
       end do
       j5 = 1
       if (jsign == -1) then
          do j4 = 1, m
             j6    =  j5
             j5    =  j5 + j5
             u     = (1.0d0, 0.0d0)
             theta =  pi/dble(j6)
             w     =  cdexp(-i*theta)
             do j2 = 1, j6
                do j1 = j2, n, j5
                   j7      = j1 + j6
                   temp    = psi(j7)*u
                   psi(j7) = psi(j1) - temp
                   psi(j1) = psi(j1) + temp
                end do
                u = u*w
             end do
          end do
       else if (jsign == 1) then
          do j4 = 1, m
             j6    =  j5
             j5    =  j5 + j5
             u     = (1.0d0, 0.0d0)
             theta =  pi/dble(j6)
             w     =  cdexp(i*theta)
             do j2 = 1, j6
                do j1 = j2, n, j5
                   j7      = j1 + j6
                   temp    = psi(j7)*u
                   psi(j7) = psi(j1) - temp
                   psi(j1) = psi(j1) + temp
                end do
                u = u*w
             end do
          end do
       else
          write(0, 20)
          stop
       end if
       if ((jsign == -1) .and. (qnormalize .eqv. .false.)) then
          renorm = 1.0d0
       else if ((jsign == 1) .and. (qnormalize .eqv. .false.)) then
          renorm = 1.0d0/dble(n)
       else if ((jsign == -1) .and. (qnormalize .eqv. .true.)) then
          renorm = ds/dsqrt(2.0d0*pi)
       else if ((jsign == 1) .and. (qnormalize .eqv. .true.)) then
          renorm = dsqrt(2.0d0*pi)/(dble(n)*ds)
       else
          write(0, 20)
          stop
       end if
       psi = renorm*psi
       return

20     format(1x, 'There was an error in calculating the FFT.')
       end subroutine                             FFT


!***********************************************************************************************************************************
! Here the graphics for propagation are intitialized, including the format for the frames displayed.
       subroutine                                 start_graphics(psi)
       use                                        plplot
       implicit                                   none
       logical, parameter ::                      qshowmovie = .false., qprob = .false.
       integer(kind = 8) ::                       j
       integer(kind = 4), parameter ::            nh = 3, nw = 3 ! unfortunately seem to have to use kind = 4 here...
       integer(kind = 8), parameter ::            m = 10, n = 2**m
       real(kind = 8) ::                          a, b, ymin, vinf, vmin, zzmax
       real(kind = 8), parameter ::               scale = 9.8d-1
       real(kind = 8), dimension(n) ::            ptntl, vv, zz
       complex(kind = 8), dimension(n) ::         psi
       common/ptl/                                ptntl
       common/wvfngraphics/                       ymin, vv, zzmax

       ymin = -1.0d0
       vinf =  ptntl(n)
       vmin =  1.0d99
       do j = 1, n
          if (ptntl(j) < vmin) then
             vmin = ptntl(j)
          end if
       end do
       a =  scale/(vinf - vmin)
       b = -a*vinf
       do j = 1, n
          vv(j) = a*ptntl(j) + b
       end do
       do j = 1, n
          if (vv(j) > -ymin) then
             vv(j) = -ymin
          end if
       end do
       zzmax  = 0.0d0
       if (qprob .eqv. .true.) then
          zz = dreal(dconjg(psi)*psi)
       else if (qprob .eqv. .false.) then
          zz = cdabs(psi)
       else
          write(0, 10)
          stop
       end if
       do j = 1, n
          if (zz(j) > zzmax) then
             zzmax = zz(j)
          end if
       end do
       call plsdev('wingcc')
       call plscolbg(255, 255, 255)
       call plscol0(15, 0, 0, 0)
       call plinit()
       if (qshowmovie .eqv. .true.) then
          call pladv(0)
       else if (qshowmovie .eqv. .false.) then
          call plssub(nh, nw)
       else
          write(0, 10)
          stop
       end if
       return

10     format(1x, 'There was an error in propagation graphics production.')
       end subroutine                             start_graphics


!***********************************************************************************************************************************
! Here each frame of the propagation is plotted individually.
       subroutine                                 plot_frame(jframe, jstep, psi)
       use                                        plplot
       implicit                                   none
       logical, parameter ::                      qshowmovie = .false., qprob = .false., qshowinit = .true.
       integer(kind = 8) ::                       jframe, jstep
       integer(kind = 4), parameter ::            nh = 3, nw = 3 ! unfortunately seem to have to use kind = 4 here...
       integer(kind = 8), parameter ::            m = 10, n = 2**m
       real(kind = 8) ::                          t, ymin, vinf, zzmax
       real(kind = 8), parameter ::               dt = 1.0d0, zztop = 5.0d-1
       real(kind = 8), dimension(n) ::            x, vv, zz, zzim, zzre, zzz, zzzim, zzzre
       real(kind = 8), dimension(nh*nw) ::        time
       complex(kind = 8), dimension(n) ::         psi, psi0
       data                                       time/0.0d0, 5.0d3, 1.0d4, 3.0d4, 5.0d4, 8.0d4, 1.0d5, 1.3d5, 1.69163d5/
       common/position/                           x
       common/wvfngraphics/                       ymin, vv, zzmax
       common/init/                               psi0

       if ((qshowmovie .eqv. .true.) .and. (qprob .eqv. .true.)) then
          zz  = zztop*dreal(dconjg(psi)*psi)/zzmax
          zzz = zztop*dreal(dconjg(psi0)*psi0)/zzmax
          call plvpor(0.0d0, 1.0d0, 0.0d0, 1.0d0)
          call plwind(x(1), x(n), ymin, -ymin)
          call plcol0(15)
          call plbop()
          call plbox('bc', 0.0d0, 0, 'bc', 0.0d0, 0)
          call plline(x, vv)
          call plline(x,  zz + vinf)
          call plcol0(1)
          if (qshowinit .eqv. .true.) then
             call plline(x,  zzz + vinf)
          end if
       else if ((qshowmovie .eqv. .false.) .and. (qprob .eqv. .true.)) then
          t = dble(jstep - 1)*dt
          if (t > time(jframe) - 5.0d-1*dt .and. t < time(jframe) + 5.0d-1*dt) then
             zz = zztop*dreal(dconjg(psi)*psi)/zzmax
             write(0, 10) jstep - 1
             call plcol0(15)
             call pladv(0)
             call plvpor(0.0d0, 1.0d0, 0.0d0, 1.0d0)
             call plwind(x(1), x(n), ymin, -ymin)
             call plbox('bc', 0.0d0, 0, 'bc', 0.0d0, 0)
             call plline(x, vv)
             call plline(x,  zz + vinf)
             jframe = jframe + 1
          else if (jframe > nh*nw) then
             return
          end if
       else if ((qshowmovie .eqv. .true.) .and. (qprob .eqv. .false.)) then
          zz    = zztop*cdabs(psi)/zzmax
          zzim  = zztop*dimag(psi)/zzmax
          zzre  = zztop*dreal(psi)/zzmax
          zzz   = zztop*cdabs(psi0)/zzmax
          zzzim = zztop*dimag(psi0)/zzmax
          zzzre = zztop*dreal(psi0)/zzmax
          call plvpor(0.0d0, 1.0d0, 0.0d0, 1.0d0)
          call plwind(x(1), x(n), ymin, -ymin)
          call plcol0(15)
          call plbop()
          call plbox('bc', 0.0d0, 0, 'bc', 0.0d0, 0)
          call plline(x, vv)
          call plcol0(3)
          call plline(x,  zz + vinf)
          call plline(x, -zz + vinf)
          call plcol0(1)
          call plline(x, zzim + vinf)
          call plcol0(9)
          if (qshowinit .eqv. .true.) then
             call plline(x, zzre + vinf)
             call plline(x,  zzz + vinf)
             call plline(x, -zzz + vinf)
             call plcol0(1)
             call plline(x, zzzim + vinf)
             call plcol0(9)
             call plline(x, zzzre + vinf)
          end if
       else if ((qshowmovie .eqv. .false.) .and. (qprob .eqv. .false.)) then
          t = dble(jstep - 1)*dt
          if (t > time(jframe) - 5.0d-1*dt .and. t < time(jframe) + 5.0d-1*dt) then
             zz   = zztop*cdabs(psi)/zzmax
             zzim = zztop*dimag(psi)/zzmax
             zzre = zztop*dreal(psi)/zzmax
             write(0, 10) jstep - 1
             call plcol0(15)
             call pladv(0)
             call plvpor(0.0d0, 1.0d0, 0.0d0, 1.0d0)
             call plwind(x(1), x(n), ymin, -ymin)
             call plbox('bc', 0.0d0, 0, 'bc', 0.0d0, 0)
             call plline(x, vv)
             call plcol0(3)
             call plline(x,  zz + vinf)
             call plline(x, -zz + vinf)
             call plcol0(1)
             call plline(x, zzim + vinf)
             call plcol0(9)
             call plline(x, zzre + vinf)
             jframe = jframe + 1
          else if (jframe > nh*nw) then
             return
          end if
       end if
       return

10     format(1x, 'At frame', 1x, i6)
       end subroutine                             plot_frame


!***********************************************************************************************************************************
! This subroutine plots the energy spectrum.
       subroutine                                 plot_EIT(spectrum)
       use                                        plplot
       implicit                                   none
       logical ::                                 qshowEIT = .true.
       integer(kind = 8) ::                       j
       integer(kind = 8), parameter ::            nEIT = 50, msteps = 18, nsteps = 2**msteps
       real(kind = 8) ::                          E_n, Smax, vmax
       real(kind = 8), dimension(nsteps) ::       energy, spectrum
       common/EITgraphics/                        energy, Smax, vmax

       if (qshowEIT .eqv. .false.) then
          return
       end if
       call plsdev('wingcc')
       call plscolbg(255, 255, 255)
       call plscol0(15, 0, 0, 0)
       call plinit()
       call pladv(0)
       call plvpor(0.0d0, 1.0d0, 0.0d0, 1.0d0)
       call plwind(0.0d0, 1.1d0*vmax, 0.0d0, 1.1d0*Smax)
       call plcol0(15)
       call plbox('bc', 0.0d0, 0, 'bc', 0.0d0, 0)
       call pllab('', '', '')
       call plcol0(10)
       do j = 1, nEIT
          call pljoin(E_n(j - 1), 0.0d0, E_n(j - 1), 1.1d0*Smax)
       end do
       call plcol0(15)
       call plline(energy, spectrum)
       call plend()
       return

       end subroutine                             plot_EIT


!***********************************************************************************************************************************
! Here the energy eigenfunctions are plotted along with the potential-energy function.
       subroutine                                 plot_EIF(E, u)
       use                                        plplot
       implicit                                   none
       integer(kind = 8) ::                       j1, j2
       integer(kind = 8), parameter ::            m = 10, n = 2**m, nEIT = 50, nshowEIT = 4
       integer(kind = 8), dimension(nshowEIT) ::  mEIT
       integer(kind = 4) ::                       jcolor ! seem to have to use kind = 4 here
       integer(kind = 4), dimension(nshowEIT) ::  jjcolor
       real(kind = 8) ::                          a, b, umax, vinf, vmin, ymin
       real(kind = 8), parameter ::               scale = 9.8d-1
       real(kind = 8), dimension(n) ::            psin, ptntl, x, vv
       real(kind = 8), dimension(nEIT) ::         E
       real(kind = 8), dimension(n, nEIT) ::      u
       data                                       jjcolor/1, 9, 3, 10/
       data                                       mEIT/1, 2, 27, 50/
       common/position/                           x
       common/ptl/                                ptntl

       ymin  = -1.0d0
       umax  =  0.0d0
       vinf  =  ptntl(n)
       vmin  =  1.0d99
       do j1 = 1, n
          if (ptntl(j1) < vmin) then
             vmin = ptntl(j1)
          end if
       end do
       a =  scale/(vinf - vmin)
       b = -a*vinf
       do j1 = 1, n
          vv(j1) = a*ptntl(j1) + b
       end do
       do j1 = 1, n
          if (vv(j1) > -ymin) then
             vv(j1) = -ymin
          end if
       end do
       do j1 = 1, nEIT
          do j2 = 1, n
             if (dabs(u(j2, j1)) > umax) then
                umax = dabs(u(j2, j1))
             end if
          end do
       end do
       umax = scale/umax
       call plsdev('wingcc')
       call plscolbg(255, 255, 255)
       call plscol0(15, 0, 0, 0)
       call plinit()
       call pladv(0)
       call plvpor(0.0d0, 1.0d0, 0.0d0, 1.0d0)
       call plwind(x(1), x(n), ymin, -ymin)
       call plcol0(15)
       call plbox('bc', 0.0d0, 0, 'bc', 0.0d0, 0)
       call pllab('', '', '')
       call plline(x, vv)
       jcolor = 1
       do j1 = 1, nshowEIT
          call plcol0(jjcolor(jcolor))
          do j2 = 1, n
             psin(j2) = umax*u(j2, mEIT(j1))
          end do
          call plline(x, psin)
          jcolor = jcolor + 1
       end do
       call plend()

       return
       end subroutine                             plot_EIF


!***********************************************************************************************************************************
! Here the bound-state eigenenergies are determined from a known analytic formula.
       function                                   E_n(j)
       implicit                                   none
       integer(kind = 8) ::                       j
       real(kind = 8) ::                          E_n
       real(kind = 8), parameter ::               D = 5.211d0*1.6021766d-19/4.35944d-18, a = 2.78d0*5.29167d-11/1.0d-10, &
                                                  mu = 5.0d-1*2.6566963d-26/9.1093837d-31

       E_n = a*dsqrt(2.0d0*D/mu)*(dble(j) + 5.0d-1) - 5.0d-1*a*a*(dble(j) + 5.0d-1)*(dble(j) + 5.0d-1)/mu

       return
       end function                               E_n
