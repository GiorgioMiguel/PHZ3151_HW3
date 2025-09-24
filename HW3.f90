PROGRAM HW3

IMPLICIT NONE

    ! declariation of variables and arrays
    integer :: P, T, i
    real :: Pr, Pc = 73.8, Tr, Tc = 304, V, Vr, Vc, Vi, Z, B, C, D, E, F, R = 83.14467, alpha = -1.66727022D-02, beta = 1.398D0, gamma = 2.96D-02, T_start = 200, T_end = 445, error, V_lower, V_upper, V_middle, f_lower, f_upper, f_middle, new_mid, old_mid, z1, z2
    real :: a1 = 8.99288497D-02, a2 = -4.94783127D-01, a3 = 4.77922245D-02, a4 = 1.03808883D-02, a5 = -2.82516861D-02, a6 = 9.49887563D-02, a7 = 5.20600880D-04, a8 = -2.93540971D-04, a9 = -1.77265112D-03, a10 = -2.51101973D-05, a11 = 8.93353441D-05, a12 = 7.88998563D-05
    real :: temps(50)
    real, parameter :: pressures(25) = (/0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0/)
    
    ! initialize constant variable
    Vc = R * (Tc/Pc)

    ! populating temps array
    do i = 1, 50
        temps(i) = T_start + (i - 1) * 5
    enddo

    ! need to solve for V at each location in the grid, 25x50 = 1,250 points on the grid. 
    do P = 1, 25
        do T = 1, 50

            ! calculating pressure and temperature dependent components
            Pr = pressures(P) / Pc
            Tr = temps(T) / Tc
            B = a1 + (a2 / Tr**2) + (a3 / Tr**3)
            C = a4 + (a5 / Tr**2) + (a6 / Tr**3)
            D = a7 + (a8 / Tr**2) + (a9 / Tr**3)
            E = a10 + (a11 / Tr**2) + (a12 / Tr**3)
            F = alpha / Tr**3

            ! guess upper and lower bounds of V for the bisection method
            Vi = R * temps(T) / pressures(P)
            V_lower = 0.1 * Vi
            V_upper = 10 * Vi

            ! solve for initial f(v_lower) and f(v_upper)
            Vr = V_lower / Vc
            z1 = (Pr * Vr) / Tr
            z2 = 1 + (B / Vr) + (C / Vr**2) + (D / Vr**4) + (E / Vr**5) + (F / Vr**2) * (beta + (gamma / Vr**2)) * exp(-gamma / Vr**2)
            f_lower = z2 - z1
            Vr = V_upper / Vc
            z1 = (Pr * Vr) / Tr
            z2 = 1 + (B / Vr) + (C / Vr**2) + (D / Vr**4) + (E / Vr**5) + (F / Vr**2) * (beta + (gamma / Vr**2)) * exp(-gamma / Vr**2)
            f_upper = z2 - z1

            ! testing to make sure bounds bracket a root, if not keep expanding until the do
            do while(f_lower * f_upper > 0) 

                ! adjust bounds
                V_lower = V_lower * 0.5
                V_upper = V_upper * 2

                ! solve for f(v_lower)
                Vr = V_lower / Vc
                z1 = (Pr * Vr) / Tr
                z2 = 1 + (B / Vr) + (C / Vr**2) + (D / Vr**4) + (E / Vr**5) + (F / Vr**2) * (beta + (gamma / Vr**2)) * exp(-gamma / Vr**2)
                f_lower = z2 - z1
                
                ! solve for f(v_upper)
                Vr = V_upper / Vc
                z1 = (Pr * Vr) / Tr
                z2 = 1 + (B / Vr) + (C / Vr**2) + (D / Vr**4) + (E / Vr**5) + (F / Vr**2) * (beta + (gamma / Vr**2)) * exp(-gamma / Vr**2)
                f_upper = z2 - z1
            enddo 

            ! initialize components
            error = 1
            do while(error > 1e-6)
                ! solve for f(v_middle)
                V_middle = (V_lower + V_upper) / 2
                old_mid = V_middle
                Vr = V_middle / Vc 
                z1 = (Pr * Vr) / Tr
                z2 = 1 + (B / Vr) + (C / Vr**2) + (D / Vr**4) + (E / Vr**5) + (F / Vr**2) * (beta + (gamma / Vr**2)) * exp(-gamma / Vr**2)
                f_middle = z2 -z1

                ! adjust bounds on location of root
                if(f_lower * f_middle < 0) then
                    V_upper = V_middle
                    f_upper = f_middle
                else if(f_lower * f_middle > 0) then
                    V_lower = V_middle
                    f_lower = f_middle  
                else if(f_lower * f_middle == 0) then
                    exit
                endif

                ! solve for new mid and calculate the current error
                new_mid = (V_lower + V_upper) / 2
                error = ABS((new_mid - old_mid) / new_mid)
            enddo

        enddo

    enddo
    
! Useful information:
! z1 = (P * V) / (R * T)
! z2 = (Pr * Vr) / Tr
! z3 = 1 + (B / Vr) + (C / Vr**2) + (D / Vr**4) + (E / Vr**5) + (F / Vr**2) * (beta + (gamma / Vr**2)) * exp(-gamma / Vr**2) 
! f(V) = z3 - z1

! a) Using the Bisection Method, solve for V in the equation of state over the entire grid of pressures and
!    temperatures. For each value of V compute Z. Finally, compute the degree of non-ideality (β = 1/Z).
!    The further away from “1” β is, the more non-ideal is the gas at those conditions.

! b) Use first order forward and second order central differencing schemes to compute the partial
!    derivatives, and , respectively, for the entire grid.

! c) Compute the heat capacity using the Shomate Equation using the above equation. Do not change
!    units. The heat capacity for a gas at constant pressure can be defined with the following integral:
!    Here Cp(Po) is the heat capacity at zero pressure (Po). You may use the values you just computed using
!    the Shomate Equation in part (c).

! d) Using your values for computed in part b, calculate Cp for the entire temperature-pressure grid.
!    For this problem, make sure that the units are consistent on both left and right hand sides. Remember,
!    your answer should be in cal/mol/K!

! e) Tabulate then plot β for the following 5 pressure levels over the entire temperature grid (0.1 bars, 1
!    bar, 10 bars, 50 bars, and 70 bars). What happens to β at different combinations of high/low pressure
!    and temperature? What does this suggest about how ideal and/or compressible the flow is as a function
!    of the pressure-temperature regime?

END PROGRAM HW3