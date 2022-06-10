program main
    
	implicit none
    ! definition of parameters
	! ep: small number
	! g: gravitational coefficient
	! R0: submerged specifuc gravity
	! na: exponent of MPM equation modified by Wong and Parker 2006
	! alpr: coefficient of MPM equation modified by Wong and Parker 2006
	! taucstar: nominal critical value of dimensionless Shields number
	! Cf: dimensionless bed resistance coefficient
	! cr: denotes 0.75
	! one: denotes 1
	! half: denotes 0.5
	! zero: denotes 0
	! div: When landslides occur every 'div' times, an output file is generated
	! num: the number that landslide cycle repeats
	real(8),parameter::     ep=1E-8_8,  g=9.8_8,   R0=1.65_8,  &
                na=1.5_8,   alpr=4.0_8, taucstar=0.0495_8,     &
                Cf=0.01_8,  cr=0.75_8,  one=1.0_8, half=0.5_8, &
                zero=0.0_8
    integer(8),parameter::  none = 0,   div = 80
    integer(8):: i, j, k, num
    character(len=64):: filename
    integer::myistat
    ! definition of initial and boundary conditions
	! B: channel width (m)
	! D: characteristic grain size (m)
	! r: porosity of alluvial deposit
	! L: computing channel length (m)
	! L0: channel length without which the normal flow assumption is implementd in decadal scale (m)
	! Lmr: bedrock macro-roughness (m)
	! v0: rock uplift rate (mm/year)
	! bss: bedrock abrasion coefficient
	! Sb0: initial bedrock slope
    real(8)::               &
        B   =   100.0_8,    &
        D   =   0.02_8,     &
        r   =   0.35_8,     &
	    L   =   60E3_8,     &
        L0  =   60E3_8,     &
        Lmr =   1.0_8,      &
        v0  =   5.0_8,      &
        bss =   0.05E-3_8,  &
        Sb0 =   0.006_8                 
    ! definition of temporal-spatial grid size
	! dt: time step in flood timescale (s)
	! dt2: time step in days to years timescale (s)
	! dt3: time step in decadal timescale (s)
	! dx: grid size in flood timescale (m)
	! dx2: grid size in days to years timescale (m)
	! dx3: grid size in decadal timescale (m)
	! dxs: grid size for calculation of backwater equation (m)
	! n, n2, n3, nb, ns: number of grids in different timescales
	! dns, dn2, dns2, dn3, dn23: ratio of different grid numbers (used for interpolation)
    real(8)::               &
        dt  =   0.5_8,      &
        dt2 =   1200.0_8,   &
        dt3 =   3600.0_8,   &
        dx  =   25.0_8,     &
        dx2 =   50.0_8,     &
        dx3 =   100.0_8,    &
        dxs =   5.0_8
    integer(8)::                    &
        n, n2, n3, nb, ns,          &
        dns, dn2, dns2, dn3, dn23
    ! definition of water and sediment supply under initial equilibrium state
	! hss: initial equilibrium alluvial thickness (m)
	! qaf: sediment feed rate per unit width at the inlet (m2/s)
	! Q0: water discharge at the inlet (m3/s)
	! qw0: water discharge per unit width at the inlet (m2/s)
	! hw0: water depth at the inlet (m)
	! t: time
	! Td: time interval of landslide dam occurrence (year)
	! deadline, deadline2, deadline3: computing deadline
    real(8)::               &
        hss =   0.2_8,      &
        qaf,                &
        Q0,                 &
        qw0,                &
        hw0,                &
        t   =   0.0_8,      &
        Td  =   100.0_8,    &
        deadline, deadline2, deadline3
	! ha: alluvial thickness
	! h_dam: dam height
	! hb: bedrock elevation
	! h: bed elevation
	! Sb: bedrock slope
	! S: bed slope
	! hw: water depth
	! hs: water depth corresponding to backwater grid size
	! Ss: bed slope corresponding to backwater grid size
	! hws: a temporary value of water depth
	! qw: discharge per unit width
	! qwl, qwr, hl, hr, ul, ur: left&right condition of shallow water equation
	! cal, car, ustar, cstar, sl, sr: wave speed of the HLL method
	! FLu, FLd, FRu, FRd, FHu, FHd, Fu, Fd: numerical flux of HLL method
	! rs, rf: source term of HLL method
	! taustar: dimensionless Shields number
	! qac: capacity volume sediment transport rate per unit width
	! p: alluvial cover fraction
	! pa: adjusted alluvial cover fraction
	! ha2, hb2, h2, S2, hw2: parameters corresponding to days to years timescale
	! ha3, hb3, h3, S3: parameters corresponding to decadal timescale
    real(8),dimension(:),allocatable::  &
        ha, h_dam, hb, h, Sb, S,        &
        hw, hs, Ss, hws,                &
        qw, qwl, qwr, hl, hr, ul, ur,   &
        cal, car, ustar, cstar, sl, sr, &
        FLu, FLd, FRu, FRd, FHu, FHd,   &
        Fu, Fd, rs, rf,                 &
        taustar, qac, p, pa,            &
        ha2, hb2, h2, S2, hw2,          &
        ha3, hb3, h3, S3
    ! Intermediate variables in the calculation of backwater equation
    real(8):: hc, hn, FR1, FR2, Sf, hp, ua, Fra, ustara, sfa, dha
    
    
    ! calculation of the initial equilibrium state
    v0  = v0/1000.0_8/(365.0_8*24.0_8*3600.0_8)
    qaf = v0/(one-hss*0.9_8/0.95_8)/bss;
    Q0  = B*SQRT((((qaf/(hss*0.9_8/0.95_8)/4.0_8/D/SQRT(R0*g*D))**(2.0_8/3.0_8)+taucstar)*R0*D/Sb0**(2.0_8/3.0_8))**3.0_8/Cf*g)
    qw0 = Q0/B
    hw0 = (Cf*qw0*qw0/g/Sb0)**(one/3.0_8)
    ! calculation of the grid numbers
    n  = L/dx
    n2 = L/dx2
    n3 = L/dx3
    nb = L0/dxs
    ns = L/dxs
    dns  = ns/n
    dn2  = n/n2
    dns2 = ns/n2
    dn3  = n/n3
    dn23 = n2/n3
    
    ! set the initial topographic data
    allocate(ha(n+1), h_dam(n+1), hb(n+1), h(n+1), Sb(n+1), S(n+1))
    ha = hss
    h_dam = (/ (0.0, i=1,1120), &
        0.0, 12.5, 25.0, 37.5, 50.0, 62.5, 75.0, 87.5, 100.0, &
        (100.0, i=1,11), &
        100.0, 87.5, 75.0, 62.5, 50.0, 37.5, 25.0, 12.5, 0.0, &
        (0.0, i=1,18052) /)
    ha = ha+h_dam
    do i = 0,n
        hb(i+1) = (n-i)*Sb0*dx
    enddo
    h  = ha+hb
    Sb = (/ ((hb(1:n)-hb(2:n+1))/dx), ((hb(n)-hb(n+1))/dx) /)
    S = (/ ((h(1:n)-h(2:n+1))/dx), ((h(n)-h(n+1))/dx) /)
    allocate(ha3(n3+1), hb3(n3+1), h3(n3+1), S3(n3+1))
    
	! Start calculation
    do num = 1, 800
        deadline  = 10*24*3600+365*24*3600*Td*(num-1)
        deadline2 = 10*365*24*3600+365*24*3600*Td*(num-1)
        deadline3 = Td*365*24*3600*num
        ! flood timescale
        ! Interpolate to fine mesh corresponding to backwater equation
        allocate(hs(ns+1), Ss(ns+1), hws(ns+1))
        hs(1) = h(1)
        do i = 1,n
            do j = 1,dns
                hs((i-1)*dns+j+1) = one*(dns-j)/dns*h(i)+one*j/dns*h(i+1)
            enddo
        enddo
        Ss(1) = (hs(1)-hs(2))/dxs
        Ss(2:ns) = (hs(1:ns-1)-hs(3:ns+1))/dxs*half
        Ss(ns+1) = (hs(ns)-hs(ns+1))/dxs
        ! calculate the initial water depth according to backwater equation
        hws(ns+1) = (Cf*qw0*qw0/g/Ss(ns+1))**(one/3.0_8)
        hc = (qw0*qw0/g/cr/cr)**(one/3.0_8)
        do k = 1,ns
            if (Ss(ns+1-k) > 0) then
                hn = (Cf*qw0*qw0/g/Ss(ns+1-k))**(one/3.0_8)
            else
                hn = hc
            endif
            FR2 = qw0/SQRT(g*hws(ns+2-k)*hws(ns+2-k)*hws(ns+2-k))
            if (FR2 <= cr+ep) then
                Sf = Cf*FR2*FR2
                hp = hws(ns+2-k)-(Ss(ns+2-k)-Sf)/(one-FR2*FR2)*dxs
                hws(ns+1-k) = hws(ns+2-k)-half*((Ss(ns+2-k)-Sf)/(one-FR2*FR2)+(Ss(ns+2-k)-Cf*qw0*qw0/(g*hp*hp*hp))/(one-qw0*qw0/(g*hp*hp*hp)))*dxs
            else
                hws(ns+1-k) = hn
            endif
            FR1 = qw0/SQRT(g*hws(ns+1-k)*hws(ns+1-k)*hws(ns+1-k))
            if (FR1<=cr+ep .AND. FR2>cr+ep) then
                ua  = qw0/hc
                Fra = ua/SQRT(g*hc)
                sfa = Cf*ua*ua/g/hc
                dha = (Ss(ns+2-k)-sfa)/(one-Fra*Fra)
                hws(ns+1-k) = hc-dha*dxs
                FR1 = qw0/SQRT((g*hws(ns+1-k)*hws(ns+1-k)*hws(ns+1-k)))
            endif
            if (Ss(ns+1-k)>=0 .AND. Ss(ns+2-k)<0) then
                Sf = Cf*FR2*FR2
                hp = hws(ns+2-k)-(Ss(ns+1-k)-Sf)/(one-FR2*FR2)*dxs
                hws(ns+1-k) = hws(ns+2-k)-half*((Ss(ns+1-k)-Sf)/(one-FR2*FR2)+(Ss(ns+1-k)-Cf*qw0*qw0/(g*hp*hp*hp))/(one-qw0*qw0/(g*hp*hp*hp)))*dxs
            endif
            if (FR1 > cr+ep) then
                hws(ns+1-k) = hn
            endif
            if (hws(ns+1-k) < ep) then
                hws(ns+1-k) = ep
            endif
        enddo
        ! Interpolate to meshes corresponding to flood timescale
        allocate(hw(n+1))
        hw(1) = hws(1)
        do i = 1,n
            hw(i+1) = hws(i*dns+1)
        enddo
        deallocate(hs, Ss, hws)
        
        ! Hydraulics£ºShallow Water Equation
        allocate(qw(n+1), hl(n+2), hr(n+2), qwl(n+2), qwr(n+2), ul(n+2), ur(n+2), &
                cal(n+2), car(n+2), ustar(n+2), cstar(n+2), sl(n+2), sr(n+2), &
                FLu(n+2), FLd(n+2), FRu(n+2), FRd(n+2), FHu(n+2), FHd(n+2), Fu(n+2), Fd(n+2), &
                rs(n+1), rf(n+1), taustar(n+1), qac(n+1), p(n+1), pa(n+1))
        qw = qw0
        do while (t<deadline)
            ! left&right conditions
            qwl = (/ qw0,   qw      /)
            qwr = (/ qw,    qw(n+1) /)
            hl  = (/ hw0,   hw      /)
            hr  = (/ hw,    hw(n+1) /)
            ul  = qwl/hl
            ur  = qwr/hr
            ! wave speed
            cal = SQRT(g*hl)
            car = SQRT(g*hr)
            ustar = half*(ul+ur)+cal-car
            cstar = half*(cal+car)+0.25_8*(ul-ur)
            sl  = min(ul-cal,ustar-cstar)
            sr  = max(ur+car,ustar+cstar)
            ! numerical flux
            FLu = qwl
            FLd = qwl*qwl/hl+half*g*hl*hl
            FRu = qwr
            FRd = qwr*qwr/hr+half*g*hr*hr
            FHu = (sr*FLu-sl*FRu+sl*sr*(hr-hl))/(sr-sl)
            FHd = (sr*FLd-sl*FRd+sl*sr*(qwr-qwl))/(sr-sl)
            do k = 1, n+2
                if (sr(k)<0) then
                    Fu(k)  = FRu(k)
                    Fd(k)  = FRd(k)
                else if (sl(k)>0) then
                    Fu(k)  = FLu(k)
                    Fd(k)  = FLd(k)
                else
                    Fu(k)  = FHu(k)
                    Fd(k)  = FHd(k)
                endif
            enddo
            ! source term
            rs = g*hw*S
            rf = Cf*abs(Fu(1:n+1))/(hw*hw)
            ! time advance
            hw = hw+(Fu(1:n+1)-Fu(2:n+2))/dx*dt
            qw = (qw+(Fd(1:n+1)-Fd(2:n+2))/dx*dt+rs*dt)/(one+rf*dt)
            
            ! calculating the variation of alluvial cover and bedrock elevation with MRSAA model
            taustar = Cf*Fu(1:n+1)*Fu(1:n+1)/(hw*hw*R0*g*D)
            p  = 0.05_8+0.9_8*ha/Lmr
            pa = 0.9_8*ha/Lmr/0.95_8
            do i = 1, n+1
                if (taustar(i) >= taucstar) then
                    qac(i) = alpr*SQRT(R0*g*D)*D*(taustar(i)-taucstar)**na
                else
                    qac(i) = zero
                endif
                if (p(i) > one) then
                    p(i) = one
                endif
                if (pa(i) > one) then
                    pa(i) = one
                endif
            enddo
            ! variation of alluvium thickness
            ha(1) = ha(1)-dt/(one-r)/p(1)*(pa(1)*qac(1)-qaf)/dx
            ha(2:n+1) = ha(2:n+1)-dt/(one-r)/p(2:n+1)*(pa(2:n+1)*qac(2:n+1)-pa(1:n)*qac(1:n))/dx
            if (ha(n+1) > Lmr) then
                ha(n+1) = Lmr
            endif
            ! variation of bedrock elevation
            hb(1:n) = hb(1:n)+dt*(v0-bss*qac(1:n)*pa(1:n)*(one-pa(1:n)))
            hb(n+1) = zero
            
            do i = 1, n+1
                if (ha(i) < zero) then
                    ha(i) = zero
                endif
            enddo
            h = ha+hb
            S = (/ ((h(1:n)-h(2:n+1))/dx), ((h(n)-h(n+1))/dx) /)
            
            t = t+dt
        enddo
        deallocate(qw, hl, hr, qwl, qwr, ul, ur,        &
                cal, car, ustar, cstar, sl, sr,         &
                FLu, FLd, FRu, FRd, FHu, FHd, Fu, Fd,   &
                rs, rf, taustar, qac, p, pa, hw)
        
        ! days to years timescale
        ! Interpolate to meshes corresponding to days to years timescale
        allocate(ha2(n2+1), hb2(n2+1), h2(n2+1))
        ha2(1) = ha(1)
        hb2(1) = hb(1)
        do i = 1,n2
            ha2(i+1) = ha(i*dn2+1)
            hb2(i+1) = hb(i*dn2+1)
        enddo
        h2 = ha2+hb2
        
        allocate(hs(ns+1), Ss(ns+1), hws(ns+1))
        allocate(hw2(n2+1), taustar(n2+1), qac(n2+1), p(n2+1), pa(n2+1))
        do while (t<deadline2)
            ! Interpolate to meshes corresponding to backwater equation
            hs(1) = h2(1)
            do i = 1,n2
                do j = 1,dns2
                    hs((i-1)*dns2+j+1) = one*(dns2-j)/dns2*h2(i)+one*j/dns2*h2(i+1)
                enddo
            enddo
            Ss(1) = (hs(1)-hs(2))/dxs
            Ss(2:ns) = (hs(1:ns-1)-hs(3:ns+1))/dxs*half
            Ss(ns+1) = (hs(ns)-hs(ns+1))/dxs
            ! calculating water depth by backwater equation
            hws(nb+1:ns+1) = (Cf*qw0*qw0/g/Ss(nb+1:ns+1))**(one/3.0_8)
            do k = 1,nb
                if (Ss(nb+1-k) > 0) then
                    hn = (Cf*qw0*qw0/g/Ss(nb+1-k))**(one/3.0_8)
                else
                    hn = hc
                endif
                FR2 = qw0/SQRT(g*hws(nb+2-k)*hws(nb+2-k)*hws(nb+2-k))
                if (FR2 <= 0.95_8) then
                    Sf = Cf*FR2*FR2
                    hp = hws(nb+2-k)-(Ss(nb+2-k)-Sf)/(one-FR2*FR2)*dxs
                    hws(nb+1-k) = hws(nb+2-k)-half*((Ss(nb+2-k)-Sf)/(one-FR2*FR2)+(Ss(nb+2-k)-Cf*qw0*qw0/(g*hp*hp*hp))/(one-qw0*qw0/(g*hp*hp*hp)))*dxs
                else
                    hws(nb+1-k) = hn
                endif
                FR1 = qw0/SQRT(g*hws(nb+1-k)*hws(nb+1-k)*hws(nb+1-k))
                if (FR1<=0.95_8 .AND. FR2>0.95_8) then
                    ua  = qw0/hc
                    Fra = ua/SQRT(g*hc)
                    ustara = ua/8.1_8*(1.75_8*D/hc)**(one/6.0_8)
                    sfa = ustara*ustara/g/hc
                    dha = (Ss(nb+1-k)-sfa)/(one-Fra*Fra)
                    hws(nb+1-k) = hc-dha*dxs
                    FR1 = qw0/SQRT((g*hws(nb+1-k)*hws(nb+1-k)*hws(nb+1-k)))
                endif
                if (FR1 > 0.95_8) then
                    hws(nb+1-k) = hn
                endif
                if (hws(nb+1-k) < ep) then
                    hws(nb+1-k) = ep
                endif
            enddo
            ! Interpolate to meshes corresponding to days to years timescale
            hw2(1) = hws(1)
            if (hw2(1) < 0.2_8) then
                hw2(1) = 0.2_8
            endif
            do i = 1,n2
                hw2(i+1) = hws(i*dns2+1)
                if (hw2(i+1) < 0.2_8) then
                    hw2(i+1) = 0.2_8
                endif
            enddo
            ! calculating the variation of alluvial cover and bedrock elevation with MRSAA model
            taustar = Cf*qw0*qw0/(hw2*hw2*R0*g*D)
            p  = 0.05_8+0.9_8*ha2/Lmr
            pa = 0.9_8*ha2/Lmr/0.95_8
            do i = 1, n2+1
                if (taustar(i) >= taucstar) then
                    qac(i) = alpr*SQRT(R0*g*D)*D*(taustar(i)-taucstar)**na
                else
                    qac(i) = zero
                endif
                if (p(i) > one) then
                    p(i) = one
                endif
                if (pa(i) > one) then
                    pa(i) = one
                endif
            enddo
            ha2(1) = ha2(1)-dt2/(one-r)/p(1)*(pa(1)*qac(1)-qaf)/dx2
            ha2(2:n2+1) = ha2(2:n2+1)-dt2/(one-r)/p(2:n2+1)*(pa(2:n2+1)*qac(2:n2+1)-pa(1:n2)*qac(1:n2))/dx2
            if (ha2(n2+1) > Lmr) then
                ha2(n2+1) = Lmr
            endif
            hb2(1:n2) = hb2(1:n2)+dt2*(v0-bss*qac(1:n2)*pa(1:n2)*(one-pa(1:n2)))
            hb2(n2+1) = zero
            
            do i = 1, n2+1
                if (ha2(i) < zero) then
                    ha2(i) = zero
                endif
            enddo
            h2 = ha2+hb2
            
            t = t+dt2
        enddo
        deallocate(hw2, taustar, qac, p, pa)
        deallocate(hs, Ss, hws)
        
        ! decadal timescale
        ! Interpolate to meshes corresponding to decadal timescale
        ha3(1) = ha2(1)
        hb3(1) = hb2(1)
        do i = 1,n3
            ha3(i+1) = ha2(i*dn23+1)
            hb3(i+1) = hb2(i*dn23+1)
        enddo
        h3 = ha3+hb3
        S3 = (/ ((h3(1:n3)-h3(2:n3+1))/dx3), ((h3(n3)-h3(n3+1))/dx3) /)
        deallocate(ha2, hb2, h2)
        allocate(taustar(n3+1), qac(n3+1), p(n3+1), pa(n3+1))
        do while (t<deadline3)
            taustar = (Cf*qw0*qw0*S3*S3/g)**(one/3.0_8)/R0/D
            p  = 0.05_8+0.9_8*ha3/Lmr
            pa = 0.9_8*ha3/Lmr/0.95_8
            do i = 1, n3+1
                if (taustar(i) >= taucstar) then
                    qac(i) = alpr*SQRT(R0*g*D)*D*(taustar(i)-taucstar)**na
                else
                    qac(i) = zero
                endif
                if (S3(i) < zero) then
                    qac(i) = zero
                endif
                if (p(i) > one) then
                    p(i) = one
                endif
                if (pa(i) > one) then
                    pa(i) = one
                endif
            enddo
            ! variation of alluvium thickness
            ha3(1) = ha3(1)-dt3/(one-r)/p(1)*(pa(1)*qac(1)-qaf)/dx3
            ha3(2:n3+1) = ha3(2:n3+1)-dt3/(one-r)/p(2:n3+1)*(pa(2:n3+1)*qac(2:n3+1)-pa(1:n3)*qac(1:n3))/dx3
            if (ha3(n3+1) > Lmr) then
                ha3(n3+1) = Lmr
            endif
            ! variation of bedrock elevation
            hb3(1:n3) = hb3(1:n3)+dt3*(v0-bss*qac(1:n3)*pa(1:n3)*(one-pa(1:n3)))
            hb3(n3+1) = zero
            
            do i = 1, n3+1
                if (ha3(i) < zero) then
                    ha3(i) = zero
                endif
            enddo
            h3 = ha3+hb3
            S3 = (/ ((h3(1:n3)-h3(2:n3+1))/dx3), ((h3(n3)-h3(n3+1))/dx3) /)
            
            t = t+dt3
        enddo
        deallocate(taustar, qac, p, pa)
        
        ! Interpolate to meshes corresponding to flood timescale
        ha(1) = ha3(1)
        hb(1) = hb3(1)
        do i = 1,n3
            do j = 1,dn3
                ha((i-1)*dn3+j+1) = one*(dn3-j)/dn3*ha3(i)+one*j/dn3*ha3(i+1)
                hb((i-1)*dn3+j+1) = one*(dn3-j)/dn3*hb3(i)+one*j/dn3*hb3(i+1)
            enddo
        enddo
        ! A landslide occurs
        ha = ha+h_dam
        h  = ha+hb
        S  = (/ ((h(1:n)-h(2:n+1))/dx), ((h(n)-h(n+1))/dx) /)
        ! output
        if ( mod(num, div) == none ) then
            write(filename,'(A,I5.5,A)')'ha_',num,'.txt'
            open(unit=2, file=trim(filename),status='replace',form='formatted',IOSTAT=myistat)
            if(myistat/=0) then
                write(*,*)'cannot open file: '//trim(filename)
            else
                do i = 1, n3+1
                    write(2,*)ha3(i)
                enddo
            endif
            close(unit=2,IOSTAT=myistat)

            write(filename,'(A,I5.5,A)')'hb_',num,'.txt'
            open(unit=3, file=trim(filename),status='replace',form='formatted',IOSTAT=myistat)
            if(myistat/=0) then
                write(*,*)'cannot open file: '//trim(filename)
            else
                do i = 1, n3+1
                    write(3,*)hb3(i)
                enddo
            endif
            close(unit=3,IOSTAT=myistat)
        endif
        
    enddo
    
end program main