program test_statistics
   implicit none
   integer::M
   integer::i
   integer :: io
   real*8::mind
   real*8::maxbd
   real*8::mintd
   real*8::maxtd
   integer::Tno
   real*8::avtd
   character, parameter :: nl = new_line('a')
   real*8, parameter :: bdef = 1.0d0

   real*8, allocatable :: Y(:)
   real*8, allocatable :: DY(:)
   open (newunit=io, file="input.txt", status="old", action="read")

   read (io, *) M
   allocate (Y(M))
   allocate (DY(M))

   do i = 1, M
      read (io, *) Y(i)
   end do

   close (io)

   call l_stat(Y, DY, M, bdef, mind, maxbd, mintd, maxtd, avtd)
   ! print *, "Y", nl, Y

   open (newunit=io, file="dy_s_fortran.txt", status="old", action="write")
      write (io, *) DY
   close (io)
   print *, nl, nl, "********Surface statistics (MSI)********", nl, nl
   print *, "M =", M, "bdef =", bdef, "mind =", mind, "maxbd =", maxbd, "mintd =", mintd, "maxtd =", maxtd, "avtd =", avtd
   print *, nl, nl, "********Surface statistics (MSII)********", nl, nl
   call h_stat(DY, M, bdef)
   deallocate (Y)
   deallocate (DY)
end program test_statistics

subroutine l_stat(Y, DY, M, bdef, mind, maxbd, mintd, maxtd, avtd)
   integer M
   Real*8 Y(*), DY(M)
!        Bunch / Terrace Definitions:
   real*8 bdef
!
   real*8 w, h
   real*8 stime
   real*8 dyi, mind, maxbd, mintd, maxtd
   real*8 avtd, tw, ntd
   integer Tno, i
   real*8 one, zero
   data one, zero/1.0d0, 0.0d0/

   character*21 tname

   real*8 xk

!
!                       ...initializations...
   w = zero
   h = zero
   tw = zero
   ntd = zero
   mind = bdef
   maxbd = zero
   mintd = M
   maxtd = zero
!

!          Compute the interstep distances, the minimal one and
!          some other quantities that characterize the bunching:
   do i = 2, M
      dyi = Y(i) - Y(i - 1)
      dy(i) = dyi
      if (dyi .le. bdef) then
!     counting the number of distances in bunches ...
         h = h + one
!                    ...and the total bunch width
         w = w + dyi

         if (dyi .lt. mind) mind = dyi
         if (dyi .gt. maxbd) maxbd = dyi
      else
         tw = tw + dyi
         ntd = ntd + one
         if (dyi .lt. mintd) mintd = dyi
         if (dyi .gt. maxtd) maxtd = dyi
      end if
   end do
   dyi = Y(1) - Y(M) + M
   dy(1) = dyi
   if (dyi .le. bdef) then
      w = w + dyi
      h = h + one
      if (dyi .lt. mind) mind = dyi
      if (dyi .gt. maxbd) maxbd = dyi
   else
      tw = tw + dyi
      ntd = ntd + one
      if (dyi .lt. mintd) mintd = dyi
      if (dyi .gt. maxtd) maxtd = dyi
   end if
   avtd = tw/ntd
   if (maxbd .gt. bdef) then
      write (*, *) 'Alert:maxbd!'
      stop
   end if
   if (mintd .lt. bdef) then
      write (*, *) 'Alert:mintd!'
      stop
   end if
!
!
!**********************************************************************************************
!      Counting the terraces and their width:

   Tno = 0
   if (dy(1) .gt. bdef) Tno = 1
   do i = 2, M
      dyi = dy(i)
      if (dyi .gt. bdef) then
         if (dy(i - 1) .le. bdef) Tno = Tno + 1
      end if
   end do
!              check if the terrace in the beginning spannes the boundaries
   if ((dy(1) .gt. bdef) .and. (dy(M) .gt. bdef)) Tno = Tno - 1
!
   if (Tno .gt. 0) then
      stime = dfloat(M)/dfloat(Tno)
      tw = tw/dfloat(Tno)
      h = h/dfloat(Tno)
      w = w/dfloat(Tno)
!*****
   else
      print *, 'terrace alert!'
   end if
3075 format(4(f16.8, 1x), e12.4)
3077 format(e15.8, 1x, 4(f12.4, 1x))
!
end subroutine l_stat

subroutine h_stat(DY, M, bdef)
   integer M
   Real*8 DY(M)
!        Bunch and Terrace Definitions:
   real*8 bdef
   real*8 xk
   integer one, zero
   data one, zero/1, 0/
   integer bno
   integer i, inc, bsz(M/2), ibw, bszi
   real*8 dyi, dyip
!
!         Average quantities:
   real*8 bw(M/2), bd(M/2), mindi(M/2), tbw, tmin, fbd(M/2), lbd(M/2), tla, stime
!          c
   include 'arrays.h'
!         integration time is:
   common/txx/xk
   real*8 par(5)
   common/pxx/par
   character*12 minsz, minst, av2, av3
   common/fn0x/minsz, minst, av2, av3
   character*9 alngdr
   common/drxx/alngdr
   character*21 tname, tname1
!     Counting the bunches and different characteristics
   Bno = zero
!
   dyi = dy(1)
   print *, dyi
   print *, bdef
   if (dyi .le. bdef) then
      print *, "tetets"
      tla = dyi
      tmin = dyi
      tbw = dyi
      Bno = one
      fbd(1) = dyi
      inc = 2
      do Ip = 2, M
         dyip = dy(Ip)
         if (dyip .le. Bdef) then
            inc = inc + 1
            if (dyip .lt. tmin) tmin = dyip
            tbw = tbw + dyip
            tla = dyip
         else
            goto 3029
         end if
      end do
   end if
3029 continue
   if (bno .eq. one) then
      bd(1) = tbw/dfloat(inc - 1)
      bw(1) = tbw
      mindi(1) = tmin
      bsz(1) = inc
      lbd(1) = tla
   end if
   print *, bsz

   do i = 2, M
      dyi = dy(i)
      if (dyi .le. Bdef) then
         if (dy(i - 1) .gt. Bdef) then
            tla = dyi
            tmin = dyi
            tbw = dyi
            Bno = Bno + one
            fbd(Bno) = dyi
            inc = 2
            do Ip = i + 1, M
               dyip = dy(Ip)
               if (dyip .le. Bdef) then
                  inc = inc + 1
                  if (dyip .lt. tmin) tmin = dyip
                  tbw = tbw + dyip
                  tla = dyip
               else
                  goto 3030
               end if
            end do
!
3030        continue
            bd(bno) = tbw/dfloat(inc - 1)
            bw(bno) = tbw
            mindi(bno) = tmin
            bsz(bno) = inc
            lbd(bno) = tla

         end if
      end if
   end do
!                  To account for the situation when the bunch crosses the
!                   boundary conditions:

   if (dy(1) .le. Bdef .and. dy(M) .le. Bdef) then
      ibw = bsz(1) + bsz(bno) - 1
      bw(1) = bw(1) + bw(bno)
      bd(1) = bw(1)/dfloat(ibw - 1)
      mindi(1) = Dmin1(mindi(1), mindi(bno))
      bsz(1) = ibw
      fbd(1) = fbd(bno)
      bsz(bno) = 0
      bd(bno) = 0.0d0
      bw(bno) = 0.0d0
      mindi(bno) = 0.0d0
      lbd(bno) = 0.0d0
      bno = bno - 1
   end if
!****

!****
   print *, bsz
   if (bno .gt. 0) then
      stime = dfloat(M)/dfloat(bno)
      do i = 1, bno
         bszi = bsz(i)
         if (bszi .gt. 300) then
            write (*, *) 'Alert:Bunch Size exceeds 300'
            goto 1035
         end if
         tmin = mindi(i)
         avmn(bszi) = avmn(bszi) + tmin
         if (tmin .lt. absmin(bszi)) absmin(bszi) = tmin
         tmin = fbd(i)
         avf(bszi) = avf(bszi) + tmin
         if (tmin .lt. fdmin(bszi)) fdmin(bszi) = tmin
         tmin = bd(i)
         avbd(bszi) = avbd(bszi) + tmin
         if (tmin .lt. bdmin(bszi)) bdmin(bszi) = tmin
         tmin = bw(i)
         avbw(bszi, 1) = avbw(bszi, 1) + tmin
         avbw(bszi, 2) = avbw(bszi, 2) + one
         if (tmin .lt. bwmin(bszi)) bwmin(bszi) = tmin
         avl(bszi) = avl(bszi) + lbd(i)
1035     continue
      end do
!***
      ! print *, "LOOP 1:"
      ! do i = 2, 300
      !    dyi = absmin(i)
      !    if (dyi .le. bdef) print *, i, dyi, fdmin(i), bwmin(i), bdmin(i)
      ! end do
!*****
      print *, "LOOP 2:"
      do i = 2, 300
         dyi = avbw(i, 2)
         if (dyi .gt. zero) then
            print *, "R1: ", i, avbw(i, 1)/dyi, avbd(i)/dyi
            print *, "R2: ", i, avmn(i)/dyi, avf(i)/dyi, avl(i)/dyi
         end if
      end do
      close (10)

   end if
!*******************************************************
!

end subroutine h_stat
