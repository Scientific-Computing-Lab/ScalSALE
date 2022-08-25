module geometry_module
   use general_utils_module, only : Deter
   implicit none
contains


   function Quad_area(x1, y1, x2, y2, x3, y3, x4, y4) result (area)
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4
      real(8)             :: area

      area = 0.5d0 * ((x3 - x1) * (y4 - y2) - (x2 - x4) * (y1 - y3))

   end function Quad_area

   function Quad_volume(x1, y1, x2, y2, x3, y3, x4, y4, cyl) result (quad_vol)
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, x4, y4, cyl
      real(8)             :: quad_vol

      quad_vol = Triangle_volume(x1, y1, x2, y2, x3, y3, cyl) + Triangle_volume(x3, y3, x4, y4, x1, y1, cyl)

   end function Quad_volume

   function Triangle_area(x1, y1, x2, y2, x3, y3) result (area_tri)
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3
      real(8)             :: area_tri

      area_tri = (x3 - x2) * (y1 - y2) - (x1 - x2) * (y3 - y2)
   end function Triangle_area

   function Triangle_volume(x1, y1, x2, y2, x3, y3, cyl) result (tri_vol)
      real(8), intent(in) :: x1, y1, x2, y2, x3, y3, cyl
      real(8)             :: tri_vol
      real(8) :: tmp
      if (cyl == 0d0) then
         tri_vol = (3d0 * (((x3 - x2) * (y1 - y2)) - ((x1 - x2) * (y3 - y2)))) / 6d0
      else
         tri_vol = (x1 + x2 + x3) * (((x3 - x2) * (y1 - y2)) - ((x1 - x2) * (y3 - y2))) / 6d0
      end if
   end function Triangle_volume

   subroutine Physical_tetraeder_volume_3d(a, b, c, i1, i2, i4, ic1, ic4, t_vol)
      implicit none
      real(8), intent(in)  :: a, b, c
      integer, intent(in)  :: i1, i2, i4, ic1, ic4
      real(8), intent(out) :: t_vol
      real(8)              :: xs, ys, zs, x_cut, y_cut, z_cut, x_cen, y_cen, z_cen, t_vol1, t_vol2, xc, yc, zc
      integer              :: in_out, in_tot

      save /cut_r/
      common /cut_r/ xs(8), ys(8), zs(8), x_cut(12), y_cut(12), z_cut(12), x_cen, y_cen, z_cen

      save /cut_i/
      common /cut_i/ in_out(8) 
      in_tot = in_out(i1) + in_out(i2) + in_out(i4)
      if (in_tot == 0) then
         t_vol = 0.
         return
      end if

      if (in_tot == 3) then
         t_vol = Tetrahederon_volume (x_cen, y_cen, z_cen, xs(i1), ys(i1), zs(i1), xs(i2), ys(i2), zs(i2), xs(i4), ys(i4), zs(i4))
         return
      end if

      if (in_out(i1) == 1) then
         if (in_out(i2) == 1) then
            t_vol1 = Tetrahederon_volume (x_cen, y_cen, z_cen, xs(i1), ys(i1), zs(i1), xs(i2), ys(i2), zs(i2), xs(i4)&
                        , ys(i4), zs(i4))

            call Cut_plane (a, b, -1d0, -c, xs(i2), ys(i2), zs(i2), xs(i4), ys(i4), zs(i4), xc, yc, zc)

            t_vol2 = Tetrahederon_volume (x_cen,y_cen, z_cen, x_cut(ic4), y_cut(ic4), z_cut(ic4), xc, yc, zc, xs(i4), ys(i4)&
            , zs(i4))
            t_vol = t_vol1 - t_vol2
            return
         else if (in_out(i4) == 1) then
            t_vol1 = Tetrahederon_volume (x_cen,y_cen, z_cen, xs(i1), ys(i1), zs(i1), xs(i2), ys(i2), zs(i2), xs(i4), ys(i4)&
            , zs(i4))

            call Cut_plane (a, b, -1d0, -c, xs(i2), ys(i2), zs(i2), xs(i4), ys(i4), zs(i4), xc, yc, zc)

            t_vol2 = Tetrahederon_volume (x_cen,y_cen, z_cen, x_cut(ic1), y_cut(ic1), z_cut(ic1), xs(i2),ys(i2),zs(i2), xc, yc, zc)

            t_vol = t_vol1 - t_vol2
            return
         else
            t_vol = Tetrahederon_volume (x_cen,y_cen, z_cen, xs(i1), ys(i1), zs(i1), x_cut(ic1), y_cut(ic1), z_cut(ic1), &
                                         x_cut(ic4), y_cut(ic4), z_cut(ic4))
            return
         end if
      end if
      if (in_out(i2) == 1) then
         if (in_out(i4) == 1) then
            t_vol1 = Tetrahederon_volume (x_cen,y_cen, z_cen, xs(i1), ys(i1), zs(i1), xs(i2), ys(i2), zs(i2), &
                                          xs(i4), ys(i4), zs(i4))

            t_vol2 = Tetrahederon_volume (x_cen,y_cen, z_cen, xs(i1), ys(i1), zs(i1), x_cut(ic1), y_cut(ic1),&
                                          z_cut(ic1), x_cut(ic4), y_cut(ic4), z_cut(ic4))

            t_vol = t_vol1 - t_vol2
            return
         else
            call Cut_plane(a, b, -1d0, -c, xs(i2), ys(i2), zs(i2), xs(i4), ys(i4), zs(i4), xc, yc, zc)
            t_vol = Tetrahederon_volume (x_cen,y_cen, z_cen, x_cut(ic1), y_cut(ic1), z_cut(ic1), xs(i2), ys(i2), zs(i2), xc, yc, zc)
            return
         end if
      end if
      call Cut_plane (a, b, -1d0, -c, xs(i2), ys(i2), zs(i2), xs(i4), ys(i4), zs(i4), xc, yc, zc)
      t_vol = Tetrahederon_volume (x_cen, y_cen, z_cen, xs(i4), ys(i4), zs(i4), x_cut(ic4), y_cut(ic4), z_cut(ic4), xc, yc, zc)

      return
   end subroutine Physical_tetraeder_volume_3d



   subroutine Cut_point_1(a, b, x1, y1, x2, y2, x_cut, y_cut)
      real(8), intent(in) :: a, b, x1, y1, x2, y2
      real(8), intent(out) :: x_cut, y_cut
      real(8) :: dist_eps  
      dist_eps = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2) * 1d-30
      x_cut = (y1 * x2 - x1 * y2 + b * (x1 - x2)) / (a * (x2 - x1) - y2 + y1 + dist_eps)
      y_cut = a * x_cut + b
   end subroutine Cut_point_1

   subroutine Cut_point_2(a, b, c, xl1, yl1, xl2, yl2, x_cut, y_cut)
      real(8), intent(inout) :: a, b, c, xl1, yl1, xl2, yl2, x_cut, y_cut
   end subroutine Cut_point_2

   subroutine Elipse_cut_line(y_rad, x_rad, y_power, x_power, yc, x_c, y1, x1, y2 , x2, yh, xh)
      implicit none
      real(8), intent(in)     :: y_rad, x_rad, y_power, x_power, yc, x_c
      real(8), intent(in out) :: y1, x1, y2 , x2                                 
      real(8), intent(out)    :: yh, xh                                          

      integer :: n, nfirst, i
      real(8) :: x_temp, y_temp, y_elipse, x_elipse, x1_temp, y1_temp, x2_temp, y2_temp
      real(8) :: dble_fh1, abs_yelipse
      dble_fh1 = dble(501 - 1)


      if (x1 == 0d0 .and. x2 == 0d0) then
         xh = 0d0
         if (y2 > y1) then
           yh = yc + y_rad
         else
           yh = yc - y_rad
         end if
         return
      end if

      do n = 1, 30
        nfirst = 0
        do i = 1, 501
          x_temp = x1 + (x2 - x1) * dble(i - 1) / dble_fh1
          y_temp = y1 + (y2 - y1) * dble(i - 1) / dble_fh1
          y_elipse = y_temp
          abs_yelipse = abs(y_elipse - yc)
          if (abs_yelipse > y_rad) then
            x_elipse = 1d10
          else
            x_elipse = x_rad * (1d0 - (abs_yelipse / y_rad) ** y_power) ** (1d0 / x_power) + x_c
          end if
          if ((x_temp > x_elipse .or. abs_yelipse > y_rad) .and. nfirst == 0) then
            x2_temp = x_temp
            y2_temp = y_temp
            x1_temp = x1 + (x2 - x1) * dble((i - 1) -1) / dble_fh1
            y1_temp = y1 + (y2 - y1) * dble((i - 1) -1) / dble_fh1
            nfirst  = 1

          end if
        end do
        if ( nfirst == 0) then
        end if
        x1 = x1_temp
        y1 = y1_temp
        x2 = x2_temp
        y2 = y2_temp
      end do
      xh = (x1 + x2) / 2d0
      yh = (y1 + y2) / 2d0
      return

   end subroutine Elipse_cut_line

   function Area(xl1, yl1, xl2, yl2, xl3, yl3, xl4, yl4) result (area_)
      real(8), intent(in) :: xl1, yl1, xl2, yl2, xl3, yl3, xl4,yl4
      real(8)             :: area_
   end function Area

   function Volume_fraction1(a, b, sd, geo_vol, x1, y1, x2, y2, x3, y3, x4, y4, cyl) result(phys_vol)
      real(8), intent(in) :: a, b, sd, geo_vol, x1, y1, x2, y2, x3, y3, x4, y4, cyl
      real(8)             :: x(4), y(4), yt, physv, xl1, xl2, yl1, yl2
      integer             :: nx1, nx2, nx3, n1, n2, i, next, last
      real(8)             :: phys_vol

      x(1) = x1
      x(2) = x2
      x(3) = x3
      x(4) = x4
      y(1) = y1
      y(2) = y2
      y(3) = y3
      y(4) = y4
      nx1 = 0
      nx2 = 0
      n1 = 0
      n2 = 0

      do i = 1, 4
         yt = a * x(i) + b
         if (sign(1.d00, y(i) - yt) /= sd ) then
            n1 = n1 + 1
            nx1 = i
         else
            n2 = n2 + 1
            if (nx2 /= 0) nx3 = nx2
            nx2 = i
         end if
      end do

      if (n1 == 4) then
         phys_vol = 0d0
         return
      end if

      if (n2 == 4) then 
         phys_vol = geo_vol
         return
      end if
      if (n1 /= 1) then 
         if (n1 /= 2) then 
            next = nx2 + 1
            if (next == 5) next = 1
            last = nx2 - 1
            if (last == 0) last = 4
            call Cut_point_1(a, b, x(nx2), y(nx2), x(next), y(next), xl1, yl1)
            call Cut_point_1(a, b, x(nx2), y(nx2), x(last), y(last), xl2, yl2)
            physv = Triangle_volume(x(nx2), y(nx2), xl1, yl1, xl2, yl2, cyl)/2.d0
            phys_vol = physv
            return
         end if
      else
         next = nx1 + 1
         if (next == 5) next = 1
         last = nx1 - 1
         if (last == 0) last = 4
         call Cut_point_1(a, b, x(nx1), y(nx1), x(next), y(next), xl1, yl1)
         call Cut_point_1(a, b, x(nx1), y(nx1), x(last), y(last), xl2, yl2)
         physv = geo_vol - Triangle_volume(x(nx1), y(nx1), xl1, yl1, xl2, yl2, cyl) / 2.d0
         phys_vol = physv
         return
      end if
      if ((nx3 == 1).and.(nx2 == 4)) then 
         nx1 = nx3
         nx3 = nx2
         nx2 = nx1
      end if
      last = nx3 - 1
      if (last == 0) last = 4
      next = nx2 + 1
      if (next == 5) next = 1
      call Cut_point_1(a, b, x(nx2), y(nx2), x(next), y(next), xl2, yl2)
      call Cut_point_1(a, b, x(nx3), y(nx3), x(last), y(last), xl1, yl1)
      physv = Quad_volume(xl1, yl1, x(nx3), y(nx3), x(nx2), y(nx2), xl2, yl2, cyl) / 2.0d0
      phys_vol = physv
   end function Volume_fraction1


   function Volume_fraction2(a, b, c, geo_vol, x1, y1, x2, y2, x3, y3, x4, y4, cyl) result(phys_vol)
      real(8), intent(in) :: a, b, c, geo_vol, x1, y1, x2, y2, x3, y3, x4, y4, cyl
      real(8)             :: x(4), y(4), ct, physv, xl1, xl2, yl1, yl2
      real(8)             :: phys_vol
      integer             :: nx1, nx2, nx3, n1, n2, i, next, last


      x(1) = x1
      y(1) = y1
      x(2) = x2
      y(2) = y2
      x(3) = x3
      y(3) = y3
      x(4) = x4
      y(4) = y4
      nx1 = 0
      nx2 = 0
      n1 = 0
      n2 = 0
      do i = 1, 4
         ct = -a * x(i) - b * y(i)
         if (ct < c) then
            n2 = n2 + 1
            if (nx2 /= 0) nx3 = nx2
            nx2 = i
         else
            n1 = n1 + 1
            nx1 = i
         end if
      end do
      if (n1 == 4) then
         phys_vol = 0.d0
         return
      end if
      if (n2 == 4) then
         phys_vol = geo_vol
         return
      end if

      if (n1 == 1) then
         next = nx1 + 1
         if (next == 5) next = 1
         last = nx1 - 1
         if (last == 0) last = 4
         call Cut_line_different(a, b, c, x(nx1), y(nx1), x(next), y(next), xl1, yl1)
         call Cut_line_different(a, b, c, x(nx1), y(nx1), x(last), y(last), xl2, yl2)
         phys_vol = geo_vol - Triangle_volume(x(nx1), y(nx1), xl1, yl1, xl2, yl2, cyl) / 2.d0
         return
      else if (n1 == 2) then
         if ((nx3 == 1) .and. (nx2 == 4)) then
            nx1 = nx3
            nx3 = nx2
            nx2 = nx1
         end if
         last = nx3 - 1
         if (last == 0) last = 4
         next = nx2 + 1
         if (next == 5) next = 1
         call Cut_line_different(a, b, c, x(nx2), y(nx2), x(next), y(next), xl2, yl2)
         call Cut_line_different(a, b, c, x(nx3), y(nx3), x(last), y(last), xl1, yl1)
         phys_vol = Quad_volume(xl1, yl1, x(nx3), y(nx3), x(nx2), y(nx2), xl2, yl2, cyl) / 2.0d0
         return
      else
         next = nx2 + 1
         if (next == 5) next = 1
         last = nx2 - 1
         if (last == 0) last = 4
         call Cut_line_different(a, b, c, x(nx2), y(nx2), x(next), y(next), xl1, yl1)
         call Cut_line_different(a, b, c, x(nx2), y(nx2), x(last), y(last), xl2, yl2)
         phys_vol = Triangle_volume(x(nx2), y(nx2), xl1, yl1, xl2, yl2, cyl) / 2.d0
         return
      end if

   end function Volume_fraction2


   function Volume_fraction_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5 ,x6, y6, z6,&
                               x7, y7, z7, x8, y8, z8, a, b, c, side, geo_vol) result(phys_vol)
      implicit none
      real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5 ,x6, y6, z6,&
                             x7, y7, z7, x8, y8, z8, a, b, c, side, geo_vol

      real(8)     :: zt,x_tot, y_tot, z_tot, tot_vol, temp, xs, ys, zs, x_cut, y_cut, z_cut, x_cen, y_cen, z_cen
      integer     :: i, in_tot, n_cut, in_out
      real(8)     :: phys_vol

      save /cut_r/
      common /cut_r/ xs(8), ys(8), zs(8), x_cut(12), y_cut(12), z_cut(12), x_cen, y_cen, z_cen
      save /cut_i/
      common /cut_i/ in_out(8)

      xs(1) = x1
      ys(1) = y1
      zs(1) = z1
      xs(2) = x2
      ys(2) = y2
      zs(2) = z2
      xs(3) = x3
      ys(3) = y3
      zs(3) = z3
      xs(4) = x4
      ys(4) = y4
      zs(4) = z4
      xs(5) = x5
      ys(5) = y5
      zs(5) = z5
      xs(6) = x6
      ys(6) = y6
      zs(6) = z6
      xs(7) = x7
      ys(7) = y7
      zs(7) = z7
      xs(8) = x8
      ys(8) = y8
      zs(8) = z8
      x_cut(1 : 12) = 1.2345d21
      in_tot = 0
      do i = 1, 8
         zt = a * xs(i) + b * ys(i) + c
         if (sign(1.d00, zs(i) - zt) == side) then
            in_out(i) = 1
         else
            in_out(i) = 0
         end if
         in_tot = in_tot + in_out(i)
      end do
      if (in_tot == 0) then
         phys_vol = 0.
         return
      end if
      if (in_tot == 8) then
         phys_vol = geo_vol
         return
      end if
      n_cut=0
      if (in_out(1) /= in_out(2)) then
         call Cut_plane( a, b, -1d0, -c, xs(1), ys(1), zs(1), xs(2), ys(2), zs(2), &
                     x_cut(1), y_cut(1), z_cut(1))
         n_cut = n_cut + 1
      end if
      if (in_out(2) /= in_out(3)) then
         call Cut_plane( a, b, -1d0, -c, xs(2), ys(2), zs(2), xs(3), ys(3), zs(3), &
                     x_cut(2), y_cut(2), z_cut(2))
         n_cut = n_cut + 1
      end if
      if (in_out(3)/=in_out(4)) then
         call Cut_plane( a, b, -1d0, -c, xs(3), ys(3), zs(3), xs(4), ys(4), zs(4), &
                     x_cut(3), y_cut(3), z_cut(3))
         n_cut = n_cut + 1
      end if
      if (in_out(4) /= in_out(1)) then
         call Cut_plane( a, b, -1d0, -c, xs(4), ys(4), zs(4), xs(1), ys(1), zs(1), &
                     x_cut(4), y_cut(4), z_cut(4))
         n_cut = n_cut + 1
      end if
      if (in_out(5) /= in_out(6)) then
         call Cut_plane( a, b, -1d0, -c, xs(5), ys(5), zs(5), xs(6), ys(6), zs(6), &
                     x_cut(5), y_cut(5), z_cut(5))
         n_cut = n_cut + 1
      end if
      if (in_out(6) /= in_out(7)) then
         call Cut_plane( a, b, -1d0, -c, xs(6), ys(6), zs(6), xs(7), ys(7), zs(7), &
                     x_cut(6), y_cut(6), z_cut(6))
         n_cut = n_cut + 1
      end if
      if (in_out(7) /= in_out(8)) then
         call Cut_plane( a, b, -1d0, -c, xs(7), ys(7), zs(7), xs(8), ys(8), zs(8), &
                     x_cut(7), y_cut(7), z_cut(7))
         n_cut = n_cut + 1
      end if
      if (in_out(8) /= in_out(5)) then
         call Cut_plane( a, b, -1d0, -c, xs(8), ys(8), zs(8), xs(5), ys(5), zs(5), &
                     x_cut(8), y_cut(8), z_cut(8))
         n_cut = n_cut + 1
      end if
      if (in_out(1) /= in_out(5)) then
         call Cut_plane( a, b, -1d0, -c, xs(1), ys(1), zs(1), xs(5), ys(5), zs(5), &
                     x_cut(9), y_cut(9), z_cut(9))
         n_cut = n_cut + 1
      end if
      if (in_out(2) /= in_out(6)) then
         call Cut_plane( a, b, -1d0, -c, xs(2), ys(2), zs(2), xs(6), ys(6), zs(6), &
                     x_cut(10), y_cut(10), z_cut(10))
         n_cut = n_cut + 1
      end if
      if (in_out(3) /= in_out(7)) then
         call Cut_plane( a, b, -1d0, -c, xs(3), ys(3), zs(3), xs(7), ys(7), zs(7), &
                     x_cut(11), y_cut(11), z_cut(11))
         n_cut = n_cut + 1
      end if
      if (in_out(4) /= in_out(8)) then
         call Cut_plane( a, b, -1d0, -c, xs(4), ys(4), zs(4), xs(8), ys(8), zs(8), &
                     x_cut(12), y_cut(12), z_cut(12))
         n_cut = n_cut + 1
      end if
      x_tot = 0.
      y_tot = 0.
      z_tot = 0.
      do i=1,12
         if (x_cut(i) /= 1.2345d21) then
            x_tot = x_tot + x_cut(i)
            y_tot = y_tot + y_cut(i)
            z_tot = z_tot + z_cut(i)
         end if
      end do
      x_cen = x_tot / n_cut
      y_cen = y_tot / n_cut
      z_cen = z_tot / n_cut
      tot_vol = 0.
      call Square_pyramid_volume(a, b, c, 1, 2, 3, 4, 1, 2, 3, 4, temp)
      tot_vol = tot_vol + temp
      call Square_pyramid_volume(a, b, c, 8, 7, 6 ,5 , 7, 6, 5, 8, temp)
      tot_vol = tot_vol + temp
      call Square_pyramid_volume(a, b, c, 1, 5, 6, 2, 9, 5, 10, 1, temp)
      tot_vol = tot_vol + temp
      call Square_pyramid_volume(a, b, c, 2, 6, 7, 3, 10 ,6, 11, 2, temp)
      tot_vol = tot_vol + temp
      call Square_pyramid_volume(a, b, c, 3, 7, 8, 4, 11, 7, 12, 3, temp)
      tot_vol = tot_vol + temp
      call Square_pyramid_volume(a, b, c, 4, 8, 5, 1, 12, 8, 9, 4, temp)
      tot_vol = tot_vol + temp
      phys_vol = tot_vol

      return
   end function Volume_fraction_3d


   subroutine Square_pyramid_volume(a, b, c, i1, i2, i3, i4, ic1, ic2, ic3, ic4, pvol)
      implicit none
      real(8), intent(in)  :: a, b, c
      integer, intent(in)  :: i1, i2, i3, i4, ic1, ic2, ic3, ic4
      real(8), intent(out) :: pvol

      real(8)              :: xs, ys, zs, x_cut, y_cut, z_cut, x_cen, y_cen, z_cen, temp, tri_vol
      integer :: in_out, in_tot
      save /cut_r/
      common /cut_r/ xs(8), ys(8), zs(8), x_cut(12), y_cut(12), z_cut(12), x_cen, y_cen, z_cen
      save /cut_i/
      common /cut_i/ in_out(8)

      in_tot = in_out(i1) + in_out(i2) + in_out(i3) + in_out(i4)
      if (in_tot == 0) then
         pvol = 0.
         return
      end if
      temp = 0.
      call Physical_tetraeder_volume_3d (a, b, c, i1, i2, i4, ic1, ic4, tri_vol)
      temp = temp + tri_vol
      call Physical_tetraeder_volume_3d (a, b, c, i2, i3, i1, ic2, ic1, tri_vol)
      temp = temp + tri_vol
      call Physical_tetraeder_volume_3d (a, b, c, i3, i4, i2, ic3, ic2, tri_vol)
      temp = temp + tri_vol
      call Physical_tetraeder_volume_3d (a, b, c, i4, i1, i3, ic4, ic3, tri_vol)
      temp = temp + tri_vol
      pvol = temp / 12.
      return
   end subroutine square_pyramid_volume


   subroutine Cut_plane(a, b, c, d, x1, y1, z1, x2, y2, z2, xc, yc, zc)
      implicit none
      real(8), intent(in)     :: a, b, c, d, x1, y1, z1, x2, y2, z2
      real(8), intent(out)    :: xc, yc, zc

      real(8)                 :: x21, y21, z21, det_xy, det_yz, det_zx, deno


      x21 = x2 - x1
      y21 = y2 - y1
      z21 = z2 - z1
      deno = a*x21 + b*y21 + c*z21

      if (abs(deno) < 1.d-12) then
         deno = 1d-12
      end if
      det_xy = x1*y2 - x2*y1
      det_yz = y1*z2 - y2*z1
      det_zx = z1*x2 - z2*x1
      xc = (d*x21 + b*det_xy - c*det_zx) / deno
      yc = (d*y21 - a*det_xy + c*det_yz) / deno
      zc = (d*z21 + a*det_zx - b*det_yz) / deno

   end subroutine Cut_plane

   function Distance_2(xl1 , yl1, xl2, yl2) result (dist_)
      real(8), intent(in) :: xl1, yl1, xl2, yl2
      real(8)             :: dist_

       dist_ = (xl2 - xl1) ** 2 + (yl2 - yl1) ** 2
   end function Distance_2


   subroutine Smooth_line(y1, x1, y2, x2, dx_dy1, dx_dy2, a, b, c, d)
      real(8), intent(in out) :: a, b, c, d
      real(8), intent(in)     :: y1, x1, y2, x2, dx_dy1, dx_dy2
      real(8)                 :: xm11, xm12, xm13, xm14, xm21, xm22, xm23, xm24, xm31, xm32, xm33, xm34, xm41, xm42, xm43, xm44

      if (y1 == y2) then
      end if
      xm11 = (3d0 * y1 - y2) * y2 ** 2 / (y1 - y2) ** 3
      xm12 = (y1 - 3d0 * y2) * y1 ** 2 / (y1 - y2) ** 3
      xm13 = -y1 * y2 ** 2 / (y1 - y2) ** 2
      xm14 = -y1 ** 2 * y2 / (y1 - y2) ** 2

      xm21 = -6d0 * y1 * y2 / (y1 - y2) ** 3
      xm22 =  6d0 * y1 * y2 / (y1 - y2) ** 3
      xm23 = y2 * (2d0 * y1 + y2) / (y1 - y2) ** 2
      xm24 = y1 * (y1 + 2d0 * y2) / (y1 - y2) ** 2

      xm31 =  3d0 * (y1 + y2) / (y1 - y2) ** 3
      xm32 = -3d0 * (y1 + y2) / (y1 - y2) ** 3
      xm33 = -(y1 + 2d0 * y2) / (y1 - y2) ** 2
      xm34 = -(2d0 * y1 + y2) / (y1 - y2) ** 2

      xm41 = -2d0 / (y1 - y2) ** 3
      xm42 =  2d0 / (y1 - y2) ** 3
      xm43 =  1d0 / (y1 - y2) ** 2
      xm44 =  1d0 / (y1 - y2) ** 2

      a = x1 * xm11 + x2 * xm12 + dx_dy1 * xm13 + dx_dy2 * xm14
      b = x1 * xm21 + x2 * xm22 + dx_dy1 * xm23 + dx_dy2 * xm24
      c = x1 * xm31 + x2 * xm32 + dx_dy1 * xm33 + dx_dy2 * xm34
      d = x1 * xm41 + x2 * xm42 + dx_dy1 * xm43 + dx_dy2 * xm44

      return
   end subroutine Smooth_line

   subroutine Cut_point_smooth_line(a, b, c, d, y1, x1 ,y2 ,x2 ,yh, xh)
      real(8), intent(in)     :: a, b, c, d
      real(8), intent(in out) :: x1, x2, y1,y2, xh, yh
      real(8)                 :: x_temp, y_temp, x1_t, x2_t, y1_t, y2_t, x_smooth, y_smooth
      integer                 :: n, n_first, i

      do n = 1, 20
         n_first = 0
         do i = 1, 501
            x_temp = x1 + (x2 - x1) * dble(i - 1) / dble(501 - 1)
            y_temp = y1 + (y2 - y1) * dble(i - 1) / dble(501 - 1)
            y_smooth = y_temp
            x_smooth = a + b * y_smooth + c * y_smooth ** 2 + d * y_smooth ** 3
            if (x_temp > x_smooth .and. n_first == 0) then
               x2_t = x_temp
               y2_t = y_temp
               x1_t = x1 + (x2 - x1) * dble((i - 1) - 1) / dble(501 - 1)
               y1_t = Y1 + (y2 - y1) * dble((i - 1) - 1) / dble(501 - 1)
               n_first = 1
            end if
        end do 
        if (n_first == 0) then
        end if
        x1 = x1_t
        y1 = y1_t
        x2 = x2_t
        y2 = y2_t
      end do 
      xh = (x1 + x2) / 2d0
      yh = (y1 + y2) / 2d0
      return
   end subroutine Cut_point_smooth_line

   subroutine Cut_line(x11, y11, x12, y12, x21, y21, x22, y22, yn_cut, xl_cut, yl_cut, ins)
       real(8), intent(in)  :: x11, y11, x12, y12, x21, y21, x22, y22
       integer, intent(in)  :: ins
       integer, intent(out) :: yn_cut
       real(8), intent(out) :: xl_cut, yl_cut
       real(8)              :: do_fiani2, dmf2, dx1, dy1, a1, b1, c1, dx2, dy2, a2, b2, c2, deter, eps





       do_fiani2 = min(Distance_2(x21, y21, x22, y22), Distance_2(x11 ,y11 ,x12 ,y12))
       dmf2 = 1d-12 * do_fiani2

       dx1 = x12 - x11
       dy1 = y12 - y11
       if (abs(dx1) < abs(dy1)) then
           a1 = 1d0
           call Line_parameters(x11, y11, x12, y12, b1, c1)
       else
           b1 = 1d0
           call Line_parameters(y11, x11, y12, x12, a1, c1)
       end if
       dx2 = x22 - x21
       dy2 = y22 - y21
       if (abs(dx2) < abs(dy2)) then
           a2 = 1d0
           call Line_parameters(x21, y21, x22, y22, b2, c2)
       else
           b2 = 1d0
           call Line_parameters(y21, x21, y22, x22, a2, c2)
       end if
       deter = b2 * a1 - a2 * b1
       if (max(abs(a1 * b2), abs(a2 * b1)) < 1d-20) then
           yn_cut = 0
           return
       end if

       if (abs(deter) / max(abs(a1 * b2), abs(a2*b1)) < 1d-8) then
           yn_cut=0
           return
       end if

       xl_cut =  (b1 * c2 - b2 * c1) / deter
       yl_cut = -(a1 * c2 - a2 * c1) / deter
       yn_cut = 1
       if (ins == 2) return


       eps = dmf2 / 1d2
       if (max((xl_cut - x11) * (xl_cut - x12), (yl_cut - y11) * (yl_cut - y12)) > eps) then
           yn_cut = 0
           return
       end if
       if (max((xl_cut - x21) * (xl_cut - x22), (yl_cut - y21) * (yl_cut - y22)) > eps) then
           yn_cut = 0
           return
       end if
       if (ins == 1) return


       if (((xl_cut - x11) ** 2 + (yl_cut - y11) ** 2) < dmf2.or.&
           ((xl_cut - x12) ** 2 + (yl_cut - y12) ** 2) < dmf2.or.&
           ((xl_cut - x21) ** 2 + (yl_cut - y21) ** 2) < dmf2.or.&
           ((xl_cut - x22) ** 2 + (yl_cut - y22) ** 2) < dmf2) then
           yn_cut = 0
           return
       end if
       return

   end subroutine Cut_line


   subroutine Cut_line_different (a, b, c, xl1, yl1, xl2, yl2, x_cut, y_cut)
      real(8), intent(in)  :: a, b, c, xl1, xl2, yl1, yl2
      real(8), intent(out) :: x_cut, y_cut
      real(8)              :: al, bl, cl, det, xleps

      al = yl2 - yl1
      bl = xl1 - xl2
      cl = (yl1 - yl2) * xl1 + (xl2 - xl1) * yl1
      xleps = sqrt((xl1 - xl2) ** 2d0 + (yl1 - yl2) ** 2d0) * 1d-30
      det = al * b - a * bl
      if (abs(det) < 1d-30) then
        x_cut = 1d30
        y_cut = 1d30
        return
      end if
      x_cut = (-cl * b + c * bl) / det
      y_cut = ( cl * a - c * al) / det
      return
   end subroutine Cut_line_different


   subroutine Rectangle_cut_line(x_rect, y_rect, xl1, yl1, xl2, yl2, is_cut, way_in)
      real(8), intent(in out) :: x_rect(4), y_rect(4), xl1, yl1, xl2, yl2
      real(8), intent(in out) :: way_in
      integer, intent(in out) :: is_cut
   end subroutine Rectangle_cut_line

   subroutine Line_parameters(x1, y1, x2, y2, b, c)
      real(8), intent(in)     :: x1, y1, x2, y2
      real(8), intent(in out) :: b, c


      if (y1 == y2) then
      end if
      b = (x2 - x1) / (y1 - y2)
      c = -(x1 + b * y1)
   end subroutine Line_parameters

   subroutine Cut_pnt(x1, y1, x2, y2, x3, y3, xc, yc ,rs ,xn ,yn)
      real(8), intent(in out) ::x1, y1, x2, y2, x3, y3, xc, yc ,rs ,xn ,yn
   end subroutine Cut_pnt

   subroutine Cut_pntn(xc1, yc1, xc2, yc2, xc3, yc3, xc, yc, rs, xn, yn)
      real(8), intent(in out) ::xc1, yc1, xc2, yc2, xc3, yc3, xc, yc, rs, xn, yn
   end subroutine Cut_pntn

   subroutine Circus(x1, y1, x2, y2, x3 ,y3 , xc, yc, rs)
      implicit real*8 (A-H,O-Z)
   end subroutine Circus

   subroutine Mirror_image(x1, y1, x2, y2, x3, y3, x_mirror, y_mirror)
      real(8), intent(out) :: x_mirror, y_mirror
      real(8), intent(in ) :: x1, y1, x2, y2, x3, y3

      real(8) :: dx, dy, dx_sq, dy_sq, dx_dy, dist, xy12, xm, ym

      dx    = x2 - x1
      dy    = y2 - y1
      dx_sq = dx * dx
      dy_sq = dy * dy
      dx_dy = dx * dy
      dist  = dx_sq + dy_sq

      xy12 = x1 * y2 - x2 * y1
      xm = (dx_sq * x3 + dx_dy * y3 + dy * xy12) / dist
      ym = (dy_sq * y3 + dx_dy * x3 - dx * xy12) / dist
      x_mirror = 2d0 * xm - x3
      y_mirror = 2d0 * ym - y3
      return
   end subroutine Mirror_image

   subroutine Point_distance_line(x1, y1, x2, y2, x0, y0, xc, yc, xl)
      real(8), intent(in out) :: x1, y1, x2, y2, x0, y0, xc, yc, xl
   end subroutine Point_distance_line

   subroutine Area_fraction_side1(a, b, side, xp, yp, xm, ym, area_fac_len, len_fact, cyl, om_cyl)
      real(8), intent(in) :: a, b, side, xp, yp, xm, ym, area_fac_len, len_fact, cyl, om_cyl
   end subroutine Area_fraction_side1

   subroutine Area_fraction_side2(a, b, c, xp, yp, xm, ym, area_fac_len, len_fact, cyl, om_cyl)
      real(8), intent(in) :: a, b, c, xp, yp, xm, ym, area_fac_len, len_fact, cyl, om_cyl
   end subroutine Area_fraction_side2

   subroutine Polygon_volume(n_vert, x_vert, y_vert, vol_pol, cyl)
      implicit none
      real(8), dimension(:),  intent(in) :: x_vert  
      real(8), dimension(:),  intent(in) :: y_vert  
      integer, intent(in)  :: n_vert     
      real(8), intent(in)  :: cyl        
      real(8), intent(out) :: vol_pol    
      real(8)              :: omcyl      
      real(8)              :: x1, x2, x3 
      real(8)              :: y1, y2, y3 
      integer              :: i          


      vol_pol = 0d0
      x1 = x_vert(1)
      y1 = y_vert(1)
      do i = 2, n_vert - 1
        x2 = x_vert(i)
        y2 = y_vert(i)
        x3 = x_vert(i + 1)
        y3 = y_vert(i + 1)
        vol_pol = vol_pol + Triangle_volume(x1, y1, x2, y2, x3, y3, cyl)
      end do

   end subroutine Polygon_volume

   subroutine Find_cut_point_line2(start_pnt, dir, r1, z1, r2, z2, cut_z, cut_r, cut_dist, second_dist, dbg, flg)
      real(8), intent(in out) :: start_pnt, dir, r1, z1, r2, z2, cut_z, cut_r, cut_dist, second_dist, dbg, flg
   end subroutine Find_cut_point_line2

   subroutine Two_cell_cut(x1, y1, x2, y2, x_cut, y_cut, n_cut)
      real(8), intent(in)     :: x1(4), y1(4), x2(4), y2(4)
      real(8), intent(in out) :: x_cut, y_cut, n_cut
   end subroutine Two_cell_cut

   function Inside1(xc, yc, x1, y1, dmf2) result (inside)
      real(8)             :: inside
      real(8), intent(in) :: xc, yc, y1, x1, dmf2
   end function Inside1

   function Inside_robust(xc, yc, x1, y1, dmf2) result (inside_r)
      logical             :: inside_r
      real(8), intent(in) :: xc, yc, x1, y1, dmf2
   end function Inside_robust

   subroutine Polygon_line_cut(n_poly, x_pol, y_pol, x1, y1, x2, y2, a, b, side, n_new,x_new, y_new)
      real(8), intent(in) :: n_poly, x_pol, y_pol, x1, y1, x2, y2, a, b, side, n_new,x_new, y_new
   end subroutine Polygon_line_cut

   subroutine Is_cut(y11, x11, y12, x12, y21, x21, y22, x22, i_res)
      real(8), intent(in) ::y11, x11, y12, x12, y21, x21, y22, x22, i_res
   end subroutine Is_cut


   subroutine Rotate_vector(vector, phi, trans_vec)
      real(8), intent(in) :: vector, phi, trans_vec
   end subroutine Rotate_vector

   subroutine Rotate_rank2_tensor(t_xx, t_xy, t_yx, t_yy, t_xx_n, t_xy_n, t_yx_n, t_yy_n, teta)
      real(8), intent (in) :: t_xx, t_xy, t_yx, t_yy, t_xx_n, t_xy_n, t_yx_n, t_yy_n, teta
   end subroutine Rotate_rank2_tensor


   subroutine Average_vof_cell1(x1, y1, x2, y2, x3, y3, x4, y4, a, b, side, xx1, yy1, xx2, yy2, x_av, y_av)
      real(8), intent(in out) :: x1, y1, x2, y2, x3, y3, x4, y4, a, b, side, xx1, yy1, xx2, yy2, x_av, y_av
   end subroutine Average_vof_cell1

   subroutine Average_vof_cell2(x1, y1, x2, y2, x3, y3, x4, y4, a, b, side, xx1, yy1, xx2, yy2, x_av, y_av)
      real(8), intent(in out) ::x1, y1, x2, y2, x3, y3, x4, y4, a, b, side, xx1, yy1, xx2, yy2, x_av, y_av
   end subroutine Average_vof_cell2

   function Root(seed, dx, dd, n) result(root_r)
      implicit none
      real(8)             :: root_r
      real(8), intent(in) :: seed,dx,dd
      integer, intent(in) :: n
      integer             :: n_loop, n_loop_max
      real(8)             :: x_old, err, fx, d_fx, x_new

      n_loop_max = 1000
      x_old = seed
      err = 1.d0
      n_loop = 1
      if (abs(dd / dx - n) < 1d-10) then
         root_r = 1d0
         return
      else
         do while ((err > 1d-15) .and. (n_loop < n_loop_max))
            n_loop = n_loop + 1
            fx     = dx / dd * ((x_old ** n - 1d0) / (x_old - 1d0)) - 1d0
            d_fx   = ((n - 1d0) * x_old ** n - n * x_old ** (n - 1) + 1d0) / ((x_old - 1d0) ** 2) * dx / dd
            x_new  = x_old - fx / d_fx
            err    = abs(x_old - x_new)
            x_old  = x_new
         end do
      end if
      root_r = x_new
   end function Root


   subroutine Layer_cut_line(y11, x11, y12, x12, layer_parm, n_part, yh, xh)
      use constants_module, only : PI, N_SEGMENT_MAX
      implicit none
      integer, intent(in)  :: n_part                
      real(8), intent(in)  :: y11, x11, y12, x12    
      real(8), intent(in)  :: layer_parm(n_part, 4) 
      real(8), intent(out) :: yh, xh                

      real(8) :: x21, y21, x22, y22                  
      real(8) :: x1, x2, y1, y2, a, b, al, be        
      real(8) :: dx_dy1, dx_dy2, x_c, y_c, x_d, y_d  
      real(8) :: tt1, tt2, c, d ,signy               
      real(8) :: xl_cut, yl_cut, x_cut, y_cut        
      integer :: n_cut                               
      integer :: is_cut                              
      integer :: n                                   



      n_cut = -1    
      do n = 1, 1 
          y21 = layer_parm(n, 1)  
          x21 = layer_parm(n, 2)  
          y22 = layer_parm(n, 3)  
          x22 = layer_parm(n, 4)  
          call Cut_line(x11, y11, x12, y12, x21, y21, x22, y22, is_cut, xl_cut, yl_cut, 1) 

          if (is_cut == 1) then 
             n_cut = n
             x_cut = xl_cut
             y_cut = yl_cut
          end if
      end do
      if (n_cut == -1) then
         write(*,*) "DOESNT INTERSECT"
         stop
      end if

         xh = x_cut
         yh = y_cut



   end subroutine Layer_cut_line

   subroutine Layer_cut_elipse(y11, x11, y12, x12, layer_parm,  n_part, yh, xh)
      use constants_module, only : PI, N_SEGMENT_MAX
      implicit none
      integer, intent(in)  :: n_part                
      real(8), intent(in)  :: y11, x11, y12, x12    
      real(8), intent(in)  :: layer_parm(n_part, 11) 
      real(8), intent(out) :: yh, xh                

      real(8) :: x21, y21, x22, y22                  
      real(8) :: x1, x2, y1, y2, a, b, al, be        
      real(8) :: dx_dy1, dx_dy2, xc, yc, x_d, y_d  
      real(8) :: tt1, tt2, c, d ,signy               
      real(8) :: xl_cut, yl_cut, x_cut, y_cut        
      integer :: n_cut                               
      integer :: is_cut                              
      integer :: n                                   



      n_cut = 1    
      call Layer_cut_line(y11,x11,y12,x12,layer_parm, n_part,y_cut, x_cut)


     x1 = 0.90*x_cut+0.1*x11
     y1 = 0.90*y_cut+0.1*y11
     x2 = x12
     y2 = y12
     a  = layer_parm(n_cut, 5)
     b  = layer_parm(n_cut, 6)
     al = layer_parm(n_cut, 7)
     be = layer_parm(n_cut, 8)
     yc = layer_parm(n_cut, 9)
     xc = layer_parm(n_cut, 10)
     call Elipse_cut_line(a, b, al, be, yc, xc, y1, x1, y2, x2,yh,xh)

   end subroutine Layer_cut_elipse




   function Hexahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8)
      implicit none
      real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8
      real(8) hexahedron_volume
      hexahedron_volume =        (-x_proj(y2, z2, y3, z3, y4, z4, y5, z5, y6, z6, y8, z8) * x1 &
                                  -x_proj(y3, z3, y4, z4, y1, z1, y6, z6, y7, z7, y5, z5) * x2 &
                                  -x_proj(y4, z4, y1, z1, y2, z2, y7, z7, y8, z8, y6, z6) * x3 &
                                  -x_proj(y1, z1, y2, z2, y3, z3, y8, z8, y5, z5, y7, z7) * x4 &
                                  +x_proj(y6, z6, y7, z7, y8, z8, y1, z1, y2, z2, y4, z4) * x5 &
                                  +x_proj(y7, z7, y8, z8, y5, z5, y2, z2, y3, z3, y1, z1) * x6 &
                                  +x_proj(y8, z8, y5, z5, y6, z6, y3, z3, y4, z4, y2, z2) * x7 &
                                  +x_proj(y5, z5, y6, z6, y7, z7, y4, z4, y1, z1, y3, z3) * x8 ) / 12d0
      return
  end function Hexahedron_volume


   function x_proj(y2, z2, y3, z3, y4, z4, y5, z5, y6, z6, y8, z8) result (x_area)
      implicit none
      real(8) :: x_area
      real(8), intent(in) :: y2, z2, y3, z3, y4, z4, y5, z5, y6, z6, y8, z8  

      x_area = y2 * z3 + y2 * z4 - y2 * z5 - y2 * z6 - y3 * z2 + y3 * z4 &
             - y4 * z2 - y4 * z3 + y4 * z5 + y4 * z8 + y5 * z2 - y5 * z4 &
             + y5 * z6 - y5 * z8 + y6 * z2 - y6 * z5 - y8 * z4 + y8 * z5

      end function x_proj

   function Tetrahederon_volume(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3) result (vol)
      implicit none
      real(8) :: vol
      real(8), intent(in) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3
      real(8) :: ax, ay, az, bx, by, bz, cx, cy, cz

      ax = x1 - x0
      ay = y1 - y0
      az = z1 - z0
      bx = x2 - x0
      by = y2 - y0
      bz = z2 - z0
      cx = x3 - x0
      cy = y3 - y0
      cz = z3 - z0
      vol = bx * (ay * cz - az * cy) - by * (ax * cz - az * cx) + bz * (ax * cy - ay * cx)
  end function Tetrahederon_volume

   function Vector_grad(u1, u2, u4, u5, y1, z1, y2, z2, y4, z4, y5, z5) result (vec_grad)
      implicit none
      real(8) :: vec_grad
      real(8), intent(in) :: u1, u2, u4, u5  
      real(8), intent(in) :: y1, y2, y4 ,y5  
      real(8), intent(in) :: z1, z2, z4, z5  

     real(8) :: u1u2u4u5, y2y1,y4y1,y5y1,z2z1,z4z1,z5z1
        
      y2y1 = y2 - y1
      y4y1 = y4 - y1
      y5y1 = y5 - y1
      z2z1 = z2 - z1    
      z4z1 = z4 - z1
      z5z1 = z5 - z1
      vec_grad = (u1 + u2 + u4) * ( y4y1 * z2z1 - z4z1 * y2y1 ) + &
                 (u1 + u4 + u5) * ( y5y1 * z4z1 - z5z1 * y4y1 ) + &
                 (u1 + u2 + u5) * ( y2y1 * z5z1 - z2z1 * y5y1 )

   end function Vector_grad


   subroutine Vector_grad_vec(u1, u2, u3, u4, u5, u6, u7, u8, u1u2u4, u1u4u5, u1u2u5, u2u3u1, u2u1u6, u2u3u6, u3u4u2, &
                                                              u3u2u7, u3u4u7, u4u1u3, u4u3u8, u4u1u8, u5u8u6, u5u6u1, &
                                                              u5u8u1, u6u5u7, u6u7u2, u6u5u2, u7u6u8, u7u8u3, u7u6u3, &
                                                              u8u7u5, u8u5u4, u8u7u4)
      implicit none
      real(8), intent(out) :: u1u2u4, u1u4u5, u1u2u5, u2u3u1, u2u1u6, u2u3u6, u3u4u2, &
                              u3u2u7, u3u4u7, u4u1u3, u4u3u8, u4u1u8, u5u8u6, u5u6u1, &
                              u5u8u1, u6u5u7, u6u7u2, u6u5u2, u7u6u8, u7u8u3, u7u6u3, &
                              u8u7u5, u8u5u4, u8u7u4
      real(8), intent(in) :: u1, u2, u3, u4, u5, u6, u7, u8  

      u1u2u4 = u1 + u2 + u4
      u1u4u5 = u1 + u4 + u5
      u1u2u5 = u1 + u2 + u5

      u2u3u1 = u2 + u3 + u1
      u2u1u6 = u2 + u1 + u6
      u2u3u6 = u2 + u3 + u6

      u3u4u2 = u3 + u4 + u2
      u3u2u7 = u3 + u2 + u7
      u3u4u7 = u3 + u4 + u7

      u4u1u3 = u4 + u1 + u3
      u4u3u8 = u4 + u3 + u8
      u4u1u8 = u4 + u1 + u8

      u5u8u6 = u5 + u8 + u6
      u5u6u1 = u5 + u6 + u1
      u5u8u1 = u5 + u8 + u1

      u6u5u7 = u6 + u5 + u7
      u6u7u2 = u6 + u7 + u2
      u6u5u2 = u6 + u5 + u2

      u7u6u8 = u7 + u6 + u8
      u7u8u3 = u7 + u8 + u3
      u7u6u3 = u7 + u6 + u3

      u8u7u5 = u8 + u7 + u5
      u8u5u4 = u8 + u5 + u4
      u8u7u4 = u8 + u7 + u4

   end subroutine Vector_grad_vec


   subroutine Vector_grad_planes(a1, a2, a3, a4, a5, a6, a7, a8, &
                                 b1, b2, b3, b4, b5, b6, b7, b8,  a4a1b2b1b4b1a2a1, a5a1b4b1b5b1a4a1, a2a1b5b1b2b1a5a1, &
                                                                  a1a2b3b2b1b2a3a2, a6a2b1b2b6b2a1a2, a3a2b6b2b3b2a6a2, &
                                                                  a2a3b4b3b2b3a4a3, a7a3b2b3b7b3a2a3, a4a3b7b3b4b3a7a3, &
                                                                  a3a4b1b4b3b4a1a4, a8a4b3b4b8b4a3a4, a1a4b8b4b1b4a8a4, &
                                                                  a6a5b8b5b6b5a8a5, a1a5b6b5b1b5a6a5, a8a5b1b5b8b5a1a5, &
                                                                  a7a6b5b6b7b6a5a6, a2a6b7b6b2b6a7a6, a5a6b2b6b5b6a2a6, &
                                                                  a8a7b6b7b8b7a6a7, a3a7b8b7b3b7a8a7, a6a7b3b7b6b7a3a7, &
                                                                  a5a8b7b8b5b8a7a8, a4a8b5b8b4b8a5a8, a7a8b4b8b7b8a4a8)
      implicit none
      real(8), intent(out) :: a4a1b2b1b4b1a2a1, a5a1b4b1b5b1a4a1, a2a1b5b1b2b1a5a1, &
                              a1a2b3b2b1b2a3a2, a6a2b1b2b6b2a1a2, a3a2b6b2b3b2a6a2, &
                              a2a3b4b3b2b3a4a3, a7a3b2b3b7b3a2a3, a4a3b7b3b4b3a7a3, &
                              a3a4b1b4b3b4a1a4, a8a4b3b4b8b4a3a4, a1a4b8b4b1b4a8a4, &
                              a6a5b8b5b6b5a8a5, a1a5b6b5b1b5a6a5, a8a5b1b5b8b5a1a5, &
                              a7a6b5b6b7b6a5a6, a2a6b7b6b2b6a7a6, a5a6b2b6b5b6a2a6, &
                              a8a7b6b7b8b7a6a7, a3a7b8b7b3b7a8a7, a6a7b3b7b6b7a3a7, &
                              a5a8b7b8b5b8a7a8, a4a8b5b8b4b8a5a8, a7a8b4b8b7b8a4a8
      real(8), intent(in) :: a1, a2, a3, a4, a5, a6, a7, a8  
      real(8), intent(in) :: b1, b2, b3, b4, b5, b6, b7, b8  

     real(8) :: a2a1, a4a1, a5a1, b2b1, b4b1, b5b1, a3a2, a1a2, a6a2, b3b2, b1b2, b6b2, &
                a4a3, a2a3, a7a3, b4b3, b2b3, b7b3, a1a4, a3a4, a8a4, b1b4, b3b4, b8b4, &
                a8a5, a6a5, a1a5, b8b5, b6b5, b1b5, a5a6, a7a6, a2a6, b5b6, b7b6, b2b6, &
                a6a7, a8a7, a3a7, b6b7, b8b7, b3b7, a7a8, a5a8, a4a8, b7b8, b5b8, b4b8

      a2a1 = a2 - a1
      a4a1 = a4 - a1
      a5a1 = a5 - a1
      b2b1 = b2 - b1
      b4b1 = b4 - b1
      b5b1 = b5 - b1

      a4a1b2b1b4b1a2a1 = a4a1 * b2b1 - b4b1 * a2a1
      a5a1b4b1b5b1a4a1 = a5a1 * b4b1 - b5b1 * a4a1
      a2a1b5b1b2b1a5a1 = a2a1 * b5b1 - b2b1 * a5a1

      a3a2 = a3 - a2
      a1a2 = a1 - a2
      a6a2 = a6 - a2
      b3b2 = b3 - b2
      b1b2 = b1 - b2
      b6b2 = b6 - b2

      a1a2b3b2b1b2a3a2 = a1a2 * b3b2 - b1b2 * a3a2
      a6a2b1b2b6b2a1a2 = a6a2 * b1b2 - b6b2 * a1a2
      a3a2b6b2b3b2a6a2 = a3a2 * b6b2 - b3b2 * a6a2

      a4a3 = a4 - a3
      a2a3 = a2 - a3
      a7a3 = a7 - a3
      b4b3 = b4 - b3
      b2b3 = b2 - b3
      b7b3 = b7 - b3

      a2a3b4b3b2b3a4a3 = a2a3 * b4b3 - b2b3 * a4a3
      a7a3b2b3b7b3a2a3 = a7a3 * b2b3 - b7b3 * a2a3
      a4a3b7b3b4b3a7a3 = a4a3 * b7b3 - b4b3 * a7a3

      a1a4 = a1 - a4
      a3a4 = a3 - a4
      a8a4 = a8 - a4
      b1b4 = b1 - b4
      b3b4 = b3 - b4
      b8b4 = b8 - b4

      a3a4b1b4b3b4a1a4 = a3a4 * b1b4 - b3b4 * a1a4
      a8a4b3b4b8b4a3a4 = a8a4 * b3b4 - b8b4 * a3a4
      a1a4b8b4b1b4a8a4 = a1a4 * b8b4 - b1b4 * a8a4

      a8a5 = a8 - a5
      a6a5 = a6 - a5
      a1a5 = a1 - a5
      b8b5 = b8 - b5
      b6b5 = b6 - b5
      b1b5 = b1 - b5

      a6a5b8b5b6b5a8a5 = a6a5 * b8b5 - b6b5 * a8a5
      a1a5b6b5b1b5a6a5 = a1a5 * b6b5 - b1b5 * a6a5
      a8a5b1b5b8b5a1a5 = a8a5 * b1b5 - b8b5 * a1a5

      a5a6 = a5 - a6
      a7a6 = a7 - a6
      a2a6 = a2 - a6
      b5b6 = b5 - b6
      b7b6 = b7 - b6
      b2b6 = b2 - b6

      a7a6b5b6b7b6a5a6 = a7a6 * b5b6 - b7b6 * a5a6
      a2a6b7b6b2b6a7a6 = a2a6 * b7b6 - b2b6 * a7a6
      a5a6b2b6b5b6a2a6 = a5a6 * b2b6 - b5b6 * a2a6

      a6a7 = a6 - a7
      a8a7 = a8 - a7
      a3a7 = a3 - a7
      b6b7 = b6 - b7
      b8b7 = b8 - b7
      b3b7 = b3 - b7

      a8a7b6b7b8b7a6a7 = a8a7 * b6b7 - b8b7 * a6a7
      a3a7b8b7b3b7a8a7 = a3a7 * b8b7 - b3b7 * a8a7
      a6a7b3b7b6b7a3a7 = a6a7 * b3b7 - b6b7 * a3a7

      a7a8 = a7 - a8
      a5a8 = a5 - a8
      a4a8 = a4 - a8
      b7b8 = b7 - b8
      b5b8 = b5 - b8
      b4b8 = b4 - b8

      a5a8b7b8b5b8a7a8 = a5a8 * b7b8 - b5b8 * a7a8
      a4a8b5b8b4b8a5a8 = a4a8 * b5b8 - b4b8 * a5a8
      a7a8b4b8b7b8a4a8 = a7a8 * b4b8 - b7b8 * a4a8

   end subroutine Vector_grad_planes

   function Line_length_in_cell(x1, y1, z1 ,x2, y2, z2 ,x3, y3, z3 ,x4, y4, z4 , &
                                x5, y5, z5 ,x6, y6, z6 ,x7, y7, z7 ,x8, y8, z8 , &
                                 n_x, n_y, n_z, x0, y0, z0) result(length)
      use general_utils_module, only : Deter
      implicit none
      real(8) :: length
      real(8), intent(in) :: x1, y1, z1 ,x2, y2, z2 ,x3, y3, z3 ,x4, y4, z4  
      real(8), intent(in) :: x5, y5, z5 ,x6, y6, z6 ,x7, y7, z7 ,x8, y8, z8  
      real(8), intent(in) :: x0, y0, z0                                      
      real(8), intent(in) :: n_x, n_y, n_z                                   

      real(8), dimension(6) :: len_from_face 
      real(8) :: n_norm                      
      real(8) :: n_xn, n_yn, n_zn            
      real(8) :: len_min_1, len_min_2        
      real(8) :: eps1, eps2                  
      integer :: i                           

      eps1 = 1d-20
      eps2 = 1d-99

      n_norm = sqrt(n_x ** 2 + n_y ** 2 + n_z ** 2)
      n_xn = n_x / n_norm
      n_yn = n_y / n_norm
      n_zn = n_z / n_norm


      len_from_face(1) = 0.25d0 * abs( &
                        (abs(Deter(x1-x0, x1-x2, x1-x3, y1-y0, y1-y2, y1-y3, z1-z0, z1-z2, z1-z3)) + eps1) / &
                        (abs(Deter(n_xn , x1-x2, x1-x3, n_yn , y1-y2, y1-y3, n_zn , z1-z2, z1-z3)) + eps2) + &
                        (abs(Deter(x2-x0, x2-x3, x2-x4, y2-y0, y2-y3, y2-y4, z2-z0, z2-z3, z2-z4)) + eps1) / &
                        (abs(Deter(n_xn , x2-x3, x2-x4, n_yn , y2-y3, y2-y4, n_zn , z2-z3, z2-z4)) + eps2) + &
                        (abs(Deter(x1-x0, x1-x3, x1-x4, y1-y0, y1-y3, y1-y4, z1-z0, z1-z3, z1-z4)) + eps1) / &
                        (abs(Deter(n_xn , x1-x3, x1-x4, n_yn , y1-y3, y1-y4, n_zn , z1-z3, z1-z4)) + eps2) + &
                        (abs(Deter(x1-x0, x1-x2, x1-x4, y1-y0, y1-y2, y1-y4, z1-z0, z1-z2, z1-z4)) + eps1) / &
                        (abs(Deter(n_xn , x1-x2, x1-x4, n_yn , y1-y2, y1-y4, n_zn , z1-z2, z1-z4)) + eps2))

      len_from_face(2) = 0.25d0 * abs( &
                        (abs(Deter(x2-x0, x2-x3, x2-x7, y2-y0, y2-y3, y2-y7, z2-z0, z2-z3, z2-z7)) + eps1) / &
                        (abs(Deter(n_xn , x2-x3, x2-x7, n_yn , y2-y3, y2-y7, n_zn , z2-z3, z2-z7)) + eps2) + &
                        (abs(Deter(x3-x0, x3-x7, x3-x6, y3-y0, y3-y7, y3-y6, z3-z0, z3-z7, z3-z6)) + eps1) / &
                        (abs(Deter(n_xn , x3-x7, x3-x6, n_yn , y3-y7, y3-y6, n_zn , z3-z7, z3-z6)) + eps2) + &
                        (abs(Deter(x2-x0, x2-x7, x2-x6, y2-y0, y2-y7, y2-y6, z2-z0, z2-z7, z2-z6)) + eps1) / &
                        (abs(Deter(n_xn , x2-x7, x2-x6, n_yn , y2-y7, y2-y6, n_zn , z2-z7, z2-z6)) + eps2) + &
                        (abs(Deter(x2-x0, x2-x3, x2-x6, y2-y0, y2-y3, y2-y6, z2-z0, z2-z3, z2-z6)) + eps1) / &
                        (abs(Deter(n_xn , x2-x3, x2-x6, n_yn , y2-y3, y2-y6, n_zn , z2-z3, z2-z6)) + eps2))

      len_from_face(3) = 0.25d0 * abs( &
                        (abs(Deter(x4-x0, x4-x3, x4-x7, y4-y0, y4-y3, y4-y7, z4-z0, z4-z3, z4-z7)) + eps1) / &
                        (abs(Deter(n_xn , x4-x3, x4-x7, n_yn , y4-y3, y4-y7, n_zn , z4-z3, z4-z7)) + eps2) + &
                        (abs(Deter(x3-x0, x3-x7, x3-x8, y3-y0, y3-y7, y3-y8, z3-z0, z3-z7, z3-z8)) + eps1) / &
                        (abs(Deter(n_xn , x3-x7, x3-x8, n_yn , y3-y7, y3-y8, n_zn , z3-z7, z3-z8)) + eps2) + &
                        (abs(Deter(x4-x0, x4-x7, x4-x8, y4-y0, y4-y7, y4-y8, z4-z0, z4-z7, z4-z8)) + eps1) / &
                        (abs(Deter(n_xn , x4-x7, x4-x8, n_yn , y4-y7, y4-y8, n_zn , z4-z7, z4-z8)) + eps2) + &
                        (abs(Deter(x4-x0, x4-x3, x4-x8, y4-y0, y4-y3, y4-y8, z4-z0, z4-z3, z4-z8)) + eps1) / &
                        (abs(Deter(n_xn , x4-x3, x4-x8, n_yn , y4-y3, y4-y8, n_zn , z4-z3, z4-z8)) + eps2))

      len_from_face(4) = 0.25d0 * abs( &
                        (abs(Deter(x5-x0, x5-x6, x5-x7, y5-y0, y5-y6, y5-y7, z5-z0, z5-z6, z5-z7)) + eps1) / &
                        (abs(Deter(n_xn , x5-x6, x5-x7, n_yn , y5-y6, y5-y7, n_zn , z5-z6, z5-z7)) + eps2) + &
                        (abs(Deter(x6-x0, x6-x7, x6-x8, y6-y0, y6-y7, y6-y8, z6-z0, z6-z7, z6-z8)) + eps1) / &
                        (abs(Deter(n_xn , x6-x7, x6-x8, n_yn , y6-y7, y6-y8, n_zn , z6-z7, z6-z8)) + eps2) + &
                        (abs(Deter(x5-x0, x5-x7, x5-x8, y5-y0, y5-y7, y5-y8, z5-z0, z5-z7, z5-z8)) + eps1) / &
                        (abs(Deter(n_xn , x5-x7, x5-x8, n_yn , y5-y7, y5-y8, n_zn , z5-z7, z5-z8)) + eps2) + &
                        (abs(Deter(x5-x0, x5-x6, x5-x8, y5-y0, y5-y6, y5-y8, z5-z0, z5-z6, z5-z8)) + eps1) / &
                        (abs(Deter(n_xn , x5-x6, x5-x8, n_yn , y5-y6, y5-y8, n_zn , z5-z6, z5-z8)) + eps2))

      len_from_face(5) = 0.25d0 * abs( &
                        (abs(Deter(x5-x0, x5-x1, x5-x4, y5-y0, y5-y1, y5-y4, z5-z0, z5-z1, z5-z4)) + eps1) / &
                        (abs(Deter(n_xn , x5-x1, x5-x4, n_yn , y5-y1, y5-y4, n_zn , z5-z1, z5-z4)) + eps2) + &
                        (abs(Deter(x1-x0, x1-x4, x1-x8, y1-y0, y1-y4, y1-y8, z1-z0, z1-z4, z1-z8)) + eps1) / &
                        (abs(Deter(n_xn , x1-x4, x1-x8, n_yn , y1-y4, y1-y8, n_zn , z1-z4, z1-z8)) + eps2) + &
                        (abs(Deter(x5-x0, x5-x4, x5-x8, y5-y0, y5-y4, y5-y8, z5-z0, z5-z4, z5-z8)) + eps1) / &
                        (abs(Deter(n_xn , x5-x4, x5-x8, n_yn , y5-y4, y5-y8, n_zn , z5-z4, z5-z8)) + eps2) + &
                        (abs(Deter(x5-x0, x5-x1, x5-x8, y5-y0, y5-y1, y5-y8, z5-z0, z5-z1, z5-z8)) + eps1) / &
                        (abs(Deter(n_xn , x5-x1, x5-x8, n_yn , y5-y1, y5-y8, n_zn , z5-z1, z5-z8)) + eps2))

      len_from_face(6) = 0.25d0 * abs( &
                        (abs(Deter(x1-x0, x1-x2, x1-x6, y1-y0, y1-y2, y1-y6, z1-z0, z1-z2, z1-z6)) + eps1) / &
                        (abs(Deter(n_xn , x1-x2, x1-x6, n_yn , y1-y2, y1-y6, n_zn , z1-z2, z1-z6)) + eps2) + &
                        (abs(Deter(x2-x0, x2-x6, x2-x5, y2-y0, y2-y6, y2-y5, z2-z0, z2-z6, z2-z5)) + eps1) / &
                        (abs(Deter(n_xn , x2-x6, x2-x5, n_yn , y2-y6, y2-y5, n_zn , z2-z6, z2-z5)) + eps2) + &
                        (abs(Deter(x1-x0, x1-x5, x1-x6, y1-y0, y1-y5, y1-y6, z1-z0, z1-z5, z1-z6)) + eps1) / &
                        (abs(Deter(n_xn , x1-x5, x1-x6, n_yn , y1-y5, y1-y6, n_zn , z1-z5, z1-z6)) + eps2) + &
                        (abs(Deter(x1-x0, x1-x2, x1-x5, y1-y0, y1-y2, y1-y5, z1-z0, z1-z2, z1-z5)) + eps1) / &
                        (abs(Deter(n_xn , x1-x2, x1-x5, n_yn , y1-y2, y1-y5, n_zn , z1-z2, z1-z5)) + eps2))



      len_min_1 = maxval(len_from_face)
      len_min_2 = len_min_1

      do i = 1, 6
         if (len_from_face(i) < len_min_1) then
            len_min_2 = len_min_1
            len_min_1 = len_from_face(i)
         else if (len_from_face(i) < len_min_2) then
            len_min_2 = len_from_face(i)
         end if
      end do

      length = len_min_1 + len_min_2
   end function Line_length_in_cell

   subroutine Mirror_image_3d(a, b, c, d, x1, y1, z1, x2, y2, z2)
      implicit none
      real(8), intent(in)  :: a, b, c, d  
      real(8), intent(in)  :: x1, y1, z1  
      real(8), intent(out) :: x2, y2, z2  
      real(8) :: l 

      l  = d - (x1 * a + y1 * b + z1 * c)
      x2 = x1 + 2d0 * l * a
      y2 = y1 + 2d0 * l * b
      z2 = z1 + 2d0 * l * c
   end subroutine Mirror_image_3d


   subroutine Plane_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3, normal_x, normal_y, normal_z)
      implicit none
      real(8), intent(in)  :: x1, y1, z1, x2, y2, z2, x3, y3, z3   
      real(8), intent(out) :: normal_x, normal_y, normal_z         
      real(8) :: dx1, dy1, dz1, dx2, dy2, dz2  
      real(8) :: norm  

      dx1 = x1 - x2
      dy1 = y1 - y2
      dz1 = z1 - z2
      dx2 = x3 - x1
      dy2 = y3 - y1
      dz2 = z3 - z1
      normal_x = dy1 * dz2 - dz1 * dy2
      normal_y = dz1 * dx2 - dx1 * dz2
      normal_z = dx1 * dy2 - dy1 * dx2
      norm = sqrt(normal_x ** 2 + normal_y ** 2 + normal_z ** 2)
      normal_x = normal_x / norm
      normal_y = normal_y / norm
      normal_z = normal_z / norm

  end subroutine Plane_normal




   subroutine Vertex_interp_3d(x0, y0, z0, i, j, k, iv1, jv1, kv1, iv2, jv2, kv2, &
        iv3, jv3, kv3, u_inter, u_orig, xint, yint, zint, w0, w1, w2, w3)
      implicit none
      real(8), intent(in)        :: x0, y0, z0, w0, w1, w2, w3
      real(8), intent(in)    :: u_orig(0:, 0:, 0:), xint(0:, 0:, 0:), yint(0:, 0:, 0:), &
                                    zint(0:, 0:, 0:)
      integer, intent(in)        :: i, j, k, iv1, jv1, kv1, iv2, jv2, kv2, iv3, jv3, kv3
      real(8), intent(out)       :: u_inter

      real(8)                    :: dx1, dx2, dx3, dy1, dy2, dy3, dz1, dz2, dz3, df1, df2, &
                                    df3, det, dx, dy, dz, beta1, beta2, beta3, beta0, fac0, fac1,&
                                    fac2, fac3, fsum


      dx1 = xint(iv1, jv1, kv1) - xint(i, j, k)
      dx2 = xint(iv2, jv2, kv2) - xint(i, j, k)
      dx3 = xint(iv3, jv3, kv3) - xint(i, j, k)
      dy1 = yint(iv1, jv1, kv1) - yint(i, j, k)
      dy2 = yint(iv2, jv2, kv2) - yint(i, j, k)
      dy3 = yint(iv3, jv3, kv3) - yint(i, j, k)
      dz1 = zint(iv1, jv1, kv1) - zint(i, j, k)
      dz2 = zint(iv2, jv2, kv2) - zint(i, j, k)
      dz3 = zint(iv3, jv3, kv3) - zint(i, j, k)

      df1 = u_orig(iv1, jv1, kv1) - u_orig(i, j, k)
      df2 = u_orig(iv2, jv2, kv2) - u_orig(i, j, k)
      df3 = u_orig(iv3, jv3, kv3) - u_orig(i, j, k)

      det = deter(dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3)

      dx = x0 - xint(i, j, k)
      dy = y0 - yint(i, j, k)
      dz = z0 - zint(i, j, k)

      if(det == 0) then
         u_inter = u_orig(i, j, k)
      else
         beta1 = deter(dx , dy , dz , dx2, dy2, dz2, dx3, dy3, dz3)
         beta2 = deter(dx1, dy1, dz1, dx , dy , dz , dx3, dy3, dz3)
         beta3 = deter(dx1, dy1, dz1, dx2, dy2, dz2, dx , dy , dz )

         beta0 = det - beta1
         beta0 = beta0 - beta2
         beta0 = beta0 - beta3

         fac0 = beta0 * w0
         fac1 = beta1 * w1
         fac2 = beta2 * w2
         fac3 = beta3 * w3

         fsum = fac0 + fac1
         fsum = fsum + fac2
         fsum = fsum + fac3


         if (fsum == 0) then
            u_inter = (beta0*u_orig(i, j, k)+ &
                       beta1*u_orig(iv1, jv1, kv1) + &
                       beta2*u_orig(iv2, jv2, kv2) + &
                       beta3*u_orig(iv3, jv3, kv3)) / det
         else
            fac0 = fac0 / fsum
            fac1 = fac1 / fsum
            fac2 = fac2 / fsum
            fac3 = fac3 / fsum
            u_inter = fac0*u_orig(i, j, k) + &
                      fac1*u_orig(iv1, jv1, kv1) + &
                      fac2*u_orig(iv2, jv2, kv2) + &
                      fac3*u_orig(iv3, jv3, kv3)
         end if
      end if

      return
      end subroutine Vertex_interp_3d

      subroutine Layer_cut_line_3d(y11, x11, y12, x12, layer_parm, n_part, yh, xh)
          use constants_module, only : PI, N_SEGMENT_MAX
          implicit none
          integer, intent(in)  :: n_part                
          real(8), intent(in)  :: y11, x11, y12, x12    
          real(8), intent(in)  :: layer_parm(n_part, 4) 
          real(8), intent(out) :: yh, xh                

          real(8) :: x21, y21, x22, y22                  
          real(8) :: x1, x2, y1, y2, a, b, al, be        
          real(8) :: dx_dy1, dx_dy2, x_c, y_c, x_d, y_d  
          real(8) :: tt1, tt2, c, d ,signy               
          real(8) :: xl_cut, yl_cut, x_cut, y_cut        
          integer :: n_cut                               
          integer :: is_cut                              
          integer :: n                                   



          n_cut = -1    
          do n = 1, 1
              y21 = layer_parm(n, 1)  
              x21 = layer_parm(n, 2)  
              y22 = layer_parm(n, 3)  
              x22 = layer_parm(n, 4)  

              y21 = 1.005d0*y21 - 0.005d0*y22   
              x21 = 1.005d0*x21 - 0.005d0*x22
              y22 = 1.005d0*y22 - 0.005d0*y21
              x22 = 1.005d0*x22 - 0.005d0*x21
              call Cut_line(x11, y11, x12, y12, x21, y21, x22, y22, is_cut, xl_cut, yl_cut, 1) 

              if (is_cut == 1) then 
                  n_cut = n
                  x_cut = xl_cut
                  y_cut = yl_cut
              end if
          end do
          if (n_cut == -1) then
              write(*,*) "geometry.f90: DOESNT INTERSECT"
              write(*,*), x11,y11,x12,y12,x21,y21,x22,y22
              stop
              return
          end if

          xh = x_cut
          yh = y_cut
          return



      end subroutine Layer_cut_line_3d

   subroutine Layer_cut_elipse_3d(y11, x11, y12, x12, layer_parm, n_part, yh, xh)
      use constants_module, only : PI, N_SEGMENT_MAX
      implicit none
      integer, intent(in)  :: n_part                
      real(8), intent(in)  :: y11, x11, y12, x12    
      real(8), intent(in)  :: layer_parm(n_part, 11) 
      real(8), intent(out) :: yh, xh                

      real(8) :: x21, y21, x22, y22                  
      real(8) :: x1, x2, y1, y2, a, b, al, be        
      real(8) :: dx_dy1, dx_dy2, xc, yc, x_d, y_d  
      real(8) :: tt1, tt2, c, d ,signy               
      real(8) :: xl_cut, yl_cut, x_cut, y_cut        
      integer :: n_cut                               
      integer :: is_cut                              
      integer :: n                                   


      n_cut = 1    

      call Layer_cut_line_3d(y11,x11,y12,x12,layer_parm, n_part, y_cut, x_cut)

      x1 = 0.90*x_cut + 0.1*x11
      y1 = 0.90*y_cut + 0.1*y11
      x2 = x12
      y2 = y12
      a  = layer_parm(n_cut, 5)
      b  = layer_parm(n_cut, 6)
      al = layer_parm(n_cut, 7)
      be = layer_parm(n_cut, 8)
      yc = layer_parm(n_cut, 9)
      xc = layer_parm(n_cut, 10)
      call Elipse_cut_line(a, b, al, be, yc, xc, y1, x1, y2, x2, yh, xh)

   end subroutine Layer_cut_elipse_3d


   subroutine Calculate_constant(arr, x1, x2, from, to, tmp)
      implicit none
      real(8), dimension(:), pointer, intent(inout)   :: arr
      integer              , intent(in)      :: from, to
      real(8)              , intent(in)      :: x1, x2
      integer, optional :: tmp 
      integer :: from_tmp
      integer                                :: i


      if (.not. present(tmp)) then
         from_tmp = from  
      else
         from_tmp = tmp   
      end if

      do i = from, to
         arr(i) = x1 + (x2 - x1) * dble(i - from_tmp) / dble(to - from_tmp)
      end do

   end subroutine Calculate_constant

   subroutine Calculate_geometry_series(arr, dr, qq, from, to, tmp)
      implicit none
      real(8), dimension(:), pointer, intent(inout)   :: arr
      integer              , intent(in)      :: from, to
      real(8)              , intent(in)      :: qq, dr
      integer, optional :: tmp 
      integer :: from_tmp
      integer                                :: i

      if (.not. present(tmp)) then
         from_tmp = from  
      else
         from_tmp = tmp   
      end if

      do i = from, to
         arr(i) = arr(i-1) + dr * qq **dble(i - from)
      end do

   end subroutine Calculate_geometry_series

   subroutine Triangle_method(x_top_1_1    , y_top_1_1    , z_top_1_1,&
                           x_top_1_nyp  , y_top_1_nyp  , z_top_1_nyp, &
                           x_top_nxp_1  , y_top_nxp_1  , z_top_nxp_1, &
                           x_top_nxp_nyp, y_top_nxp_nyp, z_top_nxp_nyp,&
                           face_i_1, face_i_nyp, face_1_j, face_nxp_j,     &
                           x_top, y_top, z_top)
      implicit none
      real(8), intent(in) :: x_top_1_1    , y_top_1_1    , z_top_1_1,&
                           x_top_1_nyp  , y_top_1_nyp  , z_top_1_nyp, &
                           x_top_nxp_1  , y_top_nxp_1  , z_top_nxp_1, &
                           x_top_nxp_nyp, y_top_nxp_nyp, z_top_nxp_nyp,&
                           face_i_1, face_i_nyp, face_1_j, face_nxp_j
      real(8), intent(out) :: x_top, y_top, z_top
      real(8) :: tet, tet_i, alpha, beta, x1, y1, z1, x2, y2, z2, &
                 x3, y3, z3, x4, y4, z4, d, amat(3, 3), bmat(3), eps, gama
      integer :: indx(3)

      eps = 1d-10

      tet = acos(x_top_1_1*x_top_nxp_1 + y_top_1_1*y_top_nxp_1 + z_top_1_1*z_top_nxp_1)
      tet_i = face_i_1 * tet

      alpha = sin(tet - tet_i) / sin(tet)
      beta  = sin(tet_i) / sin(tet)
      x1    = alpha*x_top_1_1 + beta*x_top_nxp_1
      y1    = alpha*y_top_1_1 + beta*y_top_nxp_1
      z1    = alpha*z_top_1_1 + beta*z_top_nxp_1

      tet   = acos(x_top_nxp_nyp*x_top_nxp_1 + y_top_nxp_nyp*y_top_nxp_1 + z_top_nxp_nyp*z_top_nxp_1)
      tet_i = face_nxp_j*tet
      alpha = sin(tet - tet_i) / sin(tet)
      beta  = sin(tet_i) / sin(tet)
      x2    = alpha*x_top_nxp_1 + beta*x_top_nxp_nyp
      y2    = alpha*y_top_nxp_1 + beta*y_top_nxp_nyp
      z2    = alpha*z_top_nxp_1 + beta*z_top_nxp_nyp

      tet   = acos(x_top_nxp_nyp*x_top_1_nyp + y_top_nxp_nyp*y_top_1_nyp + z_top_nxp_nyp*z_top_1_nyp)
      tet_i = face_i_nyp * tet
      alpha = sin(tet - tet_i) / sin(tet)
      beta  = sin(tet_i) / sin(tet)
      x3    = alpha*x_top_1_nyp + beta*x_top_nxp_nyp
      y3    = alpha*y_top_1_nyp + beta*y_top_nxp_nyp
      z3    = alpha*z_top_1_nyp + beta*z_top_nxp_nyp

      tet   = acos(x_top_1_1*x_top_1_nyp + y_top_1_1*y_top_1_nyp + z_top_1_1*z_top_1_nyp)
      tet_i = face_1_j * tet
      alpha = sin(tet - tet_i) / sin(tet)
      beta  = sin(tet_i) / sin(tet)
      x4    = alpha*x_top_1_1 + beta*x_top_1_nyp
      y4    = alpha*y_top_1_1 + beta*y_top_1_nyp
      z4    = alpha*z_top_1_1 + beta*z_top_1_nyp

      if (abs(face_i_1+1) < eps) then
        x_top = x4
        y_top = y4
        z_top = z4
      elseif (abs(face_i_1-1d0) < eps) then
        x_top = x2
        y_top = y2
        z_top = z2
      elseif (abs(face_1_j) < eps) then
        x_top = x1
        y_top = y1
        z_top = z1
      elseif (abs(face_1_j-1d0) < eps) then
        x_top = x3
        y_top = y3
        z_top = z3
      else
        amat(1,1) = 1d0
        amat(1,2) =   x2*x4 + y2*y4 + z2*z4
        amat(1,3) = -(x2*x3 + y2*y3 + z2*z3)
        amat(2,1) =   x2*x4 + y2*y4 + z2*z4
        amat(2,2) = 1d0
        amat(2,3) = -(x4*x3 + y4*y3 + z4*z3)
        amat(3,1) =   x2*x3 + y2*y3 + z2*z3
        amat(3,2) =   x4*x3 + y4*y3 + z4*z3
        amat(3,3) = -1d0

        bmat(1) = x2*x1 + y2*y1 + z2*z1
        bmat(2) = x4*x1 + y4*y1 + z4*z1
        bmat(3) = x3*x1 + y3*y1 + z3*z1

        call ludcmp(amat, 3, 3, indx, d)
        call lubksb(amat, 3, 3, indx, bmat)

        gama = bmat(3)
        alpha = 1d0
        x_top = x1 + gama*x3
        y_top = y1 + gama*y3
        z_top = z1 + gama*z3
      endif

      return
      end subroutine Triangle_method


   subroutine ludcmp(a, n, np, indx, d)
      implicit none

      integer, intent(in)    :: n, np
      real(8), intent(inout) :: a(np,np)
      integer, intent(out)   :: indx(n)
      real(8), intent(out)   :: d

      real(8), parameter :: eps = 1.0d-20
      integer, parameter :: nmax = 100
      integer :: i, j, k, imax
      real(8) :: vv(nmax), aamax, summ, dum

      d = 1d0
      do i = 1, n
         aamax = 0d0
         do j = 1, n
            if(abs(a(i, j)) > aamax) aamax = abs(a(i, j))
         end do
         if (aamax.eq.0.) stop 'ludcmp. singular.'
         vv(i) = 1d0 / aamax
      end do
      do j = 1, n
         if (j > 1) then
            do i = 1, j - 1
               summ = a(i, j)
               if (i > 1) then
                  do k=1,i-1
                     summ = summ - a(i, k) * a(k, j)
                  end do
                  a(i,j) = summ
               endif
            end do
         end if
         aamax = 0d0
         do i = j, n
            summ = a(i, j)
            if (j > 1) then
               do k = 1, j-1
                  summ = summ - a(i, k)*a(k, j)
               end do
               a(i,j) = summ
            end if
            dum = vv(i) * abs(summ)
            if (dum >= aamax) then
               imax = i
               aamax = dum
            end if
         end do
         if (j /= imax) then
            do k = 1, n
               dum = a(imax, k)
               a(imax, k) = a(j, k)
               a(j, k) = dum
            end do
            d = -d
            vv(imax) = vv(j)
         end if
         indx(j) = imax
         if (j /= n) then
            if (a(j, j) == 0d0 ) a(j,j) = eps
            dum = 1d0 / a(j, j)
            do i = j+1, n
               a(i, j) = a(i, j)*dum
            end do
         end if
      end do
      if ( a(n, n) == 0d0) a(n,n) = eps
      return
   end subroutine ludcmp

   subroutine lubksb(a, n, np, indx, b)
      implicit none
      real(8), intent(in)    :: a(np,np)
      real(8), intent(inout) :: b(n)
      integer, intent(in)    :: n, np, indx(n)
      integer :: ii, ll, i, j
      real(8) :: summ
      ii = 0
      do i = 1, n
         ll = indx(i)
         summ = b(ll)
         b(ll) = b(i)
         if (ii /= 0) then
            do j = ii, i-1
               summ = summ - a(i, j)*b(j)
            end do
         else if(summ /= 0d0) then
            ii = i
         end if
         b(i) = summ
      end do

      do i = n, 1, -1
         summ = b(i)
         if (i < n) then
            do j = i+1, n
               summ = summ - a(i, j)*b(j)
            end do
         end if
         b(i) = summ / a(i,i)
      end do
      return
    end subroutine lubksb

   function Create_start_index_layer(num_cells, length, index1) result (index_arr)
      implicit none
      integer      , intent (in)                                 :: length, index1
      integer, dimension(:), allocatable   , intent(in)          :: num_cells
      integer, dimension(length + 1, index1)                 :: index_arr

      integer                                           :: i, j

      do j = 1, index1 ! index = 1 always
         index_arr(1, j) = 1
         do i = 1, length
            index_arr(i + 1, j) = index_arr(i, j) + num_cells(i)
         end do
      end do
   end function Create_start_index_layer

    function Check_corner_vertex(ivirt, jvirt, kvirt, virt_nxp, virt_nyp, virt_nzp) result(is_corner)
        integer, intent(in) :: ivirt
        integer, intent(in) :: jvirt
        integer, intent(in) :: kvirt
        integer, intent(in) :: virt_nxp
        integer, intent(in) :: virt_nyp
        integer, intent(in) :: virt_nzp
        logical :: is_corner
        is_corner = .false.

        if ( (ivirt == 1 .or. ivirt == virt_nxp) .and.&
            (jvirt == 1 .or. jvirt == virt_nyp) .and.&
            (kvirt == 1 .or. kvirt == virt_nzp)) then
            is_corner = .true.
        end if
    end function
end module geometry_module
