
submodule (diagnostic_module) plot_diagnostic_module
   use hydro_step_module, only : hydro_step_t
   implicit none

contains

   module subroutine plot_diagnostic_init (this, hydro, time, diag_num)
     implicit none

     class (plot_diagnostic_t)  , intent(in out) :: this
     type(hydro_step_t), pointer, intent(in)     :: hydro           
     type(time_t), pointer      , intent(in)     :: time            
     integer                    , intent(in)     :: diag_num        

     this%hydro_step  => hydro
     this%time  => time
     this%diagnostic_counter = 0
     this%diagnostics_group_file = diag_num

    end subroutine plot_diagnostic_init

   module subroutine plot_diagnostic_apply(this)
      use advect_module, only: advect_t
      implicit none
      class (plot_diagnostic_t), intent(in out) :: this  

      CHARACTER*20 NAMFIL
      CHARACTER*30 file_name
      CHARACTER*1 CH1
      CHARACTER*2 CH2
      CHARACTER*3 CH3
      CHARACTER*20 CH4
      CHARACTER*20 CONT_NAME(20), MILA_HALF,MILA_END

      real(8) :: x0, x1, x2, y0, y1, y2, tjnp, a, b, side, uvmax, c
      real(8) :: xx, yy, si, uu, vv, vofmax, fac, emf, emf1
      integer :: i, j, ii, ij, is1, is2, is3, iii, n_inter, nvof, ncont, numb, i_inter, &
                  ivofmax, nsc, ii1, ivofmax1, ivofmax2, iim, is, iv, jv, max_iter, &
                  ii2, nxp, nyp
      real(8) :: pressure_min, pressure_max, arti_visc_min, arti_visc_max, density_min, density_max,&
                 temperature_min, temperature_max
      real(8) :: p_val, d_val, a_val, t_val

      real(8), dimension(:, :, :), pointer              :: x, y, tot_vof, density, temperature, pressure, arti_visc
      real(8), dimension(:, :, :), pointer              :: vel_x, vel_y, num_mats, mat_id
      type(advect_t)            , pointer              :: advect

      call this%hydro_step%Point_to_artificial_viscosity_data(arti_visc)
      call this%hydro_step%Point_to_density_data(density)
      call this%hydro_step%Point_to_vof_data(tot_vof)
      call this%hydro_step%Point_to_temperature_data(temperature)
      call this%hydro_step%Point_to_pressure_data(pressure)
      call this%hydro_step%Point_to_num_mat_cells_data(num_mats)
      call this%hydro_step%Point_to_mesh_data(x, y)
      call this%hydro_step%Point_to_velocity_data(vel_x, vel_y)
      call this%hydro_step%Point_to_mat_id_data(mat_id)
      call this%hydro_step%Point_to_advect(advect)


      this%diagnostic_counter = this%diagnostic_counter + 1

      if (this%diagnostic_counter > 99) then
         write(file_name, "(A5,I3)") 'plot.', this%diagnostic_counter
      elseif (this%diagnostic_counter > 9) then
      write(file_name, "(A5,I2)") 'plot.', this%diagnostic_counter
      else
         write(file_name, "(A5,I1)") 'plot.', this%diagnostic_counter
      end if
      open(unit=this%diagnostics_group_file, file=file_name, status="replace", action="write")

      emf = this%hydro_step%emf
      emf1 = 1 - emf
      nxp = this%hydro_step%nxp
      nyp = this%hydro_step%nyp

      ncont = 4
      x0 = 1d38
      x1 = -1d38
      y0 = 1d38
      y1 = -1d38
      pressure_min = 1d38
      pressure_max = -1d38
      density_min = 1d38
      density_max = -1d38
      arti_visc_min = 1d38
      arti_visc_max = -1d38
      temperature_min = 1d38
      temperature_max = -1d38
      uvmax = 0d0
      nvof = 0

      do j=1, nyp
         do i=1, nxp
            if ((i < nxp) .and. (j < nyp)) then
               if ((tot_vof(i, j, 1) > emf) .and. (tot_vof(i, j, 1) < emf1)) nvof = nvof + 1
               if (num_mats(i, j, 1) > 1) nvof = nvof + num_mats(i, j, 1) - 1
            end if
            x0 = min(x0, x(i, j, 1))
            x1 = max(x1, x(i, j, 1))
            y0 = min(y0, y(i, j, 1))
            y1 = max(y1, y(i, j, 1))
            uvmax = max(uvmax, sqrt(vel_x(i, j, 1) * vel_x(i, j, 1) + vel_y(i, j, 1) * vel_y(i, j, 1)))


            if ((tot_vof(i, j, 1) < emf) .or. (i == nxp).or.(j == nyp)) cycle

            pressure_min = min(pressure_min, pressure(i, j, 1))
            pressure_max = max(pressure_max, pressure(i, j, 1))

            arti_visc_min = min(arti_visc_min, arti_visc(i, j, 1))
            arti_visc_max = max(arti_visc_max, arti_visc(i, j, 1))

            temperature_min = min(temperature_min, temperature(i, j, 1))
            temperature_max = max(temperature_max, temperature(i, j, 1))

            density_min = min(density_min, density(i, j, 1))
            density_max = max(density_max, density(i, j, 1))
         end do
      end do
      write(this%diagnostics_group_file, 1001) "T=", this%time%time_passed, "NCYC=", this%hydro_step%ncyc, "sod"
1001  format(A5, 1PE12.5, A7, I7, 2H  , A50)
      if (abs(x0) < 1d-30) x0 = 0d0
      x0 = max(x0, -1.0D99)
      x1 = min(x1, 1.0D99)
      y0 = max(y0, -1.0D99)
      y1 = min(y1, 1.0D99)
      n_inter=0
      write(this%diagnostics_group_file, 1002) ncont, nxp, nyp, nvof, 13, 1, x0, x1, y0, y1, uvmax
1002  format(I2,2I5,I7,2I3,5(1PE14.5E3))
      write(this%diagnostics_group_file, 3333) "ARTI_VISC", arti_visc_min, arti_visc_max
      write(this%diagnostics_group_file, 3333) "DENSITY", density_min, density_max
      write(this%diagnostics_group_file, 3333) "PRESSURE", pressure_min, pressure_max
      write(this%diagnostics_group_file, 3333) "TEMPERATURE", temperature_min, temperature_max
3333  format(A20,2(1PE14.5E3))

      write(this%diagnostics_group_file, 1003) 1
      write(this%diagnostics_group_file, 1003) 2
      write(this%diagnostics_group_file, 1003) 8
      write(this%diagnostics_group_file, 1003) 5
      write(this%diagnostics_group_file, 1003) 6
      write(this%diagnostics_group_file, 1003) 3
      write(this%diagnostics_group_file, 1003) 4
      write(this%diagnostics_group_file, 1003) 7
      write(this%diagnostics_group_file, 1003) 9
      write(this%diagnostics_group_file, 1003) 10
      write(this%diagnostics_group_file, 1003) 11
      write(this%diagnostics_group_file, 1003) 12
      write(this%diagnostics_group_file, 1003) 13
1003  format(3I2)

      do j = 1, nyp
         do i = 1, nxp
            nsc = mod(int(mat_id(i, j, 1)) - 1, 14) + 1
            if (num_mats(i, j, 1) >= 2) nsc = 8
            if (tot_vof(i, j, 1) < emf1) nsc = 6
            if (tot_vof(i, j, 1) < emf) nsc = 5
            xx = x(i, j, 1)
            if (abs(xx) < 1d-30) xx = 0d0
            yy = y(i, j, 1)
            if (abs(yy) < 1d-30) yy = 0d0
            write(this%diagnostics_group_file, 1004) xx, yy, nsc
         end do
      end do
1004  format(2(1PE14.5E3), I3)

      if (nvof >= 1) then
         do i = 1, nxp -1
            do j = 1, nyp - 1
               IF ((tot_vof(i, j, 1) > emf) .and. (tot_vof(i, j, 1) > emf1)) then
                     call advect%Calculate_vof_line_new(x, y, 0d0, 0, i, j, 1, nxp - 1, nyp - 1, 0, a, b, c, x1, y1, x2, y2, max_iter)

                  write(this%diagnostics_group_file, 1005) x2, y2, x1, y1
               end if
               if (num_mats(i, j, 1) < 2) cycle
               numb = 0
!               do ii = 1, this%hydro_step%nmats
!!                  if (mat_vof(ii, i, j, 1) < emf) cycle
!                  numb = numb + 1
!                  if (numb >= num_mats(i, j, 1)) cycle
!                     call advect%Calculate_vof_line_new(x, y, 0d0, ii, i, j, 1, nxp - 1, nyp - 1, 0, a, b, c, x1, y1, x2, y2, max_iter)
!                  if (abs(x1) < 1d-30) x1 = 0d0
!                  if (abs(x2) < 1d-30) x2 = 0d0
!                  if (abs(y1) < 1d-30) y1 = 0d0
!                  if (abs(y2) < 1d-30) y2 = 0d0
!                  write(this%diagnostics_group_file, 1005) x2,y2,x1,y1
!               end do
            end do
         end do
         1005  format(4(1PE14.5E3))
      end if
      write(this%diagnostics_group_file, 1011)
1011  format('VELOCITIES')
      do j = 1, nyp
         do i = 1, nxp
            uu = vel_x(i, j, 1)
            if (abs(uu) < 1d-30) uu = 0d0
            vv = vel_y(i, j, 1)
            if (abs(vv) < 1d-30) vv = 0d0
            write(this%diagnostics_group_file, 1012) uu, vv
         end do
      end do
1012  format(2(1PE14.5E3))
      write(this%diagnostics_group_file, 1009)
1009  format('CONTOURS')
      do j = 1, nyp
         do i = 1, nxp

            a_val = min(max(arti_visc(i, j, 1), -1d99), 1d99)
            d_val = min(max(density(i, j, 1), -1d99), 1d99)
            p_val = min(max(pressure(i, j, 1), -1d99), 1d99)
            t_val = min(max(temperature(i, j, 1), -1d99), 1d99)

            write(this%diagnostics_group_file, 1010) a_val
            write(this%diagnostics_group_file, 1010) d_val
            write(this%diagnostics_group_file, 1010) p_val
            write(this%diagnostics_group_file, 1010) t_val
         end do
      end do
1010  format(30(1PE14.5E3))
      close(this%diagnostics_group_file)

   end subroutine plot_diagnostic_apply

   module subroutine plot_diagnostic_close (this)
            implicit none
            class (plot_diagnostic_t), intent(in out)               :: this

    end subroutine plot_diagnostic_close

end submodule plot_diagnostic_module
