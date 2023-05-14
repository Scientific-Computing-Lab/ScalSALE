
submodule (diagnostic_module) textual_diagnostic_hdf5_module
use hydro_step_module, only : hydro_step_t
use hdf5
implicit none

contains

    module subroutine textual_diagnostic_hdf5_init (this, hydro, time, diag_num, main_file_id)
        implicit none

        class (textual_diagnostic_hdf5_t), intent(in out) :: this
        type(hydro_step_t), pointer, intent(in)           :: hydro         
        type(time_t), pointer      , intent(in)           :: time          
        integer                    , intent(in)           :: diag_num      
        integer(hid_t)             , intent(in)           :: main_file_id
        integer                                           :: error         
   
        this%hydro_step  => hydro
        this%time  => time
        this%diagnostic_counter = 0
        this%main_file_id = main_file_id


        call h5gcreate_f(this%main_file_id, "/" // this%file_prefix, this%gid, error)




    end subroutine textual_diagnostic_hdf5_init

   module subroutine textual_diagnostic_hdf5_apply (this)
      implicit none
      class (textual_diagnostic_hdf5_t), intent(in out)               :: this
      character(3) :: num_string_3      
      real(8), dimension(:,:,:), pointer    :: ptr_x, ptr_y 
      integer     ::   error 
      integer(hsize_t), dimension(1) :: dims_1_real    
      integer(hsize_t), dimension(1) :: dims_3d_real   
      real(8), dimension(1) :: data_1_real             
      real(8), dimension(:), allocatable :: data_3d_real            
      integer     ::   rank                         

      this%diagnostic_counter = this%diagnostic_counter + 1
      write (num_string_3, "(I3)") this%diagnostic_counter
      1002  format(100000(1PE25.17))

      select case( this%file_prefix)

         case ('time')

            rank = 1

            data_1_real(1) = this%time%time_passed
            dims_1_real(1) = 1

            call h5screate_simple_f(rank, dims_1_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_1_real, dims_1_real, error)

            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

         case ('position_x')
            call this%hydro_step%Point_to_mesh_data(ptr_x, ptr_y)
            rank = 1

            dims_3d_real(1) = size(ptr_x(1:, 1, 1))
            allocate(data_3d_real(dims_3d_real(1)))
            data_3d_real = ptr_x(1:,1,1)

            call h5screate_simple_f(rank, dims_3d_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_3d_real, dims_3d_real, error)

            deallocate(data_3d_real)
            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

          case ('position_y')
            call this%hydro_step%Point_to_mesh_data(ptr_x, ptr_y)
            rank = 1

            dims_3d_real(1) = size(ptr_y(1:, 1, 1))
            allocate(data_3d_real(dims_3d_real(1)))
            data_3d_real = ptr_y(1:,1,1)

            call h5screate_simple_f(rank, dims_3d_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_3d_real, dims_3d_real, error)

            deallocate(data_3d_real)
            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

          case ('velocity_x')
            call this%hydro_step%Point_to_velocity_data(ptr_x, ptr_y)
            rank = 1

            dims_3d_real(1) = size(ptr_x(1:, 1, 1))
            allocate(data_3d_real(dims_3d_real(1)))
            data_3d_real = ptr_x(1:,1,1)

            call h5screate_simple_f(rank, dims_3d_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_3d_real, dims_3d_real, error)

            deallocate(data_3d_real)
            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

          case ('velocity_y')
            call this%hydro_step%Point_to_velocity_data(ptr_x, ptr_y)
            rank = 1

            dims_3d_real(1) = size(ptr_y(1:, 1, 1))
            allocate(data_3d_real(dims_3d_real(1)))
            data_3d_real = ptr_y(1:,1,1)

            call h5screate_simple_f(rank, dims_3d_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_3d_real, dims_3d_real, error)

            deallocate(data_3d_real)
            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

          case ('pressure')
            call this%hydro_step%Point_to_pressure_data(ptr_x)
            rank = 1

            dims_3d_real(1) = size(ptr_x(1:, 1, 1))
            allocate(data_3d_real(dims_3d_real(1)))
            data_3d_real = ptr_x(1:,1,1)

            call h5screate_simple_f(rank, dims_3d_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_3d_real, dims_3d_real, error)

            deallocate(data_3d_real)
            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

          case ('density')
            call this%hydro_step%Point_to_density_data(ptr_x)
            rank = 1

            dims_3d_real(1) = size(ptr_x(1:, 1, 1))
            allocate(data_3d_real(dims_3d_real(1)))
            data_3d_real = ptr_x(1:,1,1)

            call h5screate_simple_f(rank, dims_3d_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_3d_real, dims_3d_real, error)

            deallocate(data_3d_real)
            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

          case ('sie')
            call this%hydro_step%Point_to_sie_data(ptr_x)
            rank = 1

            dims_3d_real(1) = size(ptr_x(1:, 1, 1))
            allocate(data_3d_real(dims_3d_real(1)))
            data_3d_real = ptr_x(1:,1,1)

            call h5screate_simple_f(rank, dims_3d_real, this%dspace_id, error)

            call h5dcreate_f(this%main_file_id, "/" // trim(this%file_prefix) // "/" // trim(num_string_3), H5T_NATIVE_DOUBLE, this%dspace_id, &
                             this%dset_id, error)

            call h5dwrite_f(this%dset_id, H5T_NATIVE_DOUBLE, data_3d_real, dims_3d_real, error)

            deallocate(data_3d_real)
            call h5dclose_f(this%dset_id, error)

           call h5sclose_f(this%dspace_id, error)

          case default
      end select

   end subroutine textual_diagnostic_hdf5_apply

    module subroutine textual_diagnostic_hdf5_close (this)
            implicit none
            class (textual_diagnostic_hdf5_t), intent(in out)               :: this
            integer     ::   error 

           call h5gclose_f(this%gid, error)

    end subroutine textual_diagnostic_hdf5_close


end submodule textual_diagnostic_hdf5_module
