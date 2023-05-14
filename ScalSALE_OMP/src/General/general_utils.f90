module general_utils_module
   implicit none

contains

   function Get_axises_num(d1, d2, d3) result (return_value)
      integer , intent(in) :: d1
      integer , intent(in) :: d2
      integer , intent(in) :: d3
      integer              :: return_value

      return_value = 0
      if (d1 > 1) return_value = return_value + 1
      if (d2 > 1) return_value = return_value + 1
      if (d3 > 1) return_value = return_value + 1
   end function

   subroutine Int2d(x1, y1, x2, y2, x3, y3, xn, yn, w1, w2, w3, factor1, factor2, factor3)
      real(8), intent(in) :: x1,x2,x3,y1,y2,y3,w1,w2,w3,xn,yn
      real(8), intent(out) :: factor1, factor2, factor3
      real(8) :: det,beta1,beta2,beta3,factor_sum

      det = Z_vector_product(x1,y1,x2,y2,x3,y3)
      if(abs(det) < 1.d-40) then
        beta1 = 1.d0
        beta2 = 0.d0
        beta3 = 0.d0
      else
        beta1 = Z_vector_product(xn, yn, x2, y2, x3, y3) / det
        beta2 = Z_vector_product(x1, y1, xn, yn, x3, y3) / det
        beta3 = Z_vector_product(x1, y1, x2, y2, xn, yn) / det
      end if
      factor1 = beta1 * w1
      factor2 = beta2 * w2
      factor3 = beta3 * w3
      factor_sum = factor1 + factor2 + factor3
      if (abs(factor_sum) < 1d-30) then
        factor1 = beta1
        factor2 = beta2
        factor3 = beta3
      else
        factor1 = factor1 / factor_sum
        factor2 = factor2 / factor_sum
        factor3 = factor3 / factor_sum
      end if
      return
   end subroutine Int2d

   function Z_vector_product(xl1, yl1, xl2, yl2, xl3, yl3) result (return_value)
      real(8), intent(in) :: xl1, xl2, xl3, yl1, yl2, yl3
      real(8) return_value

      return_value = (xl3 - xl2) * (yl1 - yl2) - (xl1 - xl2) * (yl3 - yl2)
      return
   end function Z_vector_product

   function Deter(a11, a12, a13, a21, a22, a23, a31, a32, a33) result (return_value)
      real(8), intent(in) :: a11, a12, a13, a21, a22, a23, a31, a32, a33  
      real(8) return_value

      return_value = a11 * (a22 * a33 - a23 * a32) &
                    -a12 * (a21 * a33 - a23 * a31) &
                    +a13 * (a21 * a32 - a22 * a31)
   end function Deter

   function Compare_string(s1, s2) result (return_value)
      implicit none
      character(len=1), dimension(:), allocatable, intent(in) :: s1
      character(len=*),  intent(in) :: s2

      logical return_value
      integer :: i, length

      return_value = .true.

      do i = 1, size(s1)
         if (s1(i) /= s2(i:i)) then
            return_value = .false.
         end if
      end do
   end function Compare_string

    function int2str(number) result (str)
        implicit none
        integer, intent(in) :: number
        character(len=ceiling((log10(dble(number+1))))) ::str

        if (number < 10 ) then
            write(str, '(I1)') number
        elseif (number < 100) then
            write(str, '(I2)') number
        elseif (number < 1000) then
            write(str, '(I3)') number
        elseif (number < 10000) then
            write(str, '(I4)') number
        end if
    end function int2str

    function Lower(strIn) result(strOut)
     implicit none
     character(len=*), intent(in) :: strIn
     character(len=len(strIn)) :: strOut
     integer :: i,j

     do i = 1, len(strIn)
          j = iachar(strIn(i:i))
          if (j>= iachar("A") .and. j<=iachar("Z") ) then
               strOut(i:i) = achar(iachar(strIn(i:i))+32)
          else
               strOut(i:i) = strIn(i:i)
          end if
     end do

    end function Lower

    function Str_eqv(s1, s2) result(same)
        implicit none
        character(len=*), intent(in) :: s1
        character(len=*), intent(in) :: s2
        character(len=1) :: c1, c2  

        logical :: same
        integer :: i,j

        same = .false.

        if (len_trim(s1) .ne. len_trim(s2)) then
            return
        end if

        do i=1, min(len_trim(s1), len_trim(s2))
            c1 = lower(s1(i:i))
            c2 = lower(s2(i:i))
            if (c1 .ne. c2) then
                return
            end if
        end do

        same = .true.
        return
    end function str_eqv

    subroutine Error_msg(s1)
        implicit none
        character(len=*), intent(in) :: s1

        write(*,*) "Critical Error"
        write(*,*) s1
        stop
    end  subroutine Error_msg

    subroutine Warn_msg(s1)
        implicit none
        character(len=*), intent(in) :: s1

        write(*,*) "Warning"
        write(*,*) s1

    end  subroutine warn_msg

    function Str2int(str) result(int1)
        implicit none
        character(len=*),intent(in) :: str
        integer         :: int1

        read(str,*)  int1
    end function Str2int

    function Str2dble(str1) result(dble1)
        implicit none
        character(len=*),intent(in) :: str1
        real(8)         :: dble1
        read(str1,*)  dble1
    end function Str2dble

end module general_utils_module
