module assertion_mod

  use gengrid_kinds, only : dbl_kind, int_kind, real_kind

  implicit none

  private
  public :: assert_equal

  interface assert_equal
     module procedure assert_int_scalar
     module procedure assert_real_scalar
     module procedure assert_double_scalar
     module procedure assert_double_1d
  end interface assert_equal

contains
  subroutine assert_int_scalar(actual, expected, msg, rc, returnmsg)
    integer(int_kind), intent(in)  :: actual, expected
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = (actual == expected)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       write(returnmsg, '(2(a,i0))') "Fail: " // trim(msg) // " | Expected ", expected, ", got ", actual
    end if
  end subroutine assert_int_scalar

  subroutine assert_real_scalar(actual, expected, tol, msg, rc, returnmsg)
    real(real_kind),   intent(in)  :: actual, expected, tol
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = (abs(actual - expected) <= tol)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       write(returnmsg, '(2(a,g15.8))') "Fail: " // trim(msg) // " | Expected ", expected, ", got ", actual
    end if
  end subroutine assert_real_scalar

  subroutine assert_double_scalar(actual, expected, tol, msg, rc, returnmsg)
    real(dbl_kind),    intent(in)  :: actual, expected, tol
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = (abs(actual - expected) <= tol)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       write(returnmsg, '(2(a,g20.13))') "Fail: " // trim(msg) // " | Expected ", expected, ", got ", actual
    end if
  end subroutine assert_double_scalar

  subroutine assert_double_1d(actual, expected, tol, msg, rc, returnmsg)
    real(dbl_kind),    intent(in)  :: actual(:), expected(:), tol
    character(len=*),  intent(in)  :: msg
    logical,           intent(out) :: rc
    character(len=*),  intent(out) :: returnmsg

    rc = all(abs(actual - expected) <= tol)
    if (rc) then
       returnmsg = "Pass: " // trim(msg)
    else
       returnmsg = "Fail: " // trim(msg) // " | At least one element mismatched."
    end if
  end subroutine assert_double_1d
end module assertion_mod
