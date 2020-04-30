subroutine hello_world(number)
        real, intent(in) :: number

        print*, "Main Program", number
        call write_digit(number)

end subroutine hello_world