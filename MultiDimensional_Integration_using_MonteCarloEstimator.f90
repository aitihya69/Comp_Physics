program mult_dim_mc

!---THINGS TO CONSIDER:  FIRST OF ALL WE NEED TO TYPE FUNCTION THAT WE WANNA INTEGRATE, AND ACCORDING TO THE DIMENSION OF THE FUNCTION
    !---WE NEED  TO CHANGE THE FUNCTION PART OF THE PROGRAM BY ADDING CORRECT VARIABLES AND STUFF
        !---AND WHILE CALLING THE FUNCTION TO CALCULATE THE VALUE OF SUM WE NEED TO CHANGE THE ARGUMENTS OF THE FUNCTION ACCORDING TO OUR NEEDS.


    real*8,dimension(:),allocatable::x,sum1
    real*8,dimension(:,:),allocatable::limits
    real*8::result,m_vol,rand,interval,x1,f
    integer::m,n
    print*,"Type the dimension of the integration: "
    read(*,*)m
    allocate(x(m),limits(m,2),sum1(1000))
    print*,"Type the number of monte carlo samples: "
    read(*,*)n

    !---getting the limits of the integration in each dimension---!
    print*,"Type the lower and upper limits of the variables of the integration: "
    do i=1,m
        do j=1,2
            read(*,*)limits(i,j)
        enddo
    enddo

    m_vol=1.0d0

    do i=1,m
        interval=abs((limits(i,2)-limits(i,1)))
        write(*,*)interval
        m_vol=m_vol*interval
    enddo

    !---Monte Carlo estimator---!
    do i=1,1000
        sum1(i)=0.0d0
    enddo
do k=1,1000
    do i=1,n
        do j=1,m
            call random_number(rand)
            x(j)=limits(j,1)+rand+(rand*(limits(j,2)-limits(j,1)-1))
        enddo
        sum1(k)=sum1(k)+f(x(1),x(2))
    enddo
enddo

    result=m_vol*((sum(sum1)/(1000))/real(n))


    write(*,*)"The value of  the integration is: ",result,m_vol
    

end program mult_dim_mc

function f(x1,x2)
    real*8::x1,x2,f
    f=(x1**2+x2**2)
end function
