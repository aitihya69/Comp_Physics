program main
    integer::seed_size,n_mc,n=10!(size of the grid)
    real,dimension(:,:),allocatable::init_pos,energy
    real::epsa,Temp,net_spin
    integer,dimension(:),allocatable::new_seed
    integer,dimension(1:8)::d_seed
    call RANDOM_SEED(size=seed_size)
    allocate(new_seed(seed_size))
    call ran_seed(seed_size,d_seed,new_seed)
    print*,"type the number of steps: "
    read(*,*)n_mc
    print*,"Type the temperature of the bath: "
    read(*,*)Temp
    open(3,file="ising_net_energy.dat",status="unknown")
    open(4,file="ising_net_spin.dat",status="unknown")
    allocate(init_pos(n,n),energy(n,n))
    call init_cond(init_pos,n)
    call energy_calc(init_pos,n,energy,epsa)
    print*,epsa/(n*n)
    print*,(sum(init_pos))/(n*n)
    net_spin=0.0
    do i=1,n
        do j=1,n
            net_spin=net_spin+init_pos(i,j)
        enddo
    enddo
    do i=1,n_mc
        write(3,*)i,epsa/(n**2)
        write(4,*)i,net_spin/(n**2)
        net_spin=0.0
        call metropolis(n,init_pos,N_mC,Temp,epsa,net_spin)
    enddo

end program main

subroutine init_cond(init_pos,n)
    integer::n
    real::init_pos(n,n)
    open(1,file="init_pos.dat",status="unknown")
    do i=1,n
        do j=1,n
            call random_number(init_pos(i,j))
            if (init_pos(i,j)>0.50)then
                init_pos(i,j)=1.0
            else
               init_pos(i,j)=-1.0
           endif
        enddo
    enddo
    do i=1,n
        write(1,*)(init_pos(i,j),j=1,n)
    enddo
end subroutine

subroutine energy_calc(init_pos,n,energy,epsa)
    integer::n
    real::init_pos(n,n),energy(n,n)
    real::epsa,nn1,nn2,nn3,nn4
    open(2,file="ising_energy.dat",status="unknown")

    do i=1,n
        do j=1,n
            call nearest_neigh(init_pos,n,i,j,nn1,nn2,nn3,nn4)
            energy(i,j)=-init_pos(i,j)*(nn1+nn2+nn3+nn4) 
        enddo
    enddo
    do i=1,n
        write(2,*)(energy(i,j),j=1,n)
    enddo
    epsa=sum(energy)/2.0



end subroutine

subroutine nearest_neigh(init_pos,n,i,j,nn1,nn2,nn3,nn4)
    integer,intent(in)::n,i,j
    real,intent(in)::init_pos(n,n)
    real,intent(out)::nn1,nn2,nn3,nn4
    if((i-1)>=1)then
        nn1=init_pos(i-1,j)
    else
        nn1=init_pos(n,j)
    endif
    if((i+1)<=n)then
        nn2=init_pos(i+1,j)
    else
        nn2=init_pos(1,j)
    endif
    if((j-1)>=1)then
        nn3=init_pos(i,j-1)
    else
        nn3=init_pos(i,n)
    endif
    if((j+1)<=n)then
        nn4=init_pos(i,j+1)
    else
        nn4=init_pos(i,1)
    endif
    end subroutine

subroutine metropolis(n,init_pos,N_mC,Temp,epsa,net_spin)
    integer::n,n_mc,x1,y1
    real::de,init_pos(n,n),net_spin,temp,epsa,r,spin_i,spin_f,i_epsa,f_epsa,nn1,nn2,nn3,nn4
    do i=1,n
        do j=1,n
            call random_number(r)
            x1=int(r*n)
            call random_number(r)
            y1=int(r*n)
            spin_i=init_pos(x1,y1)
            spin_f=spin_i*(-1.0)
            call nearest_neigh(init_pos,n,x1,y1,nn1,nn2,nn3,nn4)
            i_epsa=-spin_i*(nn1+nn2+nn3+nn4)
            f_epsa=-spin_f*(nn1+nn2+nn3+nn4)
            de=(f_epsa-i_epsa)
            if(de>0.0)then
                call random_number(r)
                if(r<exp(-(de)/(Temp)))then
                    init_pos(x1,y1)=spin_f
                    epsa=epsa+de
                else
                    spin_f=-spin_i
                endif
            elseif(de<=0.0)then
                init_pos(x1,y1)=spin_f
                epsa=epsa+de
                endif
            enddo
        enddo
        net_spin=sum(init_pos)
    end subroutine



    subroutine ran_seed(seed_size,d_seed,new_seed)
        integer::seed_size,d_seed(8),new_seed(seed_size)
        call RANDOM_SEED(get=new_seed)
        call DATE_AND_TIME(values=d_seed)
        do i=1,8
            new_seed(i)=d_seed(8)*d_seed(7)
        enddo
        print*,new_seed(1)
        call RANDOM_SEED(put=new_seed)
    end subroutine
