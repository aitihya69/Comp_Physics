subroutine nearest_neigh(init_pos,n,i,j)
    integer::n,init_pos(n,n),nn1,nn2,nn3,nn4,i,j,a,b,c,d !init_pos is a 2d array containing the lattice points
    !i and j are the index of a single point and the nearest neighbours to this points are calculated.
    !here I've used a square lattice with n number of points in both direction.
    if((i+1)<=n)then
        nn1=init_pos(i+1,j)
    else
        nn1=init_pos(1,j)
    endif
    if((i-1)>=1)then
        nn2=init_pos(i-1,j)
    else
        nn2=init_pos(n,j)
    endif
    if((j+1)<=n)then
        nn3=init_pos(i,j+1)
    else
        nn3=init_pos(1,j)
    endif
    if((j-1)>=1)then
        nn4=init_pos(i,j-1)
    else
        nn4=init_pos(n,j)
    endif


    write(*,*)"nearest neighbours","of",init_pos(i,j),"are: ",nn1,nn2,nn3,nn4
    end subroutine
