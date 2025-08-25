program check
 implicit none
 integer,dimension(:,:), allocatable :: permuted_list !for generate_permutation
 integer,dimension(:,:), allocatable :: bond_list !for bond indices (contain permuted one)
 integer,dimension(:,:), allocatable :: bond_indices_list !same as bond indices (contain original/starting bond indices)
 integer,dimension(:), allocatable :: perm_sub_list !for sub list of permuted list (containing data for single permutation)
 integer,dimension(:), allocatable :: induced_list !for induced permuted list 
 
 !used in final list code
 integer,dimension(:,:), allocatable :: final_list !list of list which containing all induced list
 integer::mlo,mro,f
 integer::a,b,c,d,e
 integer::s,av
 integer::spl
 integer::starting_no
 
 !used in singular input file
 integer,dimension(:,:,:), allocatable :: singular_list
 integer,dimension(:,:), allocatable :: singular_matrix_list
 
!used in printing input file
 character(len=2),dimension(:), allocatable :: c_singular_matrix_list !singular_matrix_list in form of character with commas
 character(len=2)::chr


!------INPUT YOUR VALUE HERE-----------
 integer,dimension(2) :: atom_vector  !list of atom vector. Atom vector is noting but way to represent atom group 
 atom_vector=(/2,2/)               !like A2B3C4 = atom vector= (/2,3,4/) and dimension is total number of elements i.e size of list


!--------------------------------------------------------------------------

!for calulating order of final list
!final list is in form of 2d array(matrix) 
!mlo=matrix left order (this is equal to total no of permutations possible neglecting original configuration)
!mro=matrix right order (this is equal to total number of bond interactions)

mlo=0
do a=1,size(atom_vector)
 f=1
 do b=1,atom_vector(a)
   f=f*b
 end do 
 mlo=mlo+f-1
end do

f=sum(atom_vector)
mro=((f*f)-f)/2


!by using all subroutine finding final list
allocate(final_list(mlo,mro))
e=0
av=0
call bond_indices(sum(atom_vector),bond_indices_list)

starting_no=1
do a=1,size(atom_vector)
    
    if (atom_vector(a)==1) then
    cycle
    end if
    
    call generate_permutation(atom_vector(a),permuted_list,starting_no) !heree starting no is number after which we want permutations
    
    do b=2,size(permuted_list,1)
        allocate (perm_sub_list(size(permuted_list,2)))
        perm_sub_list=0
        
        !for generating permuted sub list
        do c=1,size(permuted_list,2)
            perm_sub_list(c)=permuted_list(b,c)
        end do
        
        call bond_indices(sum(atom_vector),bond_list)
        call permuted_bond_indices(bond_list,perm_sub_list)
        call induced_permutation(bond_indices_list,bond_list,induced_list)
        
        !inserting induced_list in final list
        e=e+1
        do d=1,mro
            final_list(e,d)=induced_list(d)
        end do
        
        deallocate (perm_sub_list)
        deallocate(induced_list)
        deallocate(bond_list)
    end do    
    deallocate (permuted_list)
    starting_no=starting_no+atom_vector(a)
end do


!call l_representation(final_list)

! for finding 3-matrix containing all matrix(they are baisically position of bonds same as final list but in matix form)  
allocate(singular_list(mlo,mro,mro))
singular_list=0
do a=1,mlo
    do b=1,mro
     singular_list(a,b,final_list(a,b))=1
    end do
end do

!--------------------------------------------------------------------------

!all rest part is for printing in same format like singular input file(Can be done in better way)


write(*,'(A)')'LIB "finvar.lib";'
write(*,'(A)' ,advance="no")'ring R=0,('
do a=1,mro
    write (*, '(A,i0,A)' ,advance="no")"x",a
    if (a==mro) then
        write(*,'(A)')'),dp;'
    else
        write(*,'(A)' ,advance="no") ','
    end if
end do


!line 3
do a=1,mlo
    allocate(c_singular_matrix_list(mro*mro))!its 1-d list 
    d=1
    do b=1,mro
        do c=1,mro
          write(chr,'(I1)')singular_list(a,b,c)
          if (d==mro*mro) then
          c_singular_matrix_list(d)=trim(chr)//';'
          else
          c_singular_matrix_list(d)=trim(chr)//',' 
          end if
          d=d+1
        end do
    end do
    
    write(*,'(A,i0,A,i0,A,i0,A)',advance="no")'matrix A',a,'[',mro,'][',mro,']='
    write(*,*)c_singular_matrix_list
    deallocate(c_singular_matrix_list)
end do



!line 4,5

write(*,'(A)' ,advance="no")'list L = group_reynolds('
do a=1,mlo
    write (*, '(A,i0,A)' ,advance="no")"A",a
    if (a==mlo) then
        write(*,'(A)')');'
    else
        write(*,'(A)' ,advance="no") ','
    end if
end do

write(*,'(A)')'matrix G = invariant_algebra_reynolds(L[1],1); G;'

!--------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------

!generate permuted list in form of 2D-list
!e.g for A3 it stores[[1,2,3],[1,3,2],[2,1,3],[2,3,1],[3,1,2],[3,2,1]]
subroutine generate_permutation(n,permuted_list,starting_no) 
  implicit none
  integer :: n,starting_no
  integer, dimension(:), allocatable :: perm, used
  
  integer:: m,l  !for permuted perm6uted_list
  integer,dimension(:,:), allocatable :: permuted_list

  m=1
  do l=1,n
  m=m*l
  end do
  allocate(permuted_list(m,n))
  permuted_list=0
  ! Allocate memory for arrays based on n
  allocate(perm(n), used(n))
  
  ! Initialize used array to zero (not used)
  used = 0
  perm=0
  ! Call recursive function to generate permutations
  call permute(perm, used, 1, n,starting_no)
  ! Deallocate memory (optional but recommended)
  deallocate(perm, used)
end subroutine generate_permutation 



!part of above subroutine
recursive subroutine permute(perm, used, current, n,starting_no) !used in generation of permuted list
    implicit none
    integer:: perm(:), used(:)
    integer :: current, n,i,starting_no
    !for permuted_list
    integer :: mk,k=1
    ! Base case: permutation complete (all positions filled)
    
    if (current .eq. n + 1) then
    do mk=1,n
        permuted_list(k,mk)=perm(mk)
    end do
    k=k+1
    return
    endif

    ! Loop through unused numbers
    do i = 1,n
        if (used(i) .eq. 0) then
        ! Mark number as used
        used(i) = 1
        ! Place number l_representationin current position of permutation
        perm(current) = i+starting_no-1
        ! Recursively generate permutations with remaining positions
        call permute(perm, used, current + 1, n,starting_no)
        ! Unmark number after recursion
        used(i) = 0
        endif
    enddo
    
    if (permuted_list(ubound(permuted_list,1),n)/=0) then
        k=1
    end if
end subroutine permute

!--------------------------------------------------------------------------

!2d array representing bond indices
!The number of atoms
    ! Finds the array of bond indices of the upper triangle of an interatomic distance matrix, in column wise order
    ! ( or equivalently, lower triangle of interatomic distance matrix in row wise order):
    ! [ [1,2], [1,3], [2,3], ...,[1, natom], ...,[natom-1, natom]]
    
subroutine bond_indices(total_atom,bond_list)
 implicit none
 integer,dimension(:,:), allocatable :: bond_list
 integer::n,m,total_atom
 integer::s,i,j
 n=total_atom
 m=((n*n)-n)/2
 
 allocate(bond_list(m,2))
  s=m
  j=n
 do while (j>1)
  i=j-1
  do while (i>0)
    bond_list(s,1)=i
    bond_list(s,2)=j
    s=s-1
    i=i-1
    !print*,i,j,s
  end do
  j=j-1
 end do
end subroutine bond_indices
        
!--------------------------------------------------------------------------

!apply a single permutation on bond indices list
!Permutes a bond inidice if the bond indice is affected by the permutation cycle.
subroutine permuted_bond_indices(bond_list,perm_sub_list)
    implicit none
    integer,dimension(:), allocatable :: perm_sub_list 
    integer,dimension(:,:), allocatable :: bond_list 
    integer::count1,count2,var1
    integer::i,j    
    integer::p

    do i=1,size(bond_list,1)
      count1=0
      count2=0 
      var1=0

        do j=minval(perm_sub_list),minval(perm_sub_list)+size(perm_sub_list)-1
        var1=var1+1
        if (j==perm_sub_list(j)) then
          cycle
      
        else if (bond_list(i,1)==j .and. count1==0) then
          bond_list(i,1)=perm_sub_list(var1)
          count1=count1+1
          
        else if (bond_list(i,2)==j .and. count2==0) then
          bond_list(i,2)=perm_sub_list(var1)
          count2=count2+1
        
        end if
        
      end do
        if (bond_list(i,1)>bond_list(i,2)) then
          p=bond_list(i,1)
          bond_list(i,1)=bond_list(i,2)
          bond_list(i,2)=p
        end if      
    end do    
end subroutine permuted_bond_indices

!--------------------------------------------------------------------------

!convert a given 2d permuted bond indices list into 1d permuted list according to position of bond by comparing with original bond indices 
subroutine induced_permutation(bond_indices_list,bond_list,induced_list)!bond indices_list=original ;bond_list=permuted
 
 integer,dimension(:,:), allocatable :: bond_indices_list 
 integer,dimension(:,:), allocatable :: bond_list 
 integer,dimension(:), allocatable :: induced_list
 integer::i,j
 
 allocate(induced_list(size(bond_list,1)))
 do i=1,size(bond_list,1)
  do j=1,size(bond_indices_list,1)
   if ((bond_list(i,1)==bond_indices_list(j,1)).and.(bond_list(i,2)==bond_indices_list(j,2))) then
     induced_list(i)=j
     exit
   end if
  end do
 end do
end subroutine induced_permutation


!for representing 2D-list (just if want to see any 2-d array)       
subroutine l_representation(data_list) 
  implicit none
  integer, dimension(:,:) :: data_list
  integer::i,j
  do i = 1, size(data_list, 1)  ! Loop through rows (sublists)
    write (*, *) "Sublist", i, ":"
    do j = 1, size(data_list, 2)  ! Loop through elements in each sublist
      write (*, *) "  Element", j, ":", data_list(i, j)
    end do
  end do
  return
end subroutine l_representation        


end program
