!
!     This program computes the excitation number from .fchk files 
!     Excitation number was proposed by Peter M. W. Gill to 
!     represent the number of electrons in the excited state that lie 
!     in the unoccupied space of the ground state
!     "Excitation Number: Characterizing Mutiply Excited States"
!     DOI: 10.1021/acs.jctc.7B00963
!
!     Author Andrew Bovill 
!     August 2021
!
!     Program is compiled with the LAPACK library
!
!
!     Program uses "LAPACK_Helper" to invert a GE Matrix 
!     You can find the source code here.
!     https://github.com/b-fg/LAPACK_helper/blob/master/lapack_helper.f90
!

      Program Excitation

      use lapackmod 

      implicit none 
!     First file  filename_GS is Ground state
!     Second file filename_EX is Excited state
      character(len=512) :: filename_GS, filename_EX
      integer :: N_Electrons_GS, N_Electrons_GS_Alpha, N_Electrons_GS_Beta 
      integer :: N_Electrons_EX, N_Electrons_EX_Alpha, N_Electrons_EX_Beta 
      integer :: N_Basis_GS, N_Basis_EX
      integer :: i, j
      real :: Excitation_Number, Alpha_Excitation, Beta_Excitation

!     Arrays
      real, allocatable, dimension(:) :: Alpha_MO_Array_GS, Beta_MO_Array_GS
      real, allocatable, dimension(:) :: Alpha_MO_Array_EX, Beta_MO_Array_EX

!     Matrices
      real, allocatable, dimension(:,:) :: Alpha_MO_Matrix_GS, Beta_MO_Matrix_GS
      real, allocatable, dimension(:,:) :: Alpha_MO_Matrix_EX, Beta_MO_Matrix_EX

!     Matrices of occupied orbitals.
      real, allocatable, dimension(:,:) :: Alpha_OCC_AB, Beta_OCC_AB

!     Overlap matrix (identity matrix is just test)
      real, allocatable, dimension(:,:) :: identity_matrix
      real, allocatable, dimension(:,:) :: Overlap_Alpha
      real, allocatable, dimension(:,:) :: Overlap_Beta

!     Integers for write files.
      integer,parameter :: IUnit_GS=10
      integer,parameter :: IUnit_EX=20

!     Read in Electrons and Basis number for Ground State
      write(*,*)
      write(*,*) "What's the Ground State fchk file?" 
      call Get_Command_Argument(1,filename_GS)
      write(*,*) filename_GS
      open(File=trim(filename_GS), Unit = IUnit_GS)
      call Number_of_Electrons(filename_GS,N_Electrons_GS_Alpha,IUnit_GS,1)
      call Number_of_Electrons(filename_GS,N_Electrons_GS_Beta,IUnit_GS,2)
      call Number_of_Basis_Functions(filename_GS,N_Basis_GS,IUnit_GS)

!
      N_Electrons_GS = N_Electrons_GS_Alpha + N_Electrons_GS_Beta
      write(*,*) 'Number of Alpha Electrons in Ground State = ', N_Electrons_GS_Alpha  
      write(*,*) 'Number of Beta Electrons in Ground State = ', N_Electrons_GS_Beta  
      write(*,*) 'Number of Electrons in Ground State = ', N_Electrons_GS  
      write(*,*) 'Number of Basis in Ground State = ', N_Basis_GS  

!     Read into Arrays for both Alpha and Beta Ground States
      Allocate (Alpha_MO_Array_GS(N_Basis_GS*N_Basis_GS))
      call GET_MO_Array(IUnit_GS,N_Basis_GS,Alpha_MO_Array_GS,1)
      Alpha_MO_Matrix_GS = reshape(Alpha_MO_Array_GS, (/ N_Basis_GS,N_Basis_GS /))

      Allocate (Beta_MO_Array_GS(N_Basis_GS*N_Basis_GS))
      call GET_MO_Array(IUnit_GS,N_Basis_GS,Beta_MO_Array_GS,2)
      Beta_MO_Matrix_GS = reshape(Beta_MO_Array_GS, (/ N_Basis_GS,N_Basis_GS /))

      close (IUnit_GS, status='keep') 

!     Read in Electrons and Basis number for Excited State
      write(*,*)
      write(*,*) "What's the Excited State fchk file?" 
      call Get_Command_Argument(2,filename_EX)
      write(*,*) filename_EX
      open(File=trim(filename_EX), Unit = IUnit_EX)
      call Number_of_Electrons(filename_EX,N_Electrons_EX_Alpha,IUnit_EX,1)
      call Number_of_Electrons(filename_EX,N_Electrons_EX_Beta,IUnit_EX,2)
      call Number_of_Basis_Functions(filename_EX,N_Basis_EX,IUnit_EX)

      N_Electrons_EX = N_Electrons_EX_Alpha + N_Electrons_EX_Beta
      write(*,*) 'Number of Alpha Electrons in Excited State = ', N_Electrons_EX_Alpha  
      write(*,*) 'Number of Beta Electrons in Excited State = ', N_Electrons_EX_Beta  
      write(*,*) 'Number of Electrons in Excited State = ', N_Electrons_EX  
      write(*,*) 'Number of Basis in Excited State = ', N_Basis_EX  

!     Read into Arrays for both Alpha and Beta Excited States
      Allocate (Alpha_MO_Array_EX(N_Basis_EX*N_Basis_EX))
      call GET_MO_Array(IUnit_EX,N_Basis_EX,Alpha_MO_Array_EX,1)
      Alpha_MO_Matrix_EX = reshape(Alpha_MO_Array_EX, (/ N_Basis_EX,N_Basis_EX /))

      Allocate (Beta_MO_Array_EX(N_Basis_EX*N_Basis_EX))
      call GET_MO_Array(IUnit_EX,N_Basis_EX,Beta_MO_Array_EX,2)
      Beta_MO_Matrix_EX = reshape(Beta_MO_Array_EX, (/ N_Basis_EX,N_Basis_EX /))

      close (IUnit_EX, status='keep') 
     
      if(N_Electrons_GS.ne.N_Electrons_EX) then
        write(*,*) "The number of the electrons in both files do not match!"
        STOP
      end if

      if(N_Basis_GS.ne.N_Basis_EX) then
        write(*,*) "The number of the basis in both files do not match!"
        STOP
      end if
!
!     Generate Overlap for Alpha & Beta matrices
!     

      Overlap_Alpha = matmul(inv(transpose(Alpha_MO_Matrix_GS)),inv(Alpha_MO_Matrix_GS))
      Overlap_Beta = matmul(inv(transpose(Beta_MO_Matrix_GS)),inv(Beta_MO_Matrix_GS))

!
!     Compute the Alpha & Beta occupied matrix from the generated Overlap_Alpha &
!     Overlap_Beta
!
      ALLOCATE(Alpha_OCC_AB(N_Electrons_GS_Alpha,N_Electrons_EX_Alpha))
      ALLOCATE(Beta_OCC_AB(N_Electrons_GS_Beta,N_Electrons_EX_Beta))

      Alpha_OCC_AB = MatMul(Transpose(Alpha_MO_Matrix_GS(:,1:N_Electrons_GS_Alpha)),  &
         MatMul(Overlap_Alpha,Alpha_MO_Matrix_EX(:,1:N_Electrons_EX_Alpha)))

      Beta_OCC_AB = MatMul(Transpose(Beta_MO_Matrix_GS(:,1:N_Electrons_GS_Beta)),  &
         MatMul(Overlap_Beta,Beta_MO_Matrix_EX(:,1:N_Electrons_EX_Beta)))

!     Obtain number of electrons promoted for both Alpha and Beta electrons. 

      Alpha_Excitation = N_Electrons_GS_Alpha - &
      dot_product(reshape(Alpha_OCC_AB,(/N_Electrons_GS_Alpha*N_Electrons_EX_Alpha/)),  &
      reshape(Alpha_OCC_AB,(/N_Electrons_GS_Alpha*N_Electrons_EX_Alpha/)))
      
      Beta_Excitation = N_Electrons_GS_Beta - &
      dot_product(reshape(Beta_OCC_AB,(/N_Electrons_GS_Beta*N_Electrons_EX_Beta/)),  &
      reshape(Beta_OCC_AB,(/N_Electrons_GS_Beta*N_Electrons_EX_Beta/)))

      Excitation_Number =  Alpha_Excitation + Beta_Excitation

2000  Format(/,1x,'Excitation Number: ',4x,F6.4/, &
        1x,'Alpha Contribution:',4x,F6.4,/, &
        1x,'Beta Contribution:',5x,F6.4)

      write(*,2000) Excitation_Number,Alpha_Excitation,Beta_Excitation

      End Program Excitation

      Subroutine Number_of_Electrons(filename,NElectrons,IUnit,Spin)

!     This subroutine reads the input file line by line for the number of 
!     Alpha and Beta electrons

      implicit none 
      character(len = 512) :: str1, str2, str3
      character,intent(in) :: filename
      integer, intent(in) :: Spin
      logical :: found = .false.
      integer,intent(out):: NElectrons
      integer :: io, IUnit

1000  format (A61)

!     str1 is the whole line by line for the first 61 columns
!     You can change the format statement above depending on how much
!     you want to read into the line.

!     str2 is the check string against str1. Change to whatever
!     phrase you want. Make sure to index through correct values!

      str2 = ('Number of alpha electrons')
      str3 = ('Number of beta electrons')

!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:25) == str2.and.Spin == 1) then 
          found = .true.
          write(*,*) 'Found Alpha Electrons'
        end if

        if(str1(1:24) == str3.and.Spin == 2) then 
          found = .true.
          write(*,*) 'Found Beta Electrons'
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found Electrons'
          stop
        end if

      end do
!     Convert end of str1 into integer value for Nelectrons
!     Change index above if your column in the file of the line
!     is different in your program.
      read(str1(59:61),'(I3)') NElectrons

      found =.false.
      write(*,*)
      End Subroutine Number_of_Electrons

      Subroutine Number_of_Basis_Functions(filename1,NBasis,IUnit)

!     This subroutine calls the string in question (in this case number
!     of basis sets) You can use this subroutine and modify it for any 
!     string and read the values at the end of the string.     


      implicit none 
      character(len = 512) :: str1, str2
      character,intent(in) :: filename1
      logical :: found = .false.
      integer,intent(out):: NBasis
      integer:: io, IUnit

1000  format (A61)

!     str1 is the whole line by line for the first 61 columns
!     You can change the format statement above depending on how much
!     you want to read into the line.

!     str2 is the check string against str1. Change to whatever
!     phrase you want. Make sure to index through correct values!

      str2 = ('Number of basis functions')

!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:25) == str2) then 
          found = .true.
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found # of Basis Functions.'
          stop
        end if

      end do
!     Convert end of str1 into integer value for Nelectrons
!     Change index above if your column in the file of the line
!     is different in your program.
      read(str1(59:61),'(I3)') NBasis

      found =.false.
      
      End Subroutine Number_of_Basis_Functions

      Subroutine GET_MO_Array(IUnit,NBasis,MO_Array,Spin)

!     This subroutine does the following:
!     1. Read in string to find MO coefficients for Alpha or Beta
!     2. Read in MO coefficients as an array of NBasis*NBasis length
!     3. Sends MO coefficient matrix back as array

      
      implicit none
      character(len = 512) :: str1, str2, str3
      logical :: found = .false.
      integer, intent(in) :: Spin,NBasis
      real,dimension(NBasis*NBasis), intent(out) :: MO_Array  
      integer:: io, IUnit

1000  format (A61)

      str2 = ('Alpha MO coefficients')
      str3 = ('Beta MO coefficients')

      write(*,*)
!     Resets logical for multiple subroutine runs
      found = .false.

!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:21) == str2.and.Spin == 1) then 
          found = .true.
          write(*,*) 'Found Alpha MO Coefficients'
          !Read into Alpha Array
          read (IUnit,*) MO_Array
        end if

        if(str1(1:20) == str3.and.Spin == 2) then 
          found = .true.
          write(*,*) 'Found Beta MO Coefficients'
          !Read into Beta Array
          read (IUnit,*) MO_Array
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found MO Coefficients'
          stop
        end if

      end do

      End Subroutine GET_MO_Array
