!
!     This program computes the excitation number from .fchk files 
!     Program was written and based on Perter Gills' Paper: 
!     "Excitation Number: Characterizing Mutiply Excited States"
!     DOI: 10.1021/acs.jctc.7B00963
!
!     Written by Andrew Bovill August 2021


      Program Excitation
      implicit none 
      character(len=512) :: filename_1
      integer :: N_Electrons,N_Basis,N_Basis2
      real :: Excitation_Number
      real, dimension(:,:):: Alpha_MO_Array, Beta_MO_Array
      integer,parameter :: IUnit=10

!     Obtain number of electrons and number of basis functions.

      write(*,*) "What's the fchk name?" 
      call Get_Command_Argument(1,filename_1)
      write(*,*) filename_1
      open(File=trim(filename_1), Unit = IUnit)
      call Number_of_Electrons(filename_1,N_Electrons)
      call Number_of_Basis_Functions(filename_1,N_Basis)

      call MO_Matrix(filename_1,1,Alpha_MO_Array)
      call MO_Matrix(filename_1,2,Beta_MO_Array)
      
      write(*,*) 'Number of Electrons = ', N_Electrons  
      write(*,*) 'Number of Basis = ', N_Basis  

!     Close file.
      close (IUnit, status='keep') 

      End Program Excitation


      Subroutine Number_of_Electrons(filename1,NElectrons)

!     This subroutine calls the string in question (in this case number
!     of electrons) You can use this subroutine and modify it for any string
!     and read the values at the end of the string.     



      implicit none 
      character(len = 512) :: str1, str2
      character,intent(in) :: filename1
      logical :: found = .false.
      integer,intent(out):: NElectrons
      integer:: io
      integer,parameter :: IUnit=10

1000  format (A61)

!     str1 is the whole line by line for the first 61 columns
!     You can change the format statement above depending on how much
!     you want to read into the line.

!     str2 is the check string against str1. Change to whatever
!     phrase you want. Make sure to index through correct values!

      str2 = ('Number of electrons')

!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:19) == str2) then 
          found = .true.
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found # of Electrons.'
          stop
        end if

      end do
!     Convert end of str1 into integer value for Nelectrons
!     Change index above if your column in the file of the line
!     is different in your program.
      read(str1(59:61),'(I3)') Nelectrons

!     Close file.
!     close (IUnit, status='keep') 

      End Subroutine Number_of_Electrons

      Subroutine Number_of_Basis_Functions(filename1,NBasis)

!     This subroutine calls the string in question (in this case number
!     of basis sets) You can use this subroutine and modify it for any 
!     string and read the values at the end of the string.     


      implicit none 
      character(len = 512) :: str1, str2
      character,intent(in) :: filename1
      logical :: found = .false.
      integer,intent(out):: NBasis
      integer:: io
      integer,parameter :: IUnit=10

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

!     Close file.
!     close (IUnit, status='keep') 

      End Subroutine Number_of_Basis_Functions




      Subroutine MO_Matrix(filename1,Spin,MO_Array)

!     This subroutine R the string in question has three components
!     1. Read in string to find MO coefficients for Alpha or Beta
!     2. Read in MO coefficients as an array of NBasis*NBasis length
!     3. Sends MO coefficient matrix back as array

      implicit none 
      character(len = 512) :: str1, str2, str3
      character,intent(in) :: filename1
      logical :: found = .false.
      integer, intent(in) :: Spin
      real,dimension(:,:), intent(out) :: MO_Array  
      integer:: io
      integer,parameter :: IUnit=10

1000  format (A61)

      str2 = ('Alpha MO coefficients')
      str3 = ('Beta MO coefficients')

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
          read (IUnit,*,iostat=io) MO_Array
        end if

        if(str1(1:20) == str3.and.Spin == 2) then 
          found = .true.
          write(*,*) 'Found Beta MO Coefficients'
          !Read into Beta Array
          read (IUnit,*,iostat=io) MO_Array
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found MO Coefficients'
          stop
        end if

      end do

      End Subroutine MO_Matrix
