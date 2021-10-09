      Program DDNO
      implicit none
      character(len=512) :: filename_NIO
      integer :: whole_DDNO, IUnit
      real :: total_DDNO, fractional_DDNO
      real :: total_DDNO_Alpha, total_DDNO_Beta 

      write(*,*) "Finding total DDNO eigen value."
      call Get_Command_Argument(1,filename_NIO)
      write(*,*) "Calling from file ", filename_NIO
      open(File=trim(filename_NIO), Unit = IUnit)
    
      call NIO_Eigenvalue(filename_NIO,total_DDNO,total_DDNO_Alpha, &
        total_DDNO_Beta,IUnit)

      total_DDNO = total_DDNO/2
      total_DDNO_Alpha = total_DDNO_Alpha/2
      total_DDNO_Beta = total_DDNO_Beta/2

!     Print out whole DDNo
      write(*,*) "Whole DDNO total.", whole_DDNO 
!     Print out fractional eigen values
      write(*,*) "Fractional DDNO total.", fractional_DDNO 
!     Print out total eigen values
      
      write(*,*) "Total DDNO for alphas.", total_DDNO_Alpha 
      write(*,*) "Total DDNO for betas.", total_DDNO_Beta
      write(*,*) "Total DDNO.", total_DDNO 

      End Program DDNO

      Subroutine NIO_Eigenvalue(filename,total_eig,total_eig_alpha, & 
          total_eig_beta,IUnit)

!     This subroutine calls the string in question (in this case NIO_Eigenvalue)
!     You can use this subroutine and modify it for any string
!     and read the values at the end of the string.     

      implicit none 
      character(len = 512) :: str1, str2, str3, str4, checkstring, checkstring2
      character(len = 512) :: str5, str6
      character,intent(in) :: filename
      logical :: found = .false.
      real,intent(out):: total_eig,total_eig_alpha,total_eig_beta
      real :: temp_total
      integer :: io, IUnit

1000  format (A61)
2000  format (F10.2)

!     str1 is the whole line by line for the first 61 columns
!     You can change the format statement above depending on how much
!     you want to read into the line.

!     str2 is the check string against str1. Change to whatever
!     phrase you want. Make sure to index through correct values!

      str2 = ('Alpha NIOs:')
      str4 = ('Beta NIOs:')

!     Read through Alpha NIOs
!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:11) == str2) then 
          found = .true.
          write(*,*) "Found Alpha NIOS:"
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found Alpha NIOS:'
          stop
        end if
      end do

      read(IUnit,1000) 
      read(IUnit,1000) 

!     Makes loop through NIO Eigen values until encounters empty line
!     Totals all values back into total_eig

      do while(checkstring.ne."")
        read(IUnit,1000) checkstring 
        read(checkstring(28:32),2000) temp_total
        temp_total = ABS(temp_total)
        total_eig_alpha = total_eig_alpha + temp_total
      end do

      found =.false.

!     Read through Beta NIO
!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str3
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str3(1:10) == str4) then 
          found = .true.
          write(*,*) "Found Beta NIOS:"
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found Beta NIOS:'
          stop
        end if
      end do

      read(IUnit,1000)  
      read(IUnit,1000) 

!     Makes loop through NIO Eigen values until encounters empty line
!     Totals all values back into total_eig

      do while(checkstring2.ne."")
        read(IUnit,1000) checkstring2 
        read(checkstring2(28:32),2000) temp_total
        temp_total = ABS(temp_total)
        total_eig_beta = total_eig_beta + temp_total
      end do

      total_eig = total_eig_alpha + total_eig_beta
!     do while(
!     read(str1(59:61),'(I3)') 

      found =.false.

      End Subroutine NIO_Eigenvalue
