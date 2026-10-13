PROGRAM NEB_calculation
!********************************************************
!!! Objective: To find out hopping/migration barrier from structure
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!Input: Initial and final occupation vector
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!NEB Calculation by Molecular Dynamics
!*********************************************************
  IMPLICIT NONE
   INTEGER:: N,I,J,K,NATOMS,L,NtotI,NtotF,A1,B1,C1,D1
   REAL, DIMENSION(:,:),ALLOCATABLE :: rNi,rH
   REAL, DIMENSION(:),ALLOCATABLE :: SH,SNi,OCVI,OCVF,TEMP_OCV,TEMP_OCVI
   REAL, DIMENSION(:,:),ALLOCATABLE:: coord,coordd
   REAL, ALLOCATABLE:: XLO,YLO,ZLO,XHI,YHI,ZHI
   CHARACTER(LEN=20), DIMENSION(:),ALLOCATABLE::Species,speciess
   !INTEGER :: n
   CHARACTER(LEN=32) :: initial_occupation_file,final_occupation_file
   REAL :: a,b,c
   L=10
   a=3.52e0 ! Lattice constant of Ni = 3.52 FOR SMALL H, 50%H Ni=3.66 AND 99.9%H Ni=3.74
   b=a/2e0
   c=0.0e0
   !!!! for small anount H addition
!~   XLO=0.0e0
!~   YLO=0.0e0
!~   ZLO=0.0e0
!~   XHI=53.0e0
!~   YHI=53.0e0
!~   ZHI=53.0e0
  
!~   XLO=-0.915e0
!~   YLO=-0.915e0
!~   ZLO=-0.915e0
!~   XHI=35.685e0
!~   YHI=35.685e0
!~   ZHI=35.685e0

!!!!FOR CHECKING LATTICE CONSTANT VALUE 25% H !!!!!
  
!~     XLO=-0.90e0
!~   YLO=-0.90e0
!~   ZLO=-0.90e0
!~   XHI=53.1e0
!~   YHI=53.1e0
!~   ZHI=53.1e0

!!!!! for 50% H addition (13650 Ni + 6825 H addition =20400)
!~   XLO=-0.915e0
!~   YLO=-0.915e0
!~   ZLO=-0.915e0
!~   XHI=53.985e0
!~   YHI=53.985e0
!~   ZHI=53.985e0
!!!!!!!
  !!!!FOR CHECKING LATTICE CONSTANT VALUE 75% H !!!!!
  
    XLO=-0.88e0
  YLO=-0.88e0
  ZLO=-0.88e0
  XHI=34.32e0
  YHI=34.32e0
  ZHI=34.32e0

  
!~     XLO=-0.925e0
!~   YLO=-0.925e0
!~   ZLO=-0.925e0
!~   XHI=54.575e0
!~   YHI=54.575e0
!~   ZHI=54.575e0
  !!!! for more than 12000 H addition
!~ XLO=-0.935e0
!~ YLO=-0.935e0
!~ ZLO=-0.935e0
!~ XHI=55.165e0
!~ YHI=55.165e0
!~ ZHI=55.165e0
  !!!!
  
  !StepA: Creation of Ni structure
  CALL Ni_structure()
  !StepB: Creation of H structure
  CALL H_structure()
  !StepC Reading Ni coordinates
  CALL Ni_coordinates()
  !StepD Reading H Coordinates
  CALL H_coordinates()
  !Step 1. Read input filenames
  CALL Initial_final_occupation()
  !Step 11 rewriting initial occupation
  CALL initial_occupation_reading()
  !step 12 rewriting final occupation
  CALL final_occupation_reading()
  !Step 2. Store the site occupation in array initial_occupation_vector
  CALL ReadOccupationVector_initial() !initial
  !Step 3. Store the site occupation in array final_occupation_vector
  CALL ReadOccupationVector_final() !final 
  !Step 4. Create LAMMPS file - initial configuration
  CALL CreateLAMMPSfileForRelaxation_initial()
  !Step 5. Create LAMMPS file - final configuration 
  CALL CreateLAMMPSfileForRelaxation_final()
  !Step 6. Relax the initial configuration
  !CALL RelaxUsingLAMMPS('lmp -in in_NIH_initial.eam') !initial
  CALL SYSTEM('mpirun -np 20 lmp -in in_NIH_initial.eam') !initial
  !Step 7. Relax the final configuration
  !CALL RelaxUsingLAMMPS('lmp -in in_NIH_initial1.eam') !final
  CALL SYSTEM('mpirun -np 20 lmp -in in_NIH_initial1.eam') !final
  !Step 8. Prepare LAMMPS files for NEB calculation
  CALL PrepareForNEB()
  !Step 9. Perform NEB
  !CALL NEB('mpirun -np 9 --oversubscribe lmp -partition 9x1 -in input_NEB > t5.txt')
  CALL SYSTEM('mpirun -np 9 lmp -partition 9x1 -in input_NEB')
  !Step 10. Extract the barrier from the LAMMPS output file
  !CALL ExtractNEBBarrier("grep -e 'EBF' -e '0.559' t5.txt | tail -1 | awk '{print $7}'")
  !CALL ExtractNEBBarrier("grep -e 'EBF' -e '0.559' t5.txt | tail -1")
!~   CALL SYSTEM("grep -e 'EBF' -e '0.559' t55.CSV | tail -1 | awk '{print $7}'")
  
CONTAINS
  
SUBROUTINE Ni_structure()
  IMPLICIT NONE
 OPEN(UNIT=2,FILE="NiBulk1.xyz",STATUS="REPLACE")
   
   WRITE(2,*)L*L*L*4
	WRITE(2,*)L*a,L*a,L*a,c,c,c

DO I=0,L-1
   DO J=0,L-1
      DO K=0,L-1
            WRITE(2,FMT=*)"Ni  ",REAL(K*a), REAL(J*a), REAL(I*a)
            WRITE(2,FMT=*)"Ni  ",REAL(K*a)+b, REAL(J*a)+b, REAL(I*a)
            WRITE(2,FMT=*)"Ni  ",REAL(K*a), REAL(J*a)+b, REAL(I*a)+b
            WRITE(2,FMT=*)"Ni  ",REAL(K*a)+b, REAL(J*a), REAL(I*a)+b
      END DO
   END DO
END DO
close(2)
END SUBROUTINE

!StepB: Creation of H structure
!CALL H_structure()
SUBROUTINE H_structure()
open(unit=2, file="h_bulk.xyz", status="replace")
write (2,*) L*L*L*4
write(2,*) L*a,L*a,L*a,c,c,c
   Do I=0,L-1
      Do J=0,L-1
         Do K=0,L-1
            write(2, FMT=*) "H   ", real(K*a)+b, real(J*a), real(I*a)
            write(2, FMT=*) "H   ", real(K*a), real(J*a)+b, real(I*a)
            write(2, FMT=*) "H   ", real(K*a), real(J*a), real(I*a)+b
            write(2, FMT=*) "H   ", real(K*a)+b, real(J*a)+b, real(I*a)+b
         end Do
      end Do
   end Do
   close(2)
END SUBROUTINE
!CALL Ni_coordinates()
SUBROUTINE Ni_coordinates()
IMPLICIT NONE
   character(len=2)::U
   open(unit=2,file="NiBulk1.xyz", status="OLD")
   open(unit=4,file="NiBulk1D.xyz", status="REPLACE")
   read(2,*) n
!~    print *,n
! remove space in xyz file
   read(2,*)
   allocate( rNi(n,3) )
   do I = 1,n
      read(2,*) U,rNi(I,:)
   end do
   close(2)
   J=1
   allocate( SNi(3*n) )
   do I = 1,n
      SNi(3*(I-1)+1)=rNi(I,1)
      SNi(3*(I-1)+2)=rNi(I,2)
      SNi(3*(I-1)+3)=rNi(I,3)
      J=J+3
   end do
   
    !print *, SNi(3*n)
    do i = 1,3*n
       !print *,SNi(i)
       WRITE(4,FMT=*)SNi(i)
    end do
    !WRITE(4,*)
CLOSE(4)
END SUBROUTINE

!CALL H_coordinates()
SUBROUTINE H_coordinates()

 implicit none
   character(len=20):: T
   open(unit=3, file="h_bulk.xyz", status="old")
   open (unit=4, file= "hydrogen_1D_array", status="replace")
   read(3,*) n
   !print *, n
   read(3,*)
   !read(3,*)
   allocate(rH(n,3))
   !allocate(x(n),y(n),z(n))
      Do I= 1,n
      read(3,*) T,rH(I,:)
      !read(3,*) T,x(i),y(i),z(i)
      end Do
      !print *,T
   close(3)
   J=1
   !allocate(x1(1*n),y1(1*n),z1(1*n))
   allocate( SH(3*n) )
         Do I= 1,n
         SH(3*(I-1)+1)=rH(I,1)
         SH(3*(I-1)+2)=rH(I,2)
         SH(3*(I-1)+3)=rH(I,3)
         !x1(1*(i-1)+1)=x(1)
         !y1(1*(i-1)+2)=y(2)
         !z1(1*(i-1)+3)=z(3)
         J=J+3
         !print *, SH(10)
         end Do
!~          !print *, SH(2)
         do i = 1,3*n
          !print *,SH(i)
          WRITE(4,*)SH(i)
         end do
         !WRITE(4,*)
         !print *, SH(3*n)
       close(4)
END SUBROUTINE
SUBROUTINE Initial_final_occupation()
INTEGER::M
  M=iargc()
  IF (M/=2) THEN
     WRITE(6,*) 'Please provide two input files'
     WRITE(6,*) ' Execute the code as ./neb_calculation.x file1 file2'
     STOP
  END IF
  
  !Step 1. Read input filenames
  CALL getarg(1, initial_occupation_file)
  CALL getarg(2, final_occupation_file)
  WRITE(6,*) 'Initial occupation file: ',TRIM(initial_occupation_file)
  WRITE(6,*) 'Final occupation file: ',TRIM(final_occupation_file)
  !STOP
  END SUBROUTINE
  
  !Step 2. Store the site occupation in array initial_occupation_vector
  !CALL ReadOccupationVector(initial_occupation_vector,initial_occupation_file) !initial 
  
  !Inital NIH with occupation_vector
!********************************************************************
SUBROUTINE initial_occupation_reading()
IMPLICIT NONE
!OPEN(UNIT=2, FILE="array1_initial.xyz", STATUS="OLD")
OPEN(UNIT=2, FILE="initial_file", STATUS="OLD")
OPEN(UNIT=3, FILE="array_initial_MODIFIED.xyz", STATUS="REPLACE")
ALLOCATE(OCVI(N))
OCVI=0
   DO I=1,N
   READ(2,*)OCVI(I)
   END DO
   NtotI=N+sum(OCVI)
   Do I=1,NtotI
   write(3,FMT=*)NtotI
   end Do
print *, NtotI,":initial occupation"
!~ DEALLOCATE(OCVI)
close(2)
CLOSE(3)

END SUBROUTINE
!~ subroutine generate_lammps_initial()
!~ IMPLICIT NONE
!~ !INTEGER :: NtotI
!~ open(unit=2, file="total.xyz", status="replace")
!~ NtotI= n+sum(ocvI)
!~    Do I=1,NtotI
!~    write(2,FMT=*)NtotI
!~    end Do
!~ print *, NtotI,":initial occupation"
!~ close(2)
!~ end subroutine
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!Final NIH with occupation_vector
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
SUBROUTINE final_occupation_reading()
IMPLICIT NONE
INTEGER :: II
!OPEN(UNIT=2, FILE="array1_final.xyz", STATUS="OLD")
OPEN(UNIT=2, FILE="final_file", STATUS="OLD")
OPEN(UNIT=3, FILE="array_final_MODIFIED.xyz", STATUS="REPLACE")
ALLOCATE(OCVF(N))
OCVF=0
   DO II=1,N
   READ(2,*)OCVF(II)
   END DO
   NtotF=N+sum(OCVF)
   Do II=1,NtotF
   write(3,FMT=*)NtotF
   end Do
print *, NtotF,":final occupation"
!~ DEALLOCATE(OCVF)
close(2)
CLOSE(3)
!~ subroutine generate_lammps_final()
!~ IMPLICIT NONE
!~ !INTEGER :: NtotF
!~ open(unit=2, file="total1.xyz", status="replace")
!~ NtotF= n+sum(ocvF)
!~    Do I=1,NtotF
!~    write(2,FMT=*)NtotF
!~    end Do
!~ print *, NtotF,":final occupation"
!~ close(2)
!~ end subroutine
END SUBROUTINE
  
SUBROUTINE ReadOccupationVector_initial()
implicit none
   open(unit = 2,file = "initial_NiH.xyz",status="replace")
   write(2,*)NtotI
   write(2,*)
   DO i=1,n
    write(2,*)"Ni",rNi(i,1),rNi(i,2),rNi(i,3)
     end Do
     
     
     TEMP_OCV = OCVF-OCVI
DO I=1,N
   if (TEMP_OCV(I)==-1)THEN
   A1=I
   ELSE IF(TEMP_OCV(I)==1)then
   B1=I
   END IF
END DO
WRITE(2,*)"H" ,RH(A1,:)

   TEMP_OCVI=OCVF+OCVI
      DO I=1,N
      IF (TEMP_OCVI(I)==2)THEN
      C1=I
      WRITE(2,*)"H",RH(C1,:)
      END IF
      END DO
      close(2)
      !WRITE(2,*)"H",RH(C1,:)
     
!~       write(2,*)"Ni",rNi(i,1),rNi(i,2),rNi(i,3)
!~       end Do
!~       Do I=1,N
!~       IF(OCVI(i) == 1) then
!~          write(2,*) "H",rH(i,1),rH(i,2),rH(i,3)
!~       END IF
!~       write(2,*) "Ni",rNi(i,1),rNi(i,2),rNi(i,3)
!~    END DO
!~    close(2)
!~    PRINT *,"WORKING"
   END SUBROUTINE
  
  !Step 3. Store the site occupation in array final_occupation_vector
  !CALL ReadOccupationVector(final_occupation_vector,final_occupation_file) !final 
  SUBROUTINE ReadOccupationVector_final()
  implicit none
  !INTEGER:: NtotF=13501
   open(unit = 2,file = "final_NiH.xyz",status="replace")
   write(2,*)NtotF
   write(2,*)
   DO i=1,n
      write(2,*)"Ni",rNi(i,1),rNi(i,2),rNi(i,3)
      end Do
      TEMP_OCV = OCVF-OCVI
   DO I=1,N
   if (TEMP_OCV(I)==-1)THEN
   A1=I
   ELSE IF(TEMP_OCV(I)==1)then
   B1=I
   END IF
   END DO
   WRITE(2,*)"H" ,RH(B1,:)
      TEMP_OCVI=OCVF+OCVI
      DO I=1,N
      IF (TEMP_OCVI(I)==2)THEN
      D1=I
      WRITE(2,*)"H",RH(D1,:)
      END IF
      END DO
      !WRITE(2,*)"H",RH(D1,:)
!**********************************************************
!H POSITION IN THE MIDDLE
   !DO i=1,n
!~       write(2,*)"Ni",rNi(i,1),rNi(i,2),rNi(i,3)
!~       end Do
!~       Do I=1,N
!~       IF(OCVF(i) == 1) then
!~          write(2,*)"H",rH(i,1),rH(i,2),rH(i,3)
!~       END IF
!~       write(2,*) "Ni",rNi(i,1),rNi(i,2),rNi(i,3)
!~    END DO
!**********************************************************
   close(2)
   END SUBROUTINE
  
  !Step 4. Create LAMMPS file - initial configuration
  !CALL CreateLAMMPSfileForRelaxation(initial_occupation_vector,initialLAMMPSfileForRelaxation)
  SUBROUTINE CreateLAMMPSfileForRelaxation_initial()
  IMPLICIT NONE
open(unit=2, file="initial_NiH.xyz", status="OLD")
open(unit=3, file="dump_initial2.lammps", status="REPLACE")
read(2,*)NtotI
read(2,*)
ALLOCATE(Species(NtotI),coord(NtotI,3))
   Do I=1,NtotI
   read(2,fmt=*)species(I),coord(I,1:3)
   end Do
 write(3,*)
 WRITE(3,*)NtotI, "atoms"
 write(3,*)"3 atom types"
 write(3,*)
 write(3,*)XLO, XHI,"xlo  ","xhi"
 write(3,*)YLO, YHI,"ylo  ","yhi"
 write(3,*)ZLO, ZHI,"zlo  ","zhi"
 write(3,*)
 write(3,*)"Masses"
 write(3,*)
 write(3,*)1,58.6934  
 write(3,*)2,26.9815
 write(3,*)3,1.00800
 write(3,*)
 write(3,*)"Atoms"
 WRITE(3,*)
   Do I=1,NtotI
   write(3,fmt=*)I,species(I),coord(I,:)
   end Do
CLOSE(2)
CLOSE(3)
call system(" sed -i 's\Ni\1\g' dump_initial2.lammps ")
call system(" sed -i 's\H\3\g' dump_initial2.lammps ")
END SUBROUTINE
  
  !Step 5. Create LAMMPS file - final configuration 
  !CALL CreateLAMMPSfileForRelaxation(final_occupation_vector,finalLAMMPSfileForRelaxation)
  SUBROUTINE CreateLAMMPSfileForRelaxation_final()
  IMPLICIT NONE
open(unit=2, file="final_NiH.xyz", status="OLD")
open(unit=3, file="dump_final2.lammps", status="REPLACE")
read(2,*)NtotF
read(2,*)
ALLOCATE(Speciess(NtotF),coordd(NtotF,3))
   Do I=1,NtotF
   read(2,fmt=*)speciess(I),coordd(I,1:3)
   end Do
 write(3,*)
 WRITE(3,*)NtotF, "atoms"
 write(3,*)"3 atom types"
 write(3,*)
 write(3,*)XLO, XHI,"xlo  ","xhi"
 write(3,*)YLO, YHI,"ylo  ","yhi"
 write(3,*)ZLO, ZHI,"zlo  ","zhi"
 write(3,*)
 write(3,*)"Masses"
 write(3,*)
 write(3,*)1,58.6934  
 write(3,*)2,26.9815
 write(3,*)3,1.00800
 write(3,*)
 write(3,*)"Atoms"
 WRITE(3,*)
   Do I=1,NtotF
   write(3,fmt=*)I,speciess(I),coordd(I,:)
   end Do
CLOSE(2)
CLOSE(3)
call system(" sed -i 's\Ni\1\g' dump_final2.lammps ")
call system(" sed -i 's\H\3\g' dump_final2.lammps ")
END SUBROUTINE
  
  !Step 6. Relax the initial configuration
  !CALL RelaxUsingLAMMPS(initialLAMMPSfileForRelaxation)('lmp -in in_NIH_initial.eam')
  
  !Step 7. Relax the final configuration
  !CALL RelaxUsingLAMMPS(finalLAMMPSfileForRelaxation)('lmp -in in_NIH_initial1.eam')
  
  !Step 8. Prepare LAMMPS files for NEB calculation
  !CALL PrepareForNEB()
  
  !subroutine
  
SUBROUTINE PrepareForNEB()
IMPLICIT NONE
 INTEGER:: Nhop,I
 INTEGER,DIMENSION(:),ALLOCATABLE:: X
 REAL,DIMENSION(:,:),ALLOCATABLE:: NH
 Nhop=1
OPEN(UNIT=2, FILE="neb_initial.lammps", STATUS="OLD")
!OPEN(UNIT=2, FILE="dump_initial2.lammps", STATUS="OLD")
CALL SYSTEM("sed '1,17d' neb_initial.lammps > ned_initial.lammps")
!CALL SYSTEM("sed '1,16d' dump_initial2.lammps > ned_initial.lammps")
CLOSE(2)
OPEN(UNIT=3, FILE="neb_final.lammps", STATUS="OLD")
!OPEN(UNIT=3, FILE="dump_final2.lammps", STATUS="OLD")
CALL SYSTEM("sed '1,17d' neb_final.lammps > ned_final.lammps")
!CALL SYSTEM("sed '1,16d' dump_final2.lammps > ned_final.lammps")
CLOSE(3)
OPEN(UNIT=2, FILE="ned_initial.lammps", STATUS="OLD")
CALL SYSTEM("awk '{print $1,$2}' ned_initial.lammps > ned11_initial.lammps")
CLOSE(2)
OPEN(UNIT=3, FILE="ned_final.lammps", STATUS="OLD")
CALL SYSTEM("awk '{print $1,$2,$3,$4,$5}' ned_final.lammps > ned22_final.lammps")
CLOSE(3)
OPEN(UNIT=2, FILE="ned11_initial.lammps", STATUS="OLD")
CALL SYSTEM("awk '{if($2>2) print $0}' ned11_initial.lammps > ned111_initial.lammps")
CLOSE(2)
OPEN(UNIT=3, FILE="ned22_final.lammps", STATUS="OLD")
CALL SYSTEM("awk '{if($2>2) print $0}' ned22_final.lammps > ned222_final.lammps")
CLOSE(3)
OPEN(UNIT=2, FILE="ned111_initial.lammps", STATUS="OLD")
CALL SYSTEM("awk '{print $1}' ned111_initial.lammps > ned1111_initial.lammps")
CLOSE(2)

OPEN(UNIT=3, FILE="ned222_final.lammps", STATUS="OLD")
CALL SYSTEM("awk '{print $1,$3,$4,$5}' ned222_final.lammps > ned22224_final.lammps")
!CLOSE(2)
CLOSE(3)

open(unit=2,file="ned22224_final.lammps",status="old")
call system("sort ned22224_final.lammps > ned22225_final.lammps")
close(2)
open(unit=3,file="ned22225_final.lammps",status="old")
call system("head -n 1 ned22225_final.lammps > ned66.lammps")
close(3)
OPEN (UNIT=3, FILE="ned66.lammps", status="OLD")
OPEN (UNIT=4, FILE="F1.lammps", status="replace")
!~ READ(3,*) N,X,Y,Z
!~  WRITE(4,*) 1 !NUMBER OF HYDROGEN 
!~  WRITE(4,*) N,X,Y,Z

!******************************************************
!! END 1H ATOM MOVING NEB
!*****************************************************
 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx
!HOPING OF MORE THAN 1 HATOMS
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!OPEN (UNIT=3, FILE="file3.txt", status="OLD")
!OPEN (UNIT=4, FILE="F2.lammps", status="replace")
!print *,"Nhop",Nhop
ALLOCATE(NH(Nhop,3),X(Nhop))
   DO I=1,Nhop
      READ(3,*)X(I),NH(I,1),NH(I,2),NH(I,3)
   END DO
   CLOSE(3)
 WRITE(4,*) Nhop !NUMBER OF HYDROGEN HOPING/MIGRATION
   DO I=1,Nhop
      WRITE(4,*)X(I),NH(I,1),NH(I,2),NH(I,3)
   END DO
CLOSE(4)
deallocate(NH)
deallocate(X)
END SUBROUTINE
  
!Step 9. Perform NEB
!CALL NEB('mpirun -np 9 --oversubscribe lmp -partition 9x1 -in input_NEB>t5.txt')
  
!Step 10. Extract the barrier from the LAMMPS output file
!CALL ExtractNEBBarrier("grep -e 'EBF' -e '0.559' t5.txt | tail -1 | awk '{print $7}'")
  
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!RELAXATION CALCULATION
!~ !*******************************************************
!~ CALL SYSTEM('ls')
!~ CALL SYSTEM('pwd')
!~ call system ('lmp -in in_NIH_initial.eam')
!~ call system ('lmp -in in_NIH_initial1.eam')
!!~ call system (sed -i   
!***********************************************************
! Final NEB calcilation
!***********************************************************
!call system ('mpirun -np 9 --oversubscribe lmp -partition 9x1 -in input_NEB')
!***********************************************************
!call system ('export LAMMPS_PYTHON_TOOLS=/home/sourabh/Downloads/lammps-23Jun2022/tools/python/pizza/')
!call system ('python3 /home/sourabh/Downloads/lammps-23Jun2022/tools/python/neb_final.py -o dump.hop1.final -r dumpX_neb_*')
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
END PROGRAM

