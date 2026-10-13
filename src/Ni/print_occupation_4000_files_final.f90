PROGRAM MAIN
   IMPLICIT NONE
   INTEGER :: occupation_in_1NNN1(18),occupation_vector(4000)
   INTEGER :: SiteIndex(18)=[969,970,971,974,975,976,1009,1011,1015,1369,1370,1372,1373,1375,1376,1409,1413,1415]
   INTEGER :: isite,InitialHoppingSiteIndex=1374,FinalHoppingSiteIndex=972,ipos,I
   !!!!!! Here intitial hopingsiteindex=1 and finalhoppingindex=3 this is importrant
   !!!! if initial hopping is from site 2 then Initialhopping index should be 2
   !!!!! and final hopping index should be based on choice of neighboursites then change accordingly
   !Do I=1,2048
   
   !Step 1. Read occupation_in_1NN.txt
   OPEN(UNIT=524,FILE='occupation_in_1NNN1.txt')
   READ(524,*) occupation_in_1NNN1
   CLOSE(524)
   
   !Step 2. Update occupation
   occupation_vector=0
   !!!
   !occupation_vector(1:4000:1)=1
   
     DO I=1,18
   occupation_vector(SiteIndex(I)) = 0
   END DO
   occupation_vector(InitialHoppingSiteIndex)=1
   occupation_vector(FinalHoppingSiteIndex)=0
   
   !occupation_vector(InitialHoppingSiteIndex)=1
   DO ipos=1,18
      IF (occupation_in_1NNN1(ipos)==1) THEN !update the occupation at this site
         isite=SiteIndex(ipos)
         occupation_vector(isite)=1
      END IF
   END DO
   
   !Step 3. Print initial file
   !OPEN(UNIT=524,FILE='file1')
   !OPEN(UNIT=524,FILE='array1_initial.xyz')
   OPEN(UNIT=524,FILE='initial_file')
   DO ipos=1,4000
      WRITE(524,*) occupation_vector(ipos)
   END DO
   CLOSE(524)
   
   !Step 4. Update occupation vector for file2
   occupation_vector(InitialHoppingSiteIndex)=0
   occupation_vector(FinalHoppingSiteIndex)=1
   
   
   !Step 5. Print final file
   !OPEN(UNIT=524,FILE='file2')
   !OPEN(UNIT=524,FILE='array1_final.xyz')
   OPEN(UNIT=524,FILE='final_file')
   DO ipos=1,4000
      WRITE(524,*) occupation_vector(ipos)
   END DO
   CLOSE(524)
   !END DO
END PROGRAM MAIN
