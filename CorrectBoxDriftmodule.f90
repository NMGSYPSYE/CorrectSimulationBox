module  CorrctBoxDriftModule
!***  This is a module program.
  
implicit none
   
integer, dimension(:), pointer::ityp   !the type of the atoms in box
real*8, dimension(:,:), pointer::Xp    !the position of the atoms in box
real*8, dimension(:,:), pointer::Xpb   !the position of sites in reference lattice system

!***the axial ratio(Hcpca=c/a) for HCP
real*8::Hcpca


!*** positon deviation in x, y, z directions
real*8::Boxdrift(3)  
     
!***center position of all atoms in the unit cell 
real*8::CentrPos(3)  

!***positon of site in reference lattice system
real*8::SiteRefLatt(3) 

!*** the number and position of lattice sites in a single cell,and position of atoms in a single cell
!    sitesnum=9 for BCC , sitesnum=14 for FCC, sitesnum=17 for HCP
integer::sitesnum
real*8, dimension(:,:), pointer::LatticeSite
real*8, dimension(:,:), pointer::Unitcellatoms
!***position of atoms in a single BCC cell
 
!***N: number of all atoms in the box; M: number of lattice sites in the perfect lattice system.
integer::N,M

real*8, dimension(:,:), pointer::Ncellatoms  
real*8, dimension(:,:), pointer::Ncellatomstmp 
real*8, dimension(:),   pointer::Ncelldis 
integer::Ncellnum=0
	
contains 


!***To read in the data from the input files.
subroutine Initialize(Input_AtomsPosi_name,Input_ReferLatt_name)
implicit none
   character*256, intent(in)::Input_AtomsPosi_name,Input_ReferLatt_name

   real*8::x,y,z
   integer::I,J

   open(unit=1, file=Input_AtomsPosi_name, status='old') 
   N = 0
   DO WHILE(.not.eof(1))
        read(1,*, end=100) j,X, Y, Z 
        N = N+1
     END DO
100 rewind(1)	 
    allocate(ityp(N),  xp(N,3))
   
    DO I=1, N
       read(1,*, end=100)  ityp(I),xp(I,1:3)
    END DO 		 
    close(1)
 
    open(unit=2, file=Input_ReferLatt_name, status='old') 
	 M = 0
     DO WHILE(.not.eof(2))
        read(2,*, end=200) X, Y, Z 
        M = M+1
     END DO
 200 rewind(2)	 
     allocate(Xpb(M,3))
   

	 DO I=1, M
        read(2,*, end=100)  Xpb(I,1:3)
     END DO 		 
    close(2)

    return
end subroutine Initialize
    

!***To determine the relative positions of each sites in a perfect unit cell according to crystal structure type.
subroutine CrystalType_Num_Sites(CrystalType)
implicit none
   integer, intent(in)::CrystalType
   select case (CrystalType)
	        case (1)                
                sitesnum=9
                allocate(LatticeSite(sitesnum,3))
                allocate(Unitcellatoms(sitesnum,3))
                LatticeSite(1,1:3)=0.0D0

                LatticeSite(2,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(2,2)=LatticeSite(1,2)-0.5D0    
                LatticeSite(2,3)=LatticeSite(1,3)-0.5D0

                LatticeSite(3,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(3,2)=LatticeSite(1,2)+0.5D0 
                LatticeSite(3,3)=LatticeSite(1,3)-0.5D0     

                LatticeSite(4,1)=LatticeSite(1,1)+0.5D0      
                LatticeSite(4,2)=LatticeSite(1,2)+0.5D0      
                LatticeSite(4,3)=LatticeSite(1,3)-0.5D0
     
                LatticeSite(5,1)=LatticeSite(1,1)+0.5D0     
                LatticeSite(5,2)=LatticeSite(1,2)-0.5D0    
                LatticeSite(5,3)=LatticeSite(1,3)-0.5D0
       
                LatticeSite(6,1)=LatticeSite(1,1)-0.5D0     
                LatticeSite(6,2)=LatticeSite(1,2)-0.5D0
                LatticeSite(6,3)=LatticeSite(1,3)+0.5D0
      
                LatticeSite(7,1)=LatticeSite(1,1)-0.5D0       
                LatticeSite(7,2)=LatticeSite(1,2)+0.5D0  
                LatticeSite(7,3)=LatticeSite(1,3)+0.5D0
       
                LatticeSite(8,1)=LatticeSite(1,1)+0.5D0    
                LatticeSite(8,2)=LatticeSite(1,2)+0.5D0
                LatticeSite(8,3)=LatticeSite(1,3)+0.5D0
       
                LatticeSite(9,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(9,2)=LatticeSite(1,2)-0.5D0
                LatticeSite(9,3)=LatticeSite(1,3)+0.5D0
                
            case (2)
                sitesnum=14
                allocate(LatticeSite(sitesnum,3))
                allocate(Unitcellatoms(sitesnum,3)) 
                LatticeSite(1,1:3)=0.0D0

                LatticeSite(2,1)=LatticeSite(1,1)+1.0D0
                LatticeSite(2,2)=LatticeSite(1,2)+0.0D0
                LatticeSite(2,3)=LatticeSite(1,3)+0.0D0
      
                LatticeSite(3,1)=LatticeSite(1,1)+1.0D0
                LatticeSite(3,2)=LatticeSite(1,2)+0.0D0
                LatticeSite(3,3)=LatticeSite(1,3)+1.0D0    

                LatticeSite(4,1)=LatticeSite(1,1)+0.0D0   
                LatticeSite(4,2)=LatticeSite(1,2)+0.0D0
                LatticeSite(4,3)=LatticeSite(1,3)+1.0D0      

                LatticeSite(5,1)=LatticeSite(1,1)+0.0D0
                LatticeSite(5,2)=LatticeSite(1,2)+1.0D0
                LatticeSite(5,3)=LatticeSite(1,3)+0.0D0   

                LatticeSite(6,1)=LatticeSite(1,1)+1.0D0
                LatticeSite(6,2)=LatticeSite(1,2)+1.0D0
                LatticeSite(6,3)=LatticeSite(1,3)+0.0D0   

                LatticeSite(7,1)=LatticeSite(1,1)+1.0D0
                LatticeSite(7,2)=LatticeSite(1,2)+1.0D0
                LatticeSite(7,3)=LatticeSite(1,3)+1.0D0 

                LatticeSite(8,1)=LatticeSite(1,1)+0.0D0      
                LatticeSite(8,2)=LatticeSite(1,2)+1.0D0       
                LatticeSite(8,3)=LatticeSite(1,3)+1.0D0
      
                LatticeSite(9,1)=LatticeSite(1,1)+0.5D0      
                LatticeSite(9,2)=LatticeSite(1,2)+0.5D0     
                LatticeSite(9,3)=LatticeSite(1,3)+0.0D0
   
                LatticeSite(10,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(10,2)=LatticeSite(1,2)+0.5D0
                LatticeSite(10,3)=LatticeSite(1,3)+1.0D0
      
                LatticeSite(11,1)=LatticeSite(1,1)+0.0D0
                LatticeSite(11,2)=LatticeSite(1,2)+0.5D0
                LatticeSite(11,3)=LatticeSite(1,3)+0.5D0
      
                LatticeSite(12,1)=LatticeSite(1,1)+1.0D0
                LatticeSite(12,2)=LatticeSite(1,2)+0.5D0
                LatticeSite(12,3)=LatticeSite(1,3)+0.5D0

                LatticeSite(13,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(13,2)=LatticeSite(1,2)+0.0D0
                LatticeSite(13,3)=LatticeSite(1,3)+0.5D0

                LatticeSite(14,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(14,2)=LatticeSite(1,2)+1.0D0
                LatticeSite(14,3)=LatticeSite(1,3)+0.5D0

            case (3)
                sitesnum=17
                allocate(LatticeSite(sitesnum,3))
                allocate(Unitcellatoms(sitesnum,3)) 
                LatticeSite(1,1:3)=0.0D0

                LatticeSite(2,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(2,2)=LatticeSite(1,2)-3.0D0**0.5D0/2.0D0    
                LatticeSite(2,3)=LatticeSite(1,3)+0.0D0

                LatticeSite(3,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(3,2)=LatticeSite(1,2)-3.0D0**0.5D0/2.0D0    
                LatticeSite(3,3)=LatticeSite(1,3)+0.0D0

                LatticeSite(4,1)=LatticeSite(1,1)+1.0D0
                LatticeSite(4,2)=LatticeSite(1,2)+0.0D0    
                LatticeSite(4,3)=LatticeSite(1,3)+0.0D0

                LatticeSite(5,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(5,2)=LatticeSite(1,2)+3.0D0**0.5D0/2.0D0  
                LatticeSite(5,3)=LatticeSite(1,3)+0.0D0

                LatticeSite(6,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(6,2)=LatticeSite(1,2)+3.0D0**0.5D0/2.0D0  
                LatticeSite(6,3)=LatticeSite(1,3)+0.0D0

                LatticeSite(7,1)=LatticeSite(1,1)-1.0D0
                LatticeSite(7,2)=LatticeSite(1,2)+0.0D0  
                LatticeSite(7,3)=LatticeSite(1,3)+0.0D0

                LatticeSite(8,1)=LatticeSite(1,1)+0.0D0
                LatticeSite(8,2)=LatticeSite(1,2)+0.0D0  
                LatticeSite(8,3)=LatticeSite(1,3)+Hcpca

                LatticeSite(9,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(9,2)=LatticeSite(1,2)-3.0D0**0.5D0/2.0D0    
                LatticeSite(9,3)=LatticeSite(1,3)+Hcpca

                LatticeSite(10,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(10,2)=LatticeSite(1,2)-3.0D0**0.5D0/2.0D0    
                LatticeSite(10,3)=LatticeSite(1,3)+Hcpca

                LatticeSite(11,1)=LatticeSite(1,1)+1.0D0
                LatticeSite(11,2)=LatticeSite(1,2)+0.0D0    
                LatticeSite(11,3)=LatticeSite(1,3)+Hcpca

                LatticeSite(12,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(12,2)=LatticeSite(1,2)+3.0D0**0.5D0/2.0D0  
                LatticeSite(12,3)=LatticeSite(1,3)+Hcpca

                LatticeSite(13,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(13,2)=LatticeSite(1,2)+3.0D0**0.5D0/2.0D0  
                LatticeSite(13,3)=LatticeSite(1,3)+Hcpca

                LatticeSite(14,1)=LatticeSite(1,1)-1.0D0
                LatticeSite(14,2)=LatticeSite(1,2)+0.0D0  
                LatticeSite(14,3)=LatticeSite(1,3)+Hcpca

                LatticeSite(15,1)=LatticeSite(1,1)-0.5D0
                LatticeSite(15,2)=LatticeSite(1,2)+3.0D0**0.5D0/6.0D0  
                LatticeSite(15,3)=LatticeSite(1,3)+Hcpca/2.0D0

                LatticeSite(16,1)=LatticeSite(1,1)+0.0D0
                LatticeSite(16,2)=LatticeSite(1,2)-3.0D0**0.5D0/3.0D0  
                LatticeSite(16,3)=LatticeSite(1,3)+Hcpca/2.0D0

                LatticeSite(17,1)=LatticeSite(1,1)+0.5D0
                LatticeSite(17,2)=LatticeSite(1,2)+3.0D0**0.5D0/6.0D0  
                LatticeSite(17,3)=LatticeSite(1,3)+Hcpca/2.0D0

   end select  
   return
               
end subroutine CrystalType_Num_Sites


!***To Search for atoms in the simulation box, which are able to maintain the crystal structure.
subroutine Search_For_Atoms_CrystalType(CrystalType,Radius,subatomtype) 
implicit none
   real*8, intent(in)::Radius
   integer, intent(in)::subatomtype,CrystalType
   integer::I,II,J,K,L,LL,LLL,L4,number,number1,type1,type2,MinNcellnum
   real*8::disx,disy,disz,dis,dis1,DisSitSitUCel,DisAtomAtmUniCeL
   real*8, dimension(:,:), pointer::LatticeSiteTmp
   
   call CrystalType_Num_Sites(CrystalType)
   allocate(LatticeSiteTmp(sitesnum,3))
   Ncellnum=0
   DO I=1,N
      type1=ityp(I)
      if(type1.ne.subatomtype) cycle
      DO II=1, sitesnum
         LatticeSiteTmp(II,1:3)=xp(I,1:3)+LatticeSite(II,1:3)
      END DO

      number=0

      DO J=1,sitesnum
          
         number1=0
         DO K=1,N
            type2=ityp(K)
            if(type2.ne.subatomtype) cycle

               disx=LatticeSiteTmp(J,1)-xp(K,1)
               disy=LatticeSiteTmp(J,2)-xp(K,2)
               disz=LatticeSiteTmp(J,3)-xp(K,3)
               dis=DSQRT(disx*disx+disy*disy+disz*disz)
               
               if (dis.lt.radius) then
                   Unitcellatoms(J,1:3)=xp(K,1:3)
                   number1=number1+1
               end if                             
          END DO

          if (number1.eq.1) then                     
              number=number+1
          else
              exit
          end if

       END DO

        
       if(number.eq.sitesnum) then
           Ncellnum=Ncellnum+1
           call Search_For_Many_UnCell(Ncellnum,sitesnum,Unitcellatoms)    
       end if


   END DO

   call Determine_Centrpos(Unitcellatoms,sitesnum)
 
   deallocate (LatticeSiteTmp)
   nullify(LatticeSiteTmp)
   return
end subroutine Search_For_Atoms_CrystalType


!***To search the multiple cells that meet the condition (based on the parameter 'Radius') in the simulation box.
subroutine Search_For_Many_UnCell(Ncellnum,sitesnum,Unitcellatoms)
implicit none
   !real*8, intent(in)::Radius
   real*8,intent(in)::Unitcellatoms(:,:)
   integer, intent(in)::Ncellnum,sitesnum
   integer::L,LL,L4
   if(associated(Ncellatoms))  deallocate (Ncellatoms)
   nullify(Ncellatoms)
   allocate(Ncellatoms(sitesnum*Ncellnum,3))
   IF (Ncellnum.EQ.1) THEN
       Ncellatoms=Unitcellatoms
       if(associated(Ncellatomstmp))  deallocate (Ncellatomstmp)
       nullify(Ncellatomstmp)
       allocate(Ncellatomstmp(sitesnum,3))
       Ncellatomstmp=Ncellatoms
    ELSE
       DO L=1,sitesnum*(Ncellnum-1)
          Ncellatoms(L,1:3)= Ncellatomstmp(L,1:3)                
       END DO

       DO LL=sitesnum*(Ncellnum-1)+1,sitesnum*Ncellnum
          L4=LL-sitesnum*(Ncellnum-1)
          Ncellatoms(LL,1:3)=Unitcellatoms(L4,1:3)
          if(associated(Ncellatomstmp))  deallocate (Ncellatomstmp)
          nullify(Ncellatomstmp)
          allocate(Ncellatomstmp(sitesnum*Ncellnum,3))
          Ncellatomstmp=Ncellatoms

       END DO
   END IF
end subroutine Search_For_Many_UnCell


!***To search a single cell that best maintains the crystal structure in the simulaiton box (based on the  
!   distance between atoms in a perfect cell)and determine its center position.     
subroutine Determine_Centrpos(Unitcellatoms,sitesnum)
implicit none
   real*8::DisSitSitUCel,disx,disy,disz,dis,dis1
   real*8,intent(in)::Unitcellatoms(:,:)
   integer, intent(in)::sitesnum
   integer::I,J,K,L,LL,L4,MinNcellnum
   CentrPos=0.0D0
   DisSitSitUCel=0.0D0
   if (Ncellnum.GE.1) then
       allocate(Ncelldis(Ncellnum))
       Ncelldis=0.0D0
       DO I=1,sitesnum
           DO J=I+1,sitesnum
              disx=LatticeSite(I,1)-LatticeSite(J,1)
              disy=LatticeSite(I,2)-LatticeSite(J,2)
              disz=LatticeSite(I,3)-LatticeSite(J,3)
              dis=DSQRT(disx*disx+disy*disy+disz*disz)
              DisSitSitUCel=dis+DisSitSitUCel
           END DO
       END DO
 
       DO K=1,Ncellnum
          DO I=(Ncellnum-1)*sitesnum+1,sitesnum*Ncellnum
             DO J=I+1,sitesnum*Ncellnum 
                disx=Ncellatoms(I,1)-Ncellatoms(J,1)
                disy=Ncellatoms(I,2)-Ncellatoms(J,2)
                disz=Ncellatoms(I,3)-Ncellatoms(J,3)
                dis=DSQRT(disx*disx+disy*disy+disz*disz)
                Ncelldis(K)=dis+Ncelldis(K)
             END DO
          END DO
       END DO
       
       dis=Ncelldis(1)-DisSitSitUCel
       dis=abs(dis)
       MinNcellnum=1
       DO K=2,Ncellnum
           dis1=Ncelldis(K)-DisSitSitUCel
           dis1=abs(dis1)
           if(dis1.LE.dis) then
              dis=dis1
              MinNcellnum=K
           end if
       END DO

       DO I=(MinNcellnum-1)*sitesnum+1,sitesnum*MinNcellnum
          CentrPos(1)=Ncellatoms(I,1)+CentrPos(1)
          CentrPos(2)=Ncellatoms(I,2)+CentrPos(2)
          CentrPos(3)=Ncellatoms(I,3)+CentrPos(3)
       END DO
          
       CentrPos(1)=CentrPos(1)/sitesnum
       CentrPos(2)=CentrPos(2)/sitesnum
       CentrPos(3)=CentrPos(3)/sitesnum
   else
       write(*,*) "Error, please adjust the value of parameter Radius"
       stop
   end if
end subroutine Determine_Centrpos


!*** To determine the position deviation between the reference lattice system and the simulation box.
subroutine determine_the_drift_of_box()
implicit none
   integer::I
   real*8::disx,disy,disz,dis,dis1
   Boxdrift=0.0D0
   SiteRefLatt=0.0D0
   disx=CentrPos(1)-xpb(1,1)
   disy=CentrPos(2)-xpb(1,2)
   disz=CentrPos(3)-xpb(1,3)
   dis=SQRT(disx*disx+disy*disy+disz*disz)

   SiteRefLatt(1)=xpb(1,1)
   SiteRefLatt(2)=xpb(1,2)
   SiteRefLatt(3)=xpb(1,3)

   DO I=2,M
      disx=CentrPos(1)-xpb(I,1)
      disy=CentrPos(2)-xpb(I,2)
      disz=CentrPos(3)-xpb(I,3)
      dis1=SQRT(disx*disx+disy*disy+disz*disz)
          
      if (dis1.le.dis) then
          SiteRefLatt(1)=xpb(I,1)
          SiteRefLatt(2)=xpb(I,2)
          SiteRefLatt(3)=xpb(I,3)
          dis=dis1
      end if

   END DO

   Boxdrift(1:3)=CentrPos(1:3)-SiteRefLatt(1:3)
   return
end subroutine determine_the_drift_of_box


!***To correct the postions of all atoms in the simulation box according to the boundary type.
subroutine Correct_ALL_Atoms_Position(LA,LB,LC,BoundType)
implicit none
   real*8, intent(in)::LA,LB,LC
   integer, intent(in)::BoundType
   select case (BoundType)
	       case (1)
                call Correct_ALL_Atoms_Position_BounTy1(LA,LB,LC)
           case (2)
                call Correct_ALL_Atoms_Position_BounTy2(LA,LB,LC)
   end select

   return
end subroutine  Correct_ALL_Atoms_Position



!***To correct the postions of all atoms in the simulation box according the deviation for the boundary type of 1. 
subroutine Correct_ALL_Atoms_Position_BounTy1(LA,LB,LC)
implicit none
   real*8, intent(in)::LA,LB,LC
   real*8::HalfLA,HalfLB,HalfLC
   HalfLA=LA/2
   HalfLB=LB/2
   HalfLC=LC/2    
    
   xp(1:N,1)=xp(1:N,1)-Boxdrift(1)
   xp(1:N,2)=xp(1:N,2)-Boxdrift(2)
   xp(1:N,3)=xp(1:N,3)-Boxdrift(3)

   where(dabs(xp(1:N,1)) .GT. HalfLA)                                                                      
        xp(1:N,1) = xp(1:N,1) - DSIGN(LA,XP(1:N,1))                 
   end where 
    
   where(dabs(XP(1:N,2)) .GT. HalfLB)                                                                      
        XP(1:N,2) = XP(1:N,2) - DSIGN(LB,XP(1:N,2))                 
   end where 
   
   where(dabs(XP(1:N,3)) .GT. HalfLB)                                                                      
         XP(1:N,3) = XP(1:N,3) - DSIGN(LC,XP(1:N,3))                 
   end where 
   
   return
end subroutine Correct_ALL_Atoms_Position_BounTy1


!***To correct the postions of all atoms in the simulation box according the deviation for the boundary type of 2.  
subroutine Correct_ALL_Atoms_Position_BounTy2(LA,LB,LC)
   implicit none
   real*8, intent(in)::LA,LB,LC
   real*8::HalfLA,HalfLB,HalfLC
   HalfLA=LA/2
   HalfLB=LB/2
   HalfLC=LC/2    
    
   xp(1:N,1)=xp(1:N,1)-Boxdrift(1)
   xp(1:N,2)=xp(1:N,2)-Boxdrift(2)
   xp(1:N,3)=xp(1:N,3)-Boxdrift(3)

   where(xp(1:N,1).LT.0.0D0)                                                                      
        xp(1:N,1) = xp(1:N,1) +LA                 
   end where 
    
   where(XP(1:N,2).LT. 0.0D0)                                                                      
        XP(1:N,2) = XP(1:N,2) +LB                  
   end where 
   
   where(XP(1:N,3).LT. 0.0D0)                                                                      
         XP(1:N,3) = XP(1:N,3)+LC                 
   end where 

   where(xp(1:N,1).GT.LA)                                                                      
        xp(1:N,1) = xp(1:N,1) -LA                 
   end where 
    
   where(XP(1:N,2).GT. LB)                                                                      
        XP(1:N,2) = XP(1:N,2) -LB                  
   end where 
   
   where(XP(1:N,3).GT. LC)                                                                      
         XP(1:N,3) = XP(1:N,3)-LC                 
   end where 

   return
end subroutine Correct_ALL_Atoms_Position_BounTy2


!***To write the data into the output file.
subroutine  Write_InTo_Output_File(Output_AtomsPosi_name)
implicit none
   character*256, intent(in)::Output_AtomsPosi_name
   integer::I

   open(unit=3, file=Output_AtomsPosi_name, status='unknown') 
   DO I=1, N
       write(3,FMT='(I5,3(1x,1PE14.6))') ityp(I),xp(I,1:3) 
   END DO 		 
   close(3)
   if(associated(Ncellatoms))     deallocate (Ncellatoms)
   if(associated(Ncellatomstmp))  deallocate (Ncellatomstmp)
   if(associated(ityp))           deallocate (ityp)
   if(associated(xp))             deallocate (xp)
   if(associated(xpb))            deallocate (xpb)
   if(associated(Ncelldis))       deallocate (Ncelldis)
   if(associated(LatticeSite))    deallocate (LatticeSite)
   if(associated(Unitcellatoms))  deallocate (Unitcellatoms)
    
   nullify(ityp,xp,xpb,Ncellatoms,Ncellatomstmp,Ncelldis,LatticeSite,Unitcellatoms)
   return
end subroutine Write_InTo_Output_File


end module  CorrctBoxDriftModule
