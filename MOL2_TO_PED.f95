        PROGRAM MOL2PED
        character*100 cabmol,nom_mol,cabmol2,fitxer
        character*11 cod
        character*4 orb(3000)
        character*2 str2(3000),at(3000),tipus(3000),str(3000)
        integer n,i,id(3000),inputstat,contador,inputstat2,llig1(1000),inllig(1000),sbond
        real x(3000),y(3000),z(3000),charge(3000)
        real xc,yc,zc,d(3000),u11,u12,u13,u14,u15,u16,u21,u22,u23,u24,u25,u26
        real u32,u33,u43,u44,u34,u51,u52,u53,u54,pc(20)
                
!       Obrim els fitxers de llista classificada i mol2 per tal de fer la lectura
	contador=0
        fitxer='query.mol2'
        inputstat=0
        inputstat2=0
        
	open(unit=4,file=fitxer)
        open(unit=6,file='query.shd')
        !open(unit=7,file='DB_bin.shd',form='unformatted',action='write',access='transparent')
        call cpu_time(t1) 
        
23      continue
        IF (InputStat < 0) go to 30
	IF (InputStat > 0) then
           write(*,*)"error!"
           go to 30
         end if
	 
         read(4,*,iostat=inputstat)cabmol
         IF (InputStat < 0) go to 30
	 IF (InputStat > 0) then
            write(*,*)"error!"
            go to 30
          end if
	  
          if(cabmol.eq."@<TRIPOS>MOLECULE")then
            read(4,'(a23)')cabmol(1:23)
            nom_mol=cabmol(1:23)      
            read(4,*)N,Nl
29	    continue
            read(4,*)cabmol2
            	if(cabmol2.eq."@<TRIPOS>ATOM")then
                  do i=1,N
       		    READ(4,*,iostat=inputstat2)id(I),at(i),x(i),y(i),z(i),orb(i),str(i),str2(i),charge(i)  
                    read(4,*)cabmol2
                    do ii=1,Nl*2,2
                      read(4,*)inllig(ii),llig1(ii),llig1(ii+1),tipus(ii)                     
                    enddo  		
             	  call stat8(id,at,x,y,z,orb,str,str2,charge,n,nl,inllig,llig1,tipus,pc,nom_mol)
                  write(6,'(12(f9.4,2x),a23)')pc(1),pc(2),pc(3),pc(4),pc(5),&
                  &pc(6),pc(7),pc(8),pc(9),pc(10),pc(11),pc(12),nom_mol(1:23)
                  go to 23 
	        else
                  go to 29
                end if
          else 
            go to 23
          end if         
          
30	continue    

	close(4)
        close(6)
        close(7)
        close(88)

        call cpu_time(t2)
        open(83,file='temps_ped.txt')        
        write(83,'(f9.4)')t2-t1
        close(83)
            
      	end

      subroutine stat8(id,at,x,y,z,orb,str,str2,charge,n,nl,inllig,llig1,tipus,pc,nom_mol)        
      integer id(3000),i,k
      integer llig1(1000),inllig(1000),nd,nc2,ii,nl,n
      character*2 str2(3000),at(3000),str(3000),tipus(3000)
      character*4 orb(3000),corb
      character*100 nom_mol
      real x(3000),y(3000),z(3000),charge(3000),u(500000),pc(20),dist    
      
      k=0
      do i=1,n-1
       do j=i+1,n
         dist=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)  
         if(dist.gt.1.0) then
           if((abs(charge(i)).gt.0).or.(abs(charge(j)).gt.0)) then
       	    k=k+1
            u(k)=charge(i)*charge(j)/sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)*14.4      
           end if
         end if
       enddo 
      enddo      
        
      do i=1,k-1
        do j=1,k-i
          if(u(j).gt.u(j+1)) then
            aux = u(j)
            u(j) = u(j+1)
            u(j+1)=aux
           end if
         enddo
       enddo    

        pc(1)=u(1)!+u(int(k/100)*1-1))/2d0
        ! percentil 5%
        pc(2)=u(2)!(u(floor(k/100.*5.)+1)+u(floor(k/100.*5.)+2))/2d0!+u(int(k/100)*5-1))/2d0
        ! percentil 10%
        pc(3)=u(3)!(u(floor(k/100.*10.)+1)+u(floor(k/100.*10.)+2))/2d0!u(int(k/100*10)+1)!+u(int(k/100)*10-1))/2d0        
		! percentil 20%
        pc(4)=u(4)!(u(floor(k/100.*20.)+1)+u(floor(k/100.*20.)+2))/2d0!+u(int(k/100)*20-1))/2d0
        ! percentil 30%
        pc(5)=u(5)!(u(floor(k/100.*30.)+1)+u(floor(k/100.*30.)+2))/2d0
        ! percentil 50%
        pc(6)=u(6)!(u(floor(k/100.*50.)+1)+u(floor(k/100.*50.)+2))/2d0        
        ! percentil 70%
        pc(7)=u(k-5)!(u(floor(k/100.*70.)+1)+u(floor(k/100.*70.)+2))/2d0      
        ! percentil 80%
        pc(8)=u(k-4)!(u(floor(k/100.*80.)+1)+u(floor(k/100.*80.)+2))/2d0         
        ! percentil 90%
        pc(9)=u(k-3)!(u(floor(k/100.*90.)+1)+u(floor(k/100.*90.)+2))/2d0       
        ! percentil 95%
        pc(10)=u(k-2)!(u(floor(k/100.*95.)+1)+u(floor(k/100.*98.)+1))/2d0         
        ! percentil 99%
        pc(11)=u(k-1)!(u(floor(k/100.*98.)+1)+u(k))/2d0  
        ! max
        pc(12)=u(k)
               
        end
