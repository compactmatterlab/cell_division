function dv1dt = dynamics(t,v1)
global Fmagpr dist0 storedistc storeFcph storeperipheryforce storecentforce cph N_p gammax gammay gammaz rad resize Fmagk Fmagd Fmagcr no_o no_c no_d no_k no_cr spe spa pcell  Fmagat Fmagatt  volume1 c2cd Fmagce Tfb Fmagpa Fmagpe storeFc2p storeFrepc storeFrepcr storeFcr2c storeFrpa storeFrpe storeFatp storeFcont storeF_pe storeSA storeSAt

v = reshape(v1,resize); % reshape(first size, final size)
                                              
c(1:no_c,:)=v(1:no_c,:);                                
cr(1:no_cr,:)=v(no_c+1:no_c+no_cr,:);                   
zz(1:N_p,:)=v(no_c+no_cr+1:no_c+no_cr+N_p,:);

N_p=size(zz,1);

%Centroid
Mcent=[mean(zz(:,1)),mean(zz(:,2)),mean(zz(:,3))];  

%Convex Hull
[K2,volume2] = convhulln(zz); 
                               
%P cell ratio
P_cell=pcell*volume1/volume2; % multiply pcell constant by ratio of volume initial/volume final


% Initial Triangulation
TR = triangulation(Tfb,zz);  

%Unit normal vectors to the vertices
VN = vertexNormal(TR);             

%Surface Area and pressure per unit  
                                      % code from https://www.mathworks.com/matlabcentral/answers/184234-how-do-i-determine-the-surface-area-of-a-2-d-surface-in-a-3-d-space
V1 = zz(TR(:,2), :) - zz(TR(:,1), :); % takes coordinate vectors of each triangle gives x y z of each triangle
V2 = zz(TR(:,3), :) - zz(TR(:,2), :); % Edges of each triangle
cp = 0.5*cross(V1,V2);                % cross product of each traingle of the surface vectors. gets area of each triangle
area_TR = sqrt(dot(cp, cp, 2));       % gets magnitude of area of each triangle
perp_TR = cp./area_TR;                % unit vector= vector/magnitude
F_TR = (P_cell*area_TR).*perp_TR;     %(F_TR= pressure*magnitude*unit vector). surface area of each triangle is multiplied by unit vector to indicate the direction and area so you get the magnitude)


%% Force due to internal pressure (F=P*area)
 for i=1:N_p    
     
         F_pe(i,:) = (1/3)*sum(F_TR(TR(:,1) == i | TR(:,2) == i | TR(:,3) == i,:),1); %force on each triangle. divide by three for each triangle vertex

 end
 

 %% Contractile force between periphery points (represents viscoelastic properties)
 Fcon = zeros(N_p,N_p,3);
 for n=1:N_p
     for m=1:N_p
         
         Np2Np2(n,m,:)=(zz(n,:)-zz(m,:)); %position vectors between each periphery point
         if n~=m & sum((TR(:,1) == n | TR(:,2) == n | TR(:,3) == n) & (TR(:,1) == m | TR(:,2) == m | TR(:,3) == m),1)>0
             
             %Normalize (Np2Np)
             test2(:)=Np2Np2(n,m,:);
             dist1=norm(test2);
             test2=test2(:)/dist1;
             Np2Np2(n,m,:)=test2(:);

             Fcon(n,m,:)=-(Fmagatt)./(dist0(n,m)).*(dist1-dist0(n,m)).*Np2Np2(n,m,:);
         else
             Fcon(n,m,:)=[0,0,0];
         end
     end
     Fcont(n,:)=sum(Fcon(n,:,:),2);
 end
 
 
%% Attractive force from the Substrate Plate                 
zt=0;  %threshold distance             
for w = 1:N_p
    
    Fp2cent(w,:)= zz(w,:)-Mcent(1,:);
    
    %Normalize (Fp2cent)
    test9(:)=Fp2cent(w,:);
    test9=test9(:)/norm(test9);
    Fp2cent(w,:)=test9(:);
    
    if zz(w,3)<Mcent(1,3) & zz(w,3)>zt %Pulling towards the plate
        
        Fatp(w,:)=(Fmagat*[0,0,-1]); %2/7/20
          
    elseif zz(w,3)<=zt %pushing outward after touching the plate
        
        Fatp(w,:)=1*[50*exp(-2*t)*Fp2cent(w,1), 50*exp(-2*t)*Fp2cent(w,2), 350*(zt-zz(w,3))]; % 350 represents stiffness of the plate
        
    end
    
end

Fatp2=Fatp;
Fatp2(:,3)=0; 


%% Coverplate Force. Confinement in z direction
for w = 1:N_p
    
    Fp2cph(w,:)=zz(w,:)-Mcent(1,:); %new pk edit
    

    if zz(w,3)>=cph %Pulling towards the plate
        
        Fcph(w,:)=1*[0*Fp2cph(w,1), 0*Fp2cph(w,2), 350*(cph-zz(w,3))]; % 350 represents stiffness of the coverplate
    else
        Fcph(w,:)=Fp2cph(w,:)*0;
        
    end
    
end

Fcph=Fcph*0;            %Comment this out when adding cover plate 
Fcph2=Fcph;
Fcph2(:,3)=0;

%% total forces on each periphery point
forces=Fcont+F_pe+Fatp2+Fcph2; %Fcont=attractive contractility between periphery points, F_pe= pressure, Fatp=attraction of plate, Fcph2=force from cover plate

Fav=vecnorm(forces,2,2);

fmagpe=(Fav.^3)./sum(Fav.^3); %Magnitude of periphery. external forces/ magnitude of force. Shows how large the force is compared to all the forces


%% External forces from periphery (kinesin and dynein) onto centrosome 
for r=1:no_o
    
    for q=1: no_c 
        %o2c origin to centrosomes
        o2c(r,q,:)=(c(q,:)-Mcent(r,:));
        
        %normalize o2c
        test3(:)=o2c(r,q,:);
        test3=test3(:)/norm(test3);
        o2c(r,q,:)=test3(:);
        
        for u=1:N_p
            %position of vector centrosomes to periphery
            c2p(q,u,:)=(zz(u,:)-c(q,:));
            
            %normalize c2p 
            test4(:)=c2p(q,u,:);  
            dist12=norm(test4);
            test4=test4(:)/dist12;
            c2p(q,u,:)=test4(:);
            
            if 2*pi/3 >= acos(dot(c2p(q,u,:),o2c(r,q,:))) % 2*pi/3>= angle between the two vectors
                
                force_factorp(q,u) = 1;    %either attached or not attached to centrosome 1 or 2
            else
                force_factorp(q,u) = 0; 
            end
         
             Fc2p(q,u,:)=((Fmagd*fmagpe(i,:)*no_d)>(Fmagk*rad/dist12*no_k/N_p))*(Fmagd*fmagpe(i,:)*no_d)*c2p(q,u,:)*force_factorp(q,u) - ((Fmagd*fmagpe(i,:)*no_d)<=(Fmagk*rad/dist12*no_k/N_p))*(Fmagk*rad/dist12*no_k/N_p)*c2p(q,u,:)*force_factorp(q,u); %Dynein is dependent on external forces and kinesin is evenly spread out on cortex 
        end
    end
end


%% Repulsive force centrosome to centrosome
for q=1:no_c
    for w=1:no_c
        if w~=q
            %vector of centrosome to centrosome
            ce2ce(q,w,:)=(c(q,:)-c(w,:));
            
            % normalize (Np2Np)
            test5(:)=ce2ce(q,w,:);
            dist2=norm(test5);      
            test5=test5(:)/dist2;
            ce2ce(q,w,:)=test5(:);
                        
            if dist2<=c2cd               
                Frepc(q,w,:)=(Fmagce/dist2)*ce2ce(q,w,:); % repulsion increases when distance decreases 
            else
                Frepc(q,w,:)=0*ce2ce(q,w,:);
            end
        end
    end
end



%% Force of Chromosome to Centrosome cr2c
for q=1: no_c
    
    for u= 1:no_cr
        %cr2c chromosome to centrosomes
        cr2c(q,u,:)=(c(q,:)-cr(u,:));
        
        %normalize (cr2c)
        test6(:)=cr2c(q,u,:);
        dist3=norm(test6);
        test6=test6(:)/dist3;
        cr2c(q,u,:)=test6(:);
        
        Fcr2c(q,u,:)=Fmagcr*dist3*cr2c(q,u,:); %more motors with distance
    end
end


%% Repulsive force of chromosome to periphery. Boundary check.
for q=1:N_p
      
        distp(q,:)=norm(zz(q,:)-Mcent(1,:)); %magnitude of distance from periphery to centroid
                  
end

mindist=min(distp);                          %minimum magnitude from periphery to centroid

for w=1:no_cr
    cr2mcent(w,:)=(cr(w,:)-Mcent(1,:));      %chromosome to centroid distance
    
    test10(:)=cr2mcent(w,:);
    dist5=norm(test10);
    test10=test10(:)/dist5;
    cr2mcent(w,:)=test10(:);
    
    if dist5>=(0.8*mindist)
        Frepcr(w,:)=Fmagpr*cr2mcent(w,:);
    else
        Frepcr(w,:)=0*cr2mcent(w,:);
    end
end

%% Position Vector Chromosome to chromosome and repulsive forces 
% perpandicular and parallel to the spindle axis (centrosomes)
for n=1:no_c
    for m=1:no_c
        if m~=n
            c2c(n,m,:)= c(n,:)-c(m,:);
            %normalize (c2c)
            test7(:)=c2c(n,m,:);
            distc=norm(test7);            
            test7=test7(:)/distc;
            c2c(n,m,:)=test7(:);   
        end
    end
end

for q=1: no_cr
    for u=1:no_cr
        if u~=q
            %cr2cr chromosome to chromosome
            cr2cr(q,u,:)=(cr(q,:)-cr(u,:));

            %normalize (cr2cr)
            test8(:)=cr2cr(q,u,:);          
            dist4(q,u)=norm(test8);           
            test8=test8(:)/dist4(q,u);        
            cr2cr(q,u,:)=test8(:);           

        end
    end
end

for q=1: no_cr
    for u=1:no_cr
        if u~=q
            for n=1:no_c
                for m=1:no_c
                    if m~=n
                        dd=dot(cr2cr(q,u,:),c2c(n,m,:)); 
                        if dist4(q,u)<=spa                                 
                            Frpa(q,u,:)=(Fmagpa/(dist4(q,u)))*c2c(n,m,:);         %Chromosome Repulsion Parallel
                        else                                                     
                            Frpa(q,u,:)=0*cr2cr(q,u,:);
                        end
                        if dist4(q,u)<=spe                                 
                             nn(:)=cr2cr(q,u,:)-(dd*c2c(n,m,:));                      
                            Frpe(q,u,:) =Fmagpe/(dist4(q,u))*(nn(:))/norm(nn(:)); %Chromosome Repulsion Perpandicular
                        else
                            Frpe(q,u,:)=0*cr2cr(q,u,:);
                        end
                    end
                end
            end
        end
    end
end


%% periphery force storage
for s=1:N_p
    tt=floor(t)+1;                            %floor rounds nearest integer less than or equal to t
    storeperipheryforce(1,tt,s,:)=Fatp(s,:);  %Fatp (Force category, time, periphery points, xyz)
    storeperipheryforce(2,tt,s,:)=Fcont(s,:);
    storeperipheryforce(3,tt,s,:)=F_pe(s,:);
    storeperipheryforce(4,tt,s,:)=Fcph(s,:);
end

%% centrosome force storage
    Fcc2p=sum(Fc2p(:,:,:),2);
    Fcc2p=reshape(Fcc2p,[2,3]);
    
    Fccr2c=sum(Fcr2c(:,:,:),2);
    Fccr2c=reshape(Fccr2c,[2,3]);
    
    Fcrepc=sum(Frepc(:,:,:),2);
    Fcrepc=reshape(Fcrepc,[2,3]);

for ab=1:no_c    
    tt=floor(t)+1;                           %floor rounds nearest integer less than or equal to t
    storecentforce(1,tt,ab,:)=Fcc2p(ab,:);   %Fc2p (Force category, time, centrosome 1 to 2, xyz)
    storecentforce(2,tt,ab,:)=-Fccr2c(ab,:); %Fcr2c
    storecentforce(3,tt,ab,:)=Fcrepc(ab,:);  %Frepc
end

%%store all forces
storeFc2p=Fc2p;
storeFcr2c=Fcr2c;
storeFrepc=Frepc;
storeFrepcr=Frepcr;
storeFrpa=Frpa;
storeFrpe=Frpe;
storeFatp=Fatp;
storeFcont=Fcont;
storeF_pe=F_pe;
storedistc=distc;
storeFcph=Fcph;
storeSA=sum(area_TR);
storeSAt=area_TR;

%% ODE
for i=1:no_c  %centrosomes 
    dv1dt(i) = (sum(Fc2p(i,:,1))-sum(Fcr2c(i,:,1))+sum(Frepc(i,:,1)))/(gammax*2);             %Fc2p(i,:,1)=(no_c,periphery,xyz) 
    dv1dt(i+size(v,1)) = (sum(Fc2p(i,:,2))-sum(Fcr2c(i,:,2))+sum(Frepc(i,:,2)))/(gammay*2);  
    dv1dt(i+2*size(v,1)) = (sum(Fc2p(i,:,3))-sum(Fcr2c(i,:,3))+sum(Frepc(i,:,3)))/(gammaz*2); %2*gamma to slow down centrosomes  
    1;
end

for i=no_c+1:no_c+no_cr %chromosomes (Chromosomes are only attached to centrosomes)
    j=i-no_c;           %number of chromosomes
    dv1dt(i) = (-Frepcr(j,1)+sum(Fcr2c(:,j,1))+sum(Frpa(j,:,1))+sum(Frpe(j,:,1)))/(gammax);             %Fcr2c(:,j,1)=(centrosomes, chromosomes, xyz) 
    dv1dt(i+size(v,1)) =(-Frepcr(j,2)+sum(Fcr2c(:,j,2))+sum(Frpa(j,:,2))+sum(Frpe(j,:,2)))/(gammay);    
    dv1dt(i+2*size(v,1)) = (-Frepcr(j,3)+sum(Fcr2c(:,j,3))+sum(Frpa(j,:,3))+sum(Frpe(j,:,3)))/(gammaz); 
    2;
end

for i=no_c+1+no_cr:no_c+no_cr+N_p %periphery                     
    j=i-no_c-no_cr;               %number of periphery points
    dv1dt(i) =(Fatp(j,1)+Fcph(j,1)+Fcont(j,1)+F_pe(j,1))/gammax;   
    dv1dt(i+size(v,1)) =(Fatp(j,2)+Fcph(j,2)+Fcont(j,2)+F_pe(j,2))/gammay; 
    dv1dt(i+2*size(v,1)) =(Fatp(j,3)+Fcph(j,3)+Fcont(j,3)+F_pe(j,3))/gammaz;
    3; 
end

dv1dt=dv1dt';

end
