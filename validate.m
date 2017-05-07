N=50
T=200

[x1minc,x1maxc,cgridx,c3unil,c3unilr,Ytrue, realD,mx,Y_Dtrue]=simu1(N, 1/12, 0.02, 0.03, 0.005, j*N*T+1,T, 256 ) ;            
Xc=[]
Yc=[]
for i=1:N
    Xc{i}= [ cgridx(:,i) ];
    Yc{i}=  [ c3unil(:,i) ];
end


%%Compare Eigenfunctions
%Get eigenvalues
X0mean=c3unilr -repmat( mean(c3unilr,2) ,1,N);
Xdmean=realD -repmat( mean(realD,2) ,1,N);

[VrealM0, DrealM0]=eig(  range(x)/length(x)*X0mean'*X0mean  );
[DrealM0,I]=sort(diag(DrealM0),'descend');
DrealM0=diag(DrealM0);
DDrealM0=DrealM0/N;
VrealM0=VrealM0(:,I);
pc0Real=Xdmean * VrealM0(:,1:2) * diag( diag(DrealM0(1:2,1:2)).^(-1/2) );


[hX2b, V2b, loadsb,Meansmo2bb,Db,x]=fpca(2,2,'M_0',Xc,Yc,mx',0,0,min(x1minc),max(x1maxc),'Gauss','Epan');


scatter(cgridx(:,38), c3unil(:,38))

diag( DDrealM0(1:2,1:2) )
Db(1:2)

plot( pc0Real(:,1))
hold on
plot( -V2b(:,1) )


mean( (hX2b-realD).^2 )
(diag( DrealM0(1:2,1:2) )-Db(1:2)).^2 
mse_ord(pc0Real(:,1),V2b(:,1),33) 
mse_ord(pc0Real(:,2),V2b(:,2),33) 

scatter(cgridx(:,i), c3unil(:,i))