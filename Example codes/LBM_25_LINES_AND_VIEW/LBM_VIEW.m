%% Sébastien Leclaire (2014)
clear all;close all;addpath('distFig');rng('default');NX=10; NY=10; OMEGA=1.0; rho0=1.0; deltaUX=10^-6; 
W=[4/9,1/9,1/36,1/9,1/36,1/9,1/36,1/9,1/36];
cx=[0,0,1,1, 1, 0,-1,-1,-1];cy=[0,1,1,0,-1,-1,-1, 0, 1];
SOLID=rand(NX*NY,1)>0.7;
figure(1);figure(2);distFig('Rows',1);

N=bsxfun(@times,rho0*ones(NX*NY,9),W)+0*1/36*(rand(NX*NY,9)-0.5);
[x,y]=ind2sub([NX NY],1:NX*NY);x=x';y=y';
x=repmat(x,1,9);
y=repmat(y,1,9);
cmap=jet(9)*0+0.5;

figure(1);clf
hold on;
for i=1:NX*NY
    if SOLID(i)
        p = patch([x(i)-0.5 x(i)+0.5 x(i)+0.5 x(i)-0.5 x(i)-0.5],...
            [y(i)-0.5 y(i)-0.5 y(i)+0.5 y(i)+0.5 y(i)-0.5], [0 0 0]);
    end
end
for i=2:9
    CNORM=hypot(cx(i),cy(i));
    Nx=bsxfun(@times,N(:,i),cx(i)*CNORM);%./sum(N,2);
    Ny=bsxfun(@times,N(:,i),cy(i)*CNORM);%./sum(N,2);
    hfull{i}=quiver(x(:,i),y(:,i),Nx*4,Ny*4,'color',cmap(i,:),'autoscale','off','LineWidth',2);
end
axis image
axis([0 NX+1 0 NY+1])
axis off
plot([0.5 NX+.5 NX+.5 0.5 0.5],[0.5 0.5 NY+.5 NY+.5 0.5],'k-','LineWidth',2)
title('Initial condition.','fontsize',20)
drawnow
pause(2)

for t_=1:4000    
    pauseTime=1;
    %% Collision
    
    figure(1);clf;hold on
    for i=1:NX*NY
        if SOLID(i)
            p = patch([x(i)-0.5 x(i)+0.5 x(i)+0.5 x(i)-0.5 x(i)-0.5],...
                [y(i)-0.5 y(i)-0.5 y(i)+0.5 y(i)+0.5 y(i)-0.5], [0 0 0]);
        end
    end
    for i=2:9
        CNORM=hypot(cx(i),cy(i));
        Nx=bsxfun(@times,N(:,i),cx(i)*CNORM);
        Ny=bsxfun(@times,N(:,i),cy(i)*CNORM);
        hfull{i}=quiver(x(:,i),y(:,i),Nx*4,Ny*4,'color',cmap(i,:),'autoscale','off','LineWidth',2);
    end
    axis image
    axis([0 NX+1 0 NY+1])
    axis off
    plot([0.5 NX+.5 NX+.5 0.5 0.5],[0.5 0.5 NY+.5 NY+.5 0.5],'k-','LineWidth',2)
    drawnow
    
    N_SOLID=N(SOLID,:); % No Collision
    
    rho = sum(N,2);
    ux  = sum(bsxfun(@times,N,cx*CNORM),2)./rho;ux=ux+deltaUX;
    uy  = sum(bsxfun(@times,N,cy*CNORM),2)./rho;
    workMatrix=ux*cx+uy*cy;workMatrix=(3+4.5*workMatrix).*workMatrix;
    workMatrix=bsxfun(@minus,workMatrix,1.5*(ux.^2+uy.^2));
    workMatrix=bsxfun(@times,1+workMatrix,W);
    workMatrix=bsxfun(@times,workMatrix,rho);
    N=N+(workMatrix-N)*OMEGA;
    
    N(SOLID,:)=N_SOLID;   % Dummy line only for correct visualisation of Bounce Back.
     %% Verbosity
      
    title('Collision in 3 seconds.','fontsize',20)
    fluidPlot=plot(x(~SOLID),y(~SOLID),'ro','MarkerSize',30);
    drawnow    
    pause(3*pauseTime)
    title('Collision!','fontsize',20)
    drawnow
      
    for i=2:9
        delete(hfull{i});
        CNORM=hypot(cx(i),cy(i));
        Nx=bsxfun(@times,N(:,i),cx(i)*CNORM);
        Ny=bsxfun(@times,N(:,i),cy(i)*CNORM);
        hfull{i}=quiver(x(:,i),y(:,i),Nx*4,Ny*4,'color',cmap(i,:),'autoscale','off','LineWidth',2);
    end    
    drawnow
    pause(1*pauseTime);
    delete(fluidPlot)
    title('Bounce Back in 3 seconds.','fontsize',20)
    solidPlot=plot(x(SOLID),y(SOLID),'yo','MarkerSize',30);
    drawnow
    pause(3*pauseTime);
    title('Bounce Back!','fontsize',20)         
    drawnow
        
    %% Bounce Back
    
    N(SOLID,:)=N_SOLID(:,[1 6 7 8 9 2 3 4 5]);
    
    totalStep=15;
    for nstep=1:totalStep
        for i=2:9
            delete(hfull{i})
            CNORM=hypot(cx(i),cy(i));
            Nx=bsxfun(@times,N(:,i),cx(i)*CNORM);
            Ny=bsxfun(@times,N(:,i),cy(i)*CNORM);
        
            theta=pi*nstep/totalStep;
            rotationMatrix=[cos(theta) -sin(theta);sin(theta) cos(theta)];
            NxRot=rotationMatrix(1,1)*Nx+rotationMatrix(1,2)*Ny;
            NyRot=rotationMatrix(2,1)*Nx+rotationMatrix(2,2)*Ny;
        
            NxRot(~SOLID)=Nx(~SOLID);
            NyRot(~SOLID)=Ny(~SOLID);
            
            hfull{i}=quiver(x(:,i),y(:,i),NxRot*4,NyRot*4,'color',cmap(i,:),'autoscale','off','LineWidth',2);            
        end
        axis image
        axis([0 NX+1 0 NY+1])
        axis off
        drawnow
    end
    
    pause(1*pauseTime);
    delete(solidPlot);
    title('Streaming in 3 seconds.','fontsize',20)
    drawnow
    pause(3*pauseTime);
    title('Streaming!','fontsize',20)         
    drawnow
    
    %% Streaming
    
    for i=2:9
        delete(hfull{i})
        CNORM=hypot(cx(i),cy(i));
        Nx=bsxfun(@times,N(:,i),cx(i)*CNORM);
        Ny=bsxfun(@times,N(:,i),cy(i)*CNORM);
        xMove=x;
        yMove=y;
        totalStep=15;
        for nstep=1:totalStep
            dt=1/totalStep;
                                    
            xMove(:,i)=xMove(:,i)+dt*cx(i);
            yMove(:,i)=yMove(:,i)+dt*cy(i);
           
            xMove(:,i)=(xMove(:,i)>NX+.5).*(xMove(:,i)-NX)+(1-(xMove(:,i)>NX+.5)).*xMove(:,i);
            xMove(:,i)=(xMove(:,i)<0.5).*(xMove(:,i)+NX)+(1-(xMove(:,i)<0.5)).*xMove(:,i);
            yMove(:,i)=(yMove(:,i)>NY+.5).*(yMove(:,i)-NY)+(1-(yMove(:,i)>NY+.5)).*yMove(:,i);
            yMove(:,i)=(yMove(:,i)<0.5).*(yMove(:,i)+NY)+(1-(yMove(:,i)<0.5)).*yMove(:,i);

            h=quiver(xMove(:,i),yMove(:,i),Nx*4,Ny*4,'color',cmap(i,:),'autoscale','off','LineWidth',2);
            axis image
            axis([0 NX+1 0 NY+1])  
            axis off
            drawnow            
            if nstep~=totalStep
                delete(h)
            end
        end        
    end
    
    for i=2:9
        N(:,i)=reshape(circshift(reshape(N(:,i),NX,NY),[cx(i),cy(i)]),NX*NY,1);
    end
    
    ux(SOLID)=0; uy(SOLID)=0;ux=reshape(ux,NX,NY)';uy=reshape(uy,NX,NY)';
    figure(2);clf;hold on;colormap(gray(2));image(2-reshape(SOLID,NX,NY)');
    quiver(1:NX,1:NY,ux,uy,'b','LineWidth',2);axis([0.5 NX+0.5 0.5 NY+0.5]);axis image;
    title(['Velocity field after ',num2str(t_),' time steps'],'fontsize',20);
    drawnow    
end
