%%%%%%%%%%%%%%%%%%%%%%%%%
%        2 steps Lax-Wendroff scheme for              %
%                    1D AQ system                           %
%                        4/3/2019                              %
%%%%%%%%%%%%%%%%%%%%%%%%%

global dx alpha nu rho sigma E R0 h G0 A0 L c dt T1 Xin Xout

figure(1); clf;
figure(2); clf;

% Read input data obtained using fft
Xin = dlmread('ProximalPressure_fft.dat');
Xout = dlmread('DistalPressure_fft.dat');
externalPressure =(83.34+(188.92-83.34)/3)*1333;%(Max-Min)/3 in dynes/cm^2 
disp(externalPressure)

N=26;%# of elements
alpha = 4/3;
nu = 0.035; %viscosity (cm^2/s)
rho = 1.055; %density (g/cm^3)
sigma = 0.49; %Poisson ratio
E = 2e6; %4e6; %Young's modulus (dynes/cm^2)
R0 = 0.18; %Radius (cm)
h = 0.0475; %Thickness (cm)
G0 = (E*h)/R0/(1-sigma^2);
A0 = R0^2;
L = 4.5; %Artery length (cm)
dx= L/N;
CFLfactor=0.5;
c = sqrt(G0/rho); %sound speed
dt = dx/c*CFLfactor; %CFL
% Here we used only c for the CFL condition, instead of the full 
% eigenvalues lambda_1 or lambda_2 because the other terms in the 
% eigenvalue calculation are much smaller that c and do not
% influence the CFL condition much

Ncycles=3; %# of cycles we are calculating; 
% it takes at least one cycle
% for the solution to "stabilize" (transient period)

T1=.763; %Time interval for one cycle
T = T1*Ncycles;
x = 0:dx:L;
t = 0:dt:T;
index=1; %movie frame index
index2=1;

%%
A = zeros(1,length(x))+R0^2; %A at time level n
Q = zeros(1,length(x)); %Q at time level n
A1= zeros(1,length(x)); %A at time level n+1
Q1= zeros(1,length(x)); %Q at time level n+1
U = zeros(2,length(x));
U1 = zeros(2,length(x));
U(1,:) = A;
U(2,:) = Q;

for n = 1:length(t)-1
    for i = 1:length(x)
        if i == 1 %Left boundary (x=0)
            
            %Write the code for boundary data at 0

            U1(1,i) = A1(i);
            U1(2,i) = Q1(i);
            
        elseif x(i) ==L  %Right boundary (x=L)
            
            %Write the code for boundary data at L
              
            U1(1,i) = A1(i);
            U1(2,i) = Q1(i);
        
        else %Inner elements (using the two step lax-wendroff scheme)
            
            U_midA = U(:,i); %MODIFY
            U_midB = U(:,i); %MODIFY

            U1(:,i) = 0.0; %MODIFY

            A1(i) = U1(1,i);
            Q1(i) = U1(2,i);
            
        end
    end
    A = A1; %update A, Q and U
    Q = Q1;
    U(1,:) = A;
    U(2,:) = Q;

  %plotting
   figure(1);
   width = 800; height = 600;
   set(figure(1),'Position',[15 15 width height]);
   if n==1 %For 1st layout
      subplot(2,2,1); plot(x,A);hold on;
      subplot(2,2,1); plot(x,ones(length(x))*R0*R0);hold off;
      axis([0 L 0.02 0.05]); title('Cross-sectional Area');
      xlabel('x (cm)'); ylabel('A (cm^2)');
      subplot(2,2,2); plot(x,Q); 
      axis([0 L -4 16]); title('Flow');
      xlabel('x (cm)'); ylabel('Q (cm^3/s)');
      subplot(2,2,3); plot(t,Inlet(t)/1333,'b');hold on; 
      subplot(2,2,3); plot(t,Outlet(t)/1333,'r');hold on;
      axis([0 T 50 200]); title('Proximal/Distal Pressure');
      ylabel('mmHg'); xlabel('t (sec)');
      subplot(2,2,4);plot(t(n),Q(1,N/2)/A(1,N/2)/3.14,'b-+');hold on;
      axis([0 T -20 160]); title('Velocity at mid-point');
      xlabel('t (sec)'); ylabel('V (cm/s)');
   end
   
   Nskip=200; %number of time steps to skip for plotting;
   if mod(n,Nskip)==1
     subplot(2,2,1); plot(x,A);hold on;
     subplot(2,2,1); plot(x,ones(length(x))*R0*R0);hold off;
     axis([0 L 0.02 0.05]); title('Cross-sectional Area');
      xlabel('x (cm)'); ylabel('A (cm^2)');
     subplot(2,2,2); plot(x,Q); 
     axis([0 L -4 16]); title('Flow');
     xlabel('x (cm)'); ylabel('Q (cm^3/s)');
     subplot(2,2,3);h2=plot([t(n) t(n)],[0, 250],'Color',[0 0 0]);
     subplot(2,2,4);plot(t(n),Q(1,N/2)/A(1,N/2)/3.14,'b-+');hold on;
     axis([0 T -20 160]); title('Velocity at mid-point');
      xlabel('t (sec)'); ylabel('V (cm/s)');
     drawnow;
     disp(t(n));
     
     if n > (length(t)-1)/Ncycles
     figure(2);
     width = 800; height = 600;
     set(figure(2),'Position',[15 15 width height]);
     subplot(2,1,1);
     plot(t(n),Q(1,N/2)/A(1,N/2)/3.14,'*');hold on;
     axis([0.8 1.58 -20 160]); title('Velocity at mid-point');
     xlabel('t (sec)'); ylabel('V (cm/s)');
     subplot(2,1,2);
     plot(t(n),A(1,N/2),'*');hold on;
     axis([0.8 1.58 0.02 0.05]); title('Cross-sectional Area at mid-point');
     xlabel('t (sec)'); ylabel('A (cm^2)');
     drawnow;
     MovieFrame2(index2) = getframe(figure(2)) ;
     index2=index2+1;
     end
     
     
%Next 3 lines are generating movie frames     
     MovieFrame1(index) = getframe(figure(1)) ;
     set(h2,'Visible','off')
     index=index+1;

   end
end

%Next 8 lines are saving the movie to an AVI file
writerObj = VideoWriter('Case1.avi');
writerObj.FrameRate = 15;
open(writerObj);
for i=1:length(MovieFrame1)
    frame = MovieFrame1(i) ;    
    writeVideo(writerObj, frame);
end
for i=1:length(MovieFrame2)
    frame = MovieFrame2(i) ;    
    writeVideo(writerObj, frame);
end
close(writerObj);




%End of main procedure


%%%%%%%%%%%%%%%
%         Subroutines            %
%%%%%%%%%%%%%%%

% calculating the flux term (in conservative form)
function f= flux(Utru) % F = [F1,F2]
global alpha rho R0 G0 A0 
Atru = Utru(1);
Qtru = Utru(2);
f = [Qtru; (alpha*Qtru^2/Atru)+(   G0/3/rho*((Atru/A0)^(3/2))*A0)];
end

%calculate the RHS
function s = rhs(Utru) 
global alpha nu 
Atru = Utru(1);
Qtru = Utru(2);
s = [0;(-2*nu*(alpha/(alpha-1))*(Qtru/Atru))];
end

% Return pressure boundary data at Inlet
function InP = Inlet(t)
global Xin T1
N = length(Xin);
if mod(N,2) == 0
    n = [0:(N/2 - 1) (-N/2):-1]';
else
    n = [0:floor(N/2) -floor(N/2):-1]';
end
InP = 2000+real(sum(Xin .* exp(2. * pi * 1i * n * mod(t,T1) / T1))) / N;
end

% Return pressure boundary data at Outlet
function OutP = Outlet(t)
global Xout T1
N = length(Xout);
if mod(N,2) == 0
    n = [0:(N/2 - 1) (-N/2):-1]';
else
    n = [0:floor(N/2) -floor(N/2):-1]';
end
OutP = real(sum(Xout .* exp(2. * pi * 1i * n * mod(t,T1) / T1))) / N;
end
