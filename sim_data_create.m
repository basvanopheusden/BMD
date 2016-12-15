function [C,z,x]=sim_data_create(T,d0,d1, sigm0, sigm1, sigz, sigx)

Sigmaz=sigz^2*eye(2);

lambda0=0.004;
lambda1=0.1; 
puu=0:1:4000;
ns=500;
c_vec_found=false;

pdf_dur0=((exp(-lambda0)-1)^2)/exp(-lambda0)*(puu.*exp(-lambda0*puu));
pdf_dur1=((exp(-lambda1)-1)^2)/exp(-lambda1)*(puu.*exp(-lambda1*puu));

while ~c_vec_found
    dur1=randsample(puu,ns,true,pdf_dur1);
    dur0=randsample(puu,ns,true,pdf_dur0);
    c_vec_found=ismember(T,cumsum([0 dur0(1:end-1)]+dur1));    
end

C=[];
for i=1:ns
    C=[C ones(1,dur1(i)) zeros(1,dur0(i))];
end
C=C(1:T);

z=zeros(T,2);
a=zeros(T,2);
x=zeros(T,2);

meas_noise = sigx * randn(1,2);

r_jump=sigm1*sqrt(gamrnd((d1+1)/2, 2, [1 1]));
thet_j=2*pi*rand;


v(1,:)=[r_jump*cos(thet_j) r_jump*sin(thet_j)];
z(1,:)=[0.0001 0.0001];
x(1,:)=z(1,:);

for i=2:T
    if C(i)==C(i-1)   %no changepoint (either within drift seq or saccade seq)
        v(i,:)=v(i-1,:);    %v stays the same
    elseif C(i)~=C(i-1) && C(i)==1  %changepoint: start saccade sequence
        r_jump=sigm1*sqrt(gamrnd((d1+1)/2, 2, [1 1]));
        thet_j=2*pi*rand;
        v(i,:)=[r_jump*cos(thet_j) r_jump*sin(thet_j)];
    elseif   C(i)~=C(i-1) && C(i)==0 %changepoint: start drift sequence
        r_drift=sigm0*sqrt(gamrnd((d0+1)/2, 2, [1 1]));
        thet_d=2*pi*rand;
        v(i,:)=[r_drift*cos(thet_d) r_drift*sin(thet_d)];
    end
    z(i,:)=mvnrnd(z(i-1,:)+v(i,:), Sigmaz, 1); 
 
    meas_noise = sigx * randn(1,2);
    x(i,:) = z(i,:)+ meas_noise;
       
end


end