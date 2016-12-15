function Cinf=bmd_reduced_thresh(x,sigz,sigx,filter2_size)
T=size(x,1);
xf=zeros(T,2);
xs=zeros(T,2);

%Kalman filter 1
Ez =sigz^2*eye(2);
Ex =sigx^2*eye(2);
P=(-Ez+sqrt(Ez.^2+4*Ez*Ex))/2;
K=(P+Ez)*inv(P+Ez+Ex);

xf(1,:)=x(1,:);
for t = 2:length(x)
    xf(t,:)=(xf(t-1,:)'+K*(x(t,:)-xf(t-1,:))')';
end

%backward part, Kalman smoother 1
K=P*inv(P+Ez);
xs(T,:)=xf(T,:);
for t = (T-1):-1:1
    xs(t,:)=(xf(t,:)'+K*(xs(t+1,:)-xf(t,:))')';
end

%Kalman filter 2
Ez=eye(2);
Ex=filter2_size*eye(2); 

xf=zeros(T,2);
P=(-Ez+sqrt(Ez.^2+4*Ez*Ex))/2;
K=(P+Ez)*inv(P+Ez+Ex);

xf(1,:)=xs(1,:);
for t = 2:length(x)
    xf(t,:)=xf(t-1,:)+(xs(t,:)-xf(t-1,:))*K;  
end
%backward, Kalman smoother 2
xs=zeros(T,2);
Ks=P*inv(P+Ez);

xs(T,:)=xf(T,:);
for t = (T-1):-1:1
    xs(t,:)=xf(t,:)+(xs(t+1,:)-xf(t,:))*Ks;
end


Cinf=sqrt([0; sum(diff(xs).^2,2)])>sigz+0.001;
% after matching and changing units from ms to s, we get the formula in the
% manuscript
end