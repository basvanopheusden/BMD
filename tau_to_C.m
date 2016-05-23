function Ct =tau_to_C(t01,t10, T)
Ct=zeros(1,T);
for j=1:length(t01)
Ct(t01(j):t10(j)-1)=ones(1,t10(j)-t01(j));
end
%perhaps could vectorize further

end
