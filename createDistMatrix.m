function DistMatrix =  createDistMatrix(previous,current)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N1=size(previous,1);
N2=size(current,1);
DistMatrix=zeros(N1,N2);
for i=1:N1
    for j=1:N2
        DistMatrix(i,j)=sqrt(sum((previous(i,:)-current(j,:)).^2));
    end
end

