function MinDistPos = checkMatr(DistMatr,R)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
[MinDist,MinDistPos]=min(DistMatr,[],1);
if length(MinDistPos)~=length(unique(MinDistPos))
    for i=1:length(MinDistPos)
        ElemList=i;
        for j=i+1:length(MinDistPos)
            if MinDistPos(j)==MinDistPos(i)
                ElemList=[ElemList j];
            end
        end
        if length(ElemList)>1
            Dist=MinDist(ElemList);
            [~,distpos]=min(Dist);
            for j=1:length(ElemList)
                if j~=distpos
                    MinDistPos(ElemList(j))=NaN;
                end
            end
            
        end
    end
end
MinDistPos(MinDist>R)=NaN;

