function Rmax = findMax(input,type)
%被发现是没啥卵用的算法。
maxPoint = zeros(0,size(input,2));
for i = 1:size(input,2)
    [~,Index] = max(input(:,i));
    maxPoint(i) = Index;
end
Rmax = round(sum(maxPoint)/size(input,2));
if type == 'A'
    %for lead 8 only
    extrMaxValue = input(find(diff(sign(diff(input)))==-2)+1);
    extrMaxIndex = find(diff(sign(diff(input)))==-2)+1;
    [minValue,~] = min(input);
    area = round(0.1*length(input));
    for i = 1:length(extrMaxIndex)
        [~,maxIndex] = max(extrMaxValue);
        [testMinValue,~] = min(input(maxIndex:min(extrMaxIndex(maxIndex)+area,length(input))));
        if testMinValue ~= minValue
            extrMaxValue(maxIndex) = -99999;
        else
            Rmax = maxIndex;
            break;
        end
    end
end
end