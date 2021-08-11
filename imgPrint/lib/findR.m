%20201104 R finding MARK II
function Rloc = findR(input)
Rloc = zeros(1,size(input,1));
ysum = zeros(1,size(input,1));
for i = 1:size(input,2)
     ecgsig = input(:,i);
    [~,qrs_i_raw,~]=pan_tompkin(ecgsig,1000,0);
    for j = -5:5
        p = min(max(1,qrs_i_raw+j),length(ecgsig));
        ysum(p)=ysum(p) + max(mean(ecgsig(qrs_i_raw))-mean(ecgsig),0.05);
    end
end
    MinPeakHeight = mean(maxk(ysum,round(size(ecgsig,1)/2000))) * 0.5;
    [~,Rloc] = findpeaks(ysum,1:length(ecgsig),'MinPeakHeight',MinPeakHeight,...
        'MinPeakDistance',200);
    %%
%     plot(1:length(ecgsig),y)
%     hold on
%     plot(Rloc,ecgsig(Rloc),'ro')
%     plot(ecgsig)
%     xlabel('Seconds')
%     title('R Peaks Localized by Wavelet Transform with Automatic Annotations')
    %%
end