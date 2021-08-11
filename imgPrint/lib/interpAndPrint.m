%20201102 Cosplay contour finished
function output = interpAndPrint(sinput,cm)
iinput = interp2(sinput,(1:.01:size(sinput,2))',1:0.3:size(sinput,1),'cubic');
iinput = imresize(iinput,[512 512]);
iinput = normalize(iinput);
cmax = max(max(iinput,[],'all'),-min(iinput,[],'all'));
%https://www.mathworks.com/matlabcentral/fileexchange/20922-stlwrite-write-ascii-or-binary-stl-files
output = iinput +cmax;
output = (output / (cmax*2)) * 255;
output = ind2rgb(round(output'),cm);
end

%%
