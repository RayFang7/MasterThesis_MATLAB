function output = threeDprint(sinput)
iinput = interp2(sinput,(1:.01:size(sinput,2))',1:0.3:size(sinput,1),'cubic');
iinput = imresize(iinput,[128 128]);
cmax = max(max(iinput,[],'all'),-min(iinput,[],'all'));
%  output = iinput - min(iinput,[],'all');
%  output = (output / (max(output,[],'all'))) * 256;
output = iinput +cmax; 
output = (output / (cmax*2)) * 255;
end