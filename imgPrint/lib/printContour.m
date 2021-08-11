function rgb = printContour(M,cm)
        [~,c]=contour(1:length(M),1:6,M',256);
        cmax = max(max(M,[],'all'),-min(M,[],'all'));
        set(c,'LineColor','none')
        set(c,'Fill','on')
        colormap(cm);
        axis tight
        axis off
        caxis([-cmax,cmax]);
        F = getframe;
        rgb = F.cdata;
end