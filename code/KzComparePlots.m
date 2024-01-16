for ii = 1:11
    figure
    subplot(1,3,1)
    plot(CTD.sig0(:,ii),CTD.de(:,ii)); axis ij;
    hold on
    subplot(1,3,2)
    plot(CTD.bvfrq(:,ii),CTD.de(:,ii)); axis ij;
    hold on; plot(XX(:,ii),CTD.de(:,ii), "-.");
    subplot(1,3,3)
    plot(Kz_bfrq(:,ii),CTD.de(:,ii)); axis ij;
    hold on
    plot(Kzxx(:,ii),CTD.de(:,ii))
end

figure
plot(log10(Kzxx(:,1)), CTD.de(:,1));
hold on
plot(log10(Kzfilt(:,1)), CTD.de(:,1), "--")
axis ij

