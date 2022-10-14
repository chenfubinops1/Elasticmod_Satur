function ksat=gassm(k0,kdry,kf,por)
ksat=kdry+(1-kdry./k0).^2./(por./kf+(1-por)./k0-kdry./k0^2);

end