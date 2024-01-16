function Kz = KzInterp(Kz1, Kz2,ti, tf, t)
if t==ti
    Kz = Kz1;
elseif t==tf
    Kz = Kz2;
else
    Kz = Kz2.*(t/(tf-ti)) + Kz1.*(1 - (t/(tf-ti)));
end

Kz(isnan(Kz)) = 1e-2;

end