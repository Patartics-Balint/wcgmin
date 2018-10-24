function [M, blk_str, ne, nd, nz, nw] = components(sys, Kt)
    [M, Delta, blk_str] = lftdata(sys);
    M = prescale(M);
    blk_str = cat(1, blk_str.Size);
    [ny_sys, nu_sys] = size(sys);
    [nu, ny] = size(Kt);
    ne = ny_sys - ny;
    nd = nu_sys - nu;
    [nw, nz] = size(Delta);
    blk_str = [blk_str; nd, ne];
end