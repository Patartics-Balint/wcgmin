function [K, m] = stabmax_syn(M, Kt, Dr_ss, Dc_ss, opt)
    % STABMAX_SYN Solve the synthesis in the stability maximisation
    % problem
    %
    % [K, m] = STABMAX_SYN(M, Kt, Dr_ss, Dc_ss, opt)
    %   Solves the optimisation below with structured controller Kt.
    % 
    %                ||                                           ||
    %                ||     +----+   +---------+   +---------+    ||
    %                || <---+ Dr <---+         <---+ inv(Dc) <--+ ||
    %                ||     +----+   |         |   +---------+    ||
    %                ||              +    M    +                  ||
    %     1/m = min  ||              x         x                  ||
    %                ||              +         +                  ||
    %           Kt   ||      +-------+         <-------+          ||
    %                ||      |       +---------+       |          ||
    %                ||      |                         |          ||
    %                ||      |          +----+         |          ||
    %                ||      +----------> Kt +---------+          ||
    %                ||                 +----+                    ||
    %                ||                                           ||
    %
    % [K, m] = STABMAX_SYN(M, ny, nu, Dr_ss, Dc_ss, opt)
    %   NOT IMPLEMENTED YET
    %   Solves the optimisation above with the unctructured controller Kt.    
    clt = close_stab_max_loop(M, Kt, Dr_ss, Dc_ss);
    [cl, m_inv] = hinfstruct(clt, opt.hinfstructOptions);
    K = getValue(Kt, cl);
    m = 1 / m_inv;
end