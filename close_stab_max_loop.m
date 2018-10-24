function cl = close_stab_max_loop(M, K, Dr, Dc, m, g_sqrt)
    % CLOSE_STAB_MAX_LOOP Forms the closed loop interconnection used in the
    % stability maximasation.
    %
    % Syntax
    %   cl = CLOSE_STAB_MAX_LOOP(M, K, Dr, Dc)
    %   cl = CLOSE_STAB_MAX_LOOP(M, K, Dr, Dc, m, g_sqrt)
    %
    % Description
    %   cl = CLOSE_STAB_MAX_LOOP(M, K, Dr, Dc) forms the closed loop
    %   interconnection below. The performance channels of the plant are
    %   deleted.
    %
    %        +----+   +-----------+   +---------+
    %    <---+ Dr <---+           <---+ inv(Dc) <---
    %        +----+   |           |   +---------+
    %                 |     M     |     
    %                 x           x
    %                 |           |
    %        +--------+           <--------+
    %        |        +-----------+        |
    %        |                             |
    %        |            +---+            |
    %        +------------> K +------------+
    %                     +---+
    %   cl = CLOSE_STAB_MAX_LOOP(M, K, Dr, Dc, m, g_sqrt) forms the
    %   closed loop interconnection seen below.
    %
    %        +----------+   +-----------+   +---------------+
    %    <---+sqrt(m)*Dr<---+           <---+inv(Dc)*sqrt(m)<---
    %        +----------+   |           |   +---------------+
    %                       |           |
    %        +----------+   |     M     |     +----------+
    %    <---+ 1/g_sqrt <---+           <-----+ 1/g_sqrt <------
    %        +----------+   |           |     +----------+
    %                       |           |
    %              +--------+           <--------+
    %              |        +-----------+        |
    %              |                             |
    %              |            +---+            |
    %              +------------> K +------------+
    %                           +---+
    %
    
    if nargin == 4
        m = 1;        
        g_sqrt = inf;
    end
    [nu, ny] = size(K);
    nz = size(Dr, 1);
    nw = size(Dc, 1);
    [nyM, nuM] = size(M);
    ne = nyM - nz - ny;
    nd = nuM - nw - nu;    
    cl = lft(M, K);
    left_scale = blkdiag(sqrt(m) * Dr, eye(ne) / g_sqrt);
    right_scale = blkdiag(sqrt(m) * inv(Dc), eye(nd) / g_sqrt);
    cl = left_scale * cl * right_scale;
    if isinf(g_sqrt)
        cl = cl(1:nz, 1:nw);
    end
end