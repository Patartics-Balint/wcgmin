function cl = close_wcg_min_loop(M, K, Dr, Dc, g_sqrt)
    % CLOSE_WCG_MIN_LOOP Forms the closed loop interconnection used in the
    % worst case minimisation.
    %
    % Syntax
    %   cl = CLOSE_WCG_MIN_LOOP(M, K, Dr, Dc)
    %   cl = CLOSE_WCG_MIN_LOOP(M, K, Dr, Dc, g_sqrt)
    %
    % Description
    %
    %   cl = CLOSE_WCG_MIN_LOOP(M, K, Dr, Dc, g_sqrt) forms the closed    
    %   loop interconnection below. The default value of g_sqrt is 1.
    %
    %            +----+      +-----------+    +---------+
    %     <------+ Dr <------+           <----+ inv(Dc) <------
    %            +----+      |           |    +---------+
    %                        |           |
    %         +----------+   |     M     |    +----------+
    %     <---+ 1/g_sqrt <---+           <----+ 1/g_sqrt <-----
    %         +----------+   |           |    +----------+
    %                        |           |
    %               +--------+           <--------+
    %               |        +-----------+        |
    %               |                             |
    %               |            +---+            |
    %               +------------> K +------------+
    %                            +---+
        
    if nargin == 4
        g_sqrt = 1;
    end
    nz = size(Dr, 1);
    nw = size(Dc, 1);    
    cl = lft(M, K);
    [ny_cl, nu_cl] = size(cl);
    ne = ny_cl - nz;
    nd = nu_cl - nw;
    left_scale = blkdiag(Dr, eye(ne) / g_sqrt);
    right_scale = blkdiag(Dc, eye(nd) * g_sqrt);
    cl = left_scale * cl / right_scale;
end
