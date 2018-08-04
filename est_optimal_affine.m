function O =est_optimal_affine(p2,p1)
A = [p2(1,1),p2(1,2),1,0,0,0;
        0,0,0,p2(1,1),p2(1,2),1;
        p2(2,1),p2(2,2),1,0,0,0;
        0,0,0,p2(2,1),p2(2,2),1;
        p2(3,1),p2(3,2),1,0,0,0;
        0,0,0,p2(3,1),p2(3,2),1];
b = [p1(1,1);p1(1,2);p1(2,1);p1(2,2);p1(3,1);p1(3,2)];
B = A\b;
O = [B(1),B(2),B(3);
        B(4),B(5),B(6)];
end