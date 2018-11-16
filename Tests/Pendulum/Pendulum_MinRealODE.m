function qDot = Pendulum_MinRealODE(t,q,tData)

q1 = q(1);
q2 = q(2);

g = 9.806;

A = ...
[ 1,     0;
  0, -25/3];
 
b = ... 
[              q2;
 (5*g*sin(q1))/2];

if(rank(A)<2)
    keyboard
end
qDot = [A\b];
end
