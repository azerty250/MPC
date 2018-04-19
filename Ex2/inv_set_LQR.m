function [X_lqr Olqr] = inv_set_LQR(X,U,A,B,N,Q,R)

G = U.A;
g = U.b;

%% LQR
K = [];
H = Q;

for i = 1:N
    K = [K; -inv(R+B'*H*B)*B'*H*A];
    H = Q+K(i,:)'*R*K(i,:)+(A+B*K(i,:))'*H*(A+B*K(i,:));
end

K = K(N,:);

%% Compute Oinf

X_lqr = Polyhedron([X.A;G*K],[X.b;g]);

Omega = X_lqr;

H = X_lqr.A;
h = X_lqr.b;

A_lqr = A+B*K;

while 1
    Omega_new = Omega;   
    F = Omega_new.A;
    f = Omega_new.b;
    
    Omega = Polyhedron([F;F*A_lqr],[f;f]);
    
    if(Omega == Omega_new)
        Olqr = Omega_new;
        return
    end
    
end

end