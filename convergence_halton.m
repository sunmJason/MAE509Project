function indices = convergence_halton(M)
    phi = -18.41;
    h = 0.00191;
    
    param1 = [0.5 1.5]*phi;
    param2 = [0.5 1.5]*h;
    
    p = 2;
    A = zeros(M, 2);
    B = zeros(M, 2);
    halton_set = net(haltonset(4), M);
    
    A(:,1) = param1(1) + (param1(2) - param1(1)).*halton_set(:, 1);
    A(:,2) = param2(1) + (param2(2) - param2(1)).*halton_set(:, 2);
    
    B(:,1) = param1(1) + (param1(2) - param1(1)).*halton_set(:, 3);
    B(:,2) = param2(1) + (param2(2) - param2(1)).*halton_set(:, 4);
    
    C = zeros(M,p,p);
    for i = 1:p
        C(:,:,i) = B;
        C(:,i,i) = A(:,i);
    end

    % Run the model and compute selected model output at sampled parameter
    for  j = 1:M
        yA(j,1) = ss_sltn(A(j,:));
        yB(j,1) = ss_sltn(B(j,:));
        for i = 1:p
            yC(j,i) = ss_sltn(C(j,:,i));
        end
    end

    % Compute sensitivity indices
    f0  = mean(yA) ;
    VARy = mean(yA.^2) - f0^2 ;
    
    for i = 1:p
        yCi = yC(:,i);
    
	    % fist order indices	
        Si(i)  = ( 1/M*sum(yA.*yCi) - f0^2 ) / VARy ; 
        % total effects indices
        STi(i) = 1 -  ( 1/M*sum(yB.*yCi) - f0^2 ) / VARy ;
    end
    indices = [Si' STi'];
end


