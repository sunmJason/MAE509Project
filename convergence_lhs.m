function indices = convergence_lhs(M)

    p = 2;
    
    phi = -18.41;
    h = 0.00191;
    
    param1 = [0.5 1.5]*phi;
    param2 = [0.5 1.5]*h;
    
    D = zeros(M, 2);
    E = zeros(M, 2);
    
    D(:,1) = param1(1) + (param1(2) - param1(1)).*lhsdesign(M, 1);
    D(:,2) = param2(1) + (param2(2) - param2(1)).*lhsdesign(M, 1);
    
    E(:,1) = param1(1) + (param1(2) - param1(1)).*lhsdesign(M, 1);
    E(:,2) = param2(1) + (param2(2) - param2(1)).*lhsdesign(M, 1);
    
    F = zeros(M,p,p);
    for i = 1:p
        F(:,:,i) = E;
        F(:,i,i) = D(:,i);
    end
    
    % Run the model and compute selected model output at sampled parameter
    for  j = 1:M
        yD(j,1) = ss_sltn(D(j,:));
        yE(j,1) = ss_sltn(E(j,:));
        for i = 1:p
            yF(j,i) = ss_sltn(F(j,:,i));
        end
    end
    
    % Compute sensitivity indices
    f0  = mean(yD) ;
    VARy = mean(yD.^2) - f0^2 ;
    
    for i = 1:p
        yFi = yF(:,i);
    
	    % fist order indices	
        Si(i)  = ( 1/M*sum(yD.*yFi) - f0^2 ) / VARy ; 
        % total effects indices
        STi(i) = 1 -  ( 1/M*sum(yE.*yFi) - f0^2 ) / VARy ;
    end
    
    % Plot results
    % sensitivity indices
    indices = [Si' STi'];
end

