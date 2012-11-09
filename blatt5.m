function blatt5()
A = [3 -6; 4 -8; 0 1];
[QR, alpha] = myQR( A );
x = solveWithQR( QR, alpha, b )
end

function [ QR, alpha ] = myQR( A )
    [m n QR R alpha] = initializeVariables(A);
        
    for i=1:n
        % calculate and register v
        columnAndRowErased = createSubmatrix(R, i);
        v = calculateReflectionPlane( columnAndRowErased );
        QR(i:m, i) = v;
        
        R = calculateR( R, v, i);
        alpha(i) = R(i,i);       
    end   
    
    QR = registerR( QR, R );
end

function [m n QR R alpha] = initializeVariables( A )
    [m n] = size(A);
    QR = zeros(m, n);
    R = A;
    alpha = zeros( min(m, n), 1);
    
end

function v = calculateReflectionPlane( M )
    sigma = calculateSigma( M );
    euclidNorm = norm(M(:,1),2);
    v = M(:,1);
    v(1) = v(1) + sigma * euclidNorm;
end

function sigma = calculateSigma( M )
    sigma = sign(M(1,1));
    if (sigma == 0)
        sigma = 1;
    end
end

function A = createSubmatrix( A, i )
    if (i ~= 0)
        [m n] = size(A);
        A = A((i):m, (i):n);
    end
end

function R = calculateR( R, v, insertPosition )
    [m n] = size(R);

    rightSizeV = zeros(m, 1);
    rightSizeV( insertPosition:m, 1) = v;
    R = R - (2/(rightSizeV'*rightSizeV)) * (rightSizeV*(rightSizeV'*R));
end

function QR = registerR( QR, R )
    [m n] = size(R);
    
    for row = 1 : (min(m, n)-1)
        for column = (row + 1) : min(m, n)
            QR(row, column) = R(row, column);
        end
    end
end



