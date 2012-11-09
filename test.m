function test()
A = [3 -6; 4 -8; 0 1];
[QR, alpha] = qrDecomposition( A )
end

function [ QR, a ] = qrDecomposition( A )
    [m n] = size(A);
    QR = zeros(m, n);
    
    [QR, Q] = fillQRWithReflectionPlanesAndCalculateQ(QR, A);
    
    R=Q*A;
    QR = fillQRWithR(QR, R)
    
    a = 1;
end

function [QR, Q] = fillQRWithReflectionPlanesAndCalculateQ( QR, A )
    [m n] = size(A);

    % erster schritt auﬂerhalb der schleife
    v = calculateReflectionPlane(A);
    QR(1:m, 1) = v;
    Q = calculateSubQ(v, 0);
    
    for i=2:n 
        QA= Q*A
        v = calculateReflectionPlane(createSubmatrix(QA, i));
        QR(i:m, i) = v;
       
        % Q_i werden aufmultipliziert
        Q = calculateSubQ(v, i-1) * Q;
    end
end

function QR = fillQRWithR(QR, R)
    
    
end

function v = calculateReflectionPlane( M )
    sgn = sign(M(1,1));
    if (sgn == 0)
        sgn = 1;
    end
    euclidNorm = norm(M(:,1),2);
    v = M(:,1);
    v(1) = v(1) + sgn * euclidNorm;
end

function A = createSubmatrix( A, i )
    [m n] = size(A);
    A = A((i):m, (i):n);
end

function subQ = calculateSubQ( v, additionalUnitMatrixDimension )
    [unitMatrixDimension n] = size(v);
    unitMatrix = eye(unitMatrixDimension);
    
    tempQ = unitMatrix - (2/(v'*v))*(v*v');
    
    subQDimension = unitMatrixDimension + additionalUnitMatrixDimension;
    subQ = eye(subQDimension);
    
    subQ((additionalUnitMatrixDimension+1):subQDimension, (additionalUnitMatrixDimension+1):subQDimension) = tempQ;
end

