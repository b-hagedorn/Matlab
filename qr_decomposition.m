function [ QR, a ] = qrDecomposition( A )
    [m n] = size(A);
    submatrix = A;
    QR = zeros(m, n);
    for i=1:n 
        v = calculateReflectionPlane(submatrix);
        QR(i:m, i) = v;
        submatrix = createSubmatrix(submatrix);
        calculateSubQ(v, i-1);
    end

end

function v = calculateReflectionPlane( A )
    sgn = sign(A(1,1));
    euclidNorm = norm(A(:,1),2);
    v = A(:,1);
    v(1) = v(1) + sgn * euclidNorm;
end

function A = createSubmatrix( A )
    [m n] = size(A);
    A = A(2:m, 2:n)
end

function subQ = calculateSubQ( v, additionalUnitMatrixDimension )
    unitMatrixDimension = size(v);
    unitMatrix = eye(unitMatrixDimension);
    tempQ = unitMatrix - (2/(v'*v))*(v*v');
    
    subQDimension = unitMatrixDimension + additionalUnitMatrixDimension
    subQ = eye(subQDimension);
    subQ(additionalUnitMatrixDimension:subQDimension, additionalUnitMatrixDimension:subQDimension) = tempQ;
end

