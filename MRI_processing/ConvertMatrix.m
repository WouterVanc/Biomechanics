function [List] = ConvertMatrix(TestMatrix)

[m,n] = size(TestMatrix);
logicmatrix = TestMatrix > 0;
n_list = nnz(logicmatrix);
List = zeros(n_list,1);

c = 0;
for i = 1:m
    for j = 1:n
        if TestMatrix(i,j) ~= 0
            c = c+1;
            List(c,1) = TestMatrix(i,j);
        end
    end
end

List = List';             