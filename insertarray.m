function x12 = insertarray(dim,x1,x2,i)
% Insert data at specific dimension and locations in a multidimensional array
% <Syntax>
% x12 = insertarray(dim,x1,x2,i)
% 
% <Input>
% dim - Dimension to operate along, specified as a positive integer scalar.
% x1 - Original array, specified as a vector, matrix, or multidimensional 
%     array.
% x2 - Array to insert, x2 and x1 should have the same data type and 
%     compatible sizes (the lengths of the dimensions match except for the
%     operating dimension dim).
% i - Locations in x1 to insert x2, specified as a scalar or a vector. If 
%     i is a scaler, x2 is inserted after position floor(i) in x1. For
%     example, specify i = 3.5 to insert x2 between page 3 and 4 in x1. You
%     can also specify inserted locations for each page of x2 by providing
%     i as a vector. This vector must be the same size as the operating
%     dimension of x2.
% 
% <Output>
% x12 - New data after inserting x2 into x1.
% 
% <Examples>
% example 1: inserting two rows of x2 after row 2 and row 3 of x1
% x1 = char(randi(9,4,3)+65);
% x2 = char(randi(9,2,3)+65);
% x12 = insertarray(1,x1,x2,[2 3]);
% example 2: inserting two pages of x2 before page 1 and page 3 of x1
% x1 = randi(9,2,3,3);
% x2 = randi(9,2,3,2);
% x12 = insertarray(3,x1,x2,[0.5 2.5]);
% 
% <Copyright>
% Author:   Yang Liu
% Contact:  liuyang-y2003@foxmail.com
% Update:   2022-06-01
% Version:  1.0.0

if isscalar(i)
    i = ones(1,size(x2,dim))*i;
else if iscolumn(i)
        i = i';
    end
end

idx = [(1:size(x1,dim)),i];
idx = [idx;[zeros(1,size(x1,dim)),ones(1,size(x2,dim))]];
idx = sortrows(idx')';
x12 = cat(dim,x1,x2);
a = cell(1,ndims(x12));

for j = 1:ndims(x12) 
    a(j) = {1:size(x12,j)};
end

a(dim) = {idx(2,:)==0}; 
x12(a{:}) = x1;
a(dim) = {idx(2,:)==1}; 
x12(a{:}) = x2;
