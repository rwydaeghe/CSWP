function gridplot(C) %matrix to faces (for p matrix)
C = [C; zeros(1,size(C,1))];
B = [transpose(C); zeros(1,size(C,1))];
C = transpose(B);
pcolor(C)
colormap(gray(2))
axis ij
axis square
