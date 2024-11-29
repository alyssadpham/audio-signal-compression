load CA_01;
%%Problem 1
%load sample data
figure(1);
stem(x);hold on;plot(x);grid;

%then I find orthonormal basis vectors
figure(2);
hold on;grid;
fplot(@(x) cos(1*pi*(x-1)/7), [1,8], '-.');
fplot(@(x) cos(2*pi*(x-1)/7), [1,8], '-hexagram');
fplot(@(x) cos(3*pi*(x-1)/7), [1,8], 'r--o');
fplot(@(x) cos(4*pi*(x-1)/7), [1,8], 'pentagram-');
fplot(@(x) cos(5*pi*(x-1)/7), [1,8], '--*');
fplot(@(x) cos(6*pi*(x-1)/7), [1,8], '--o');
fplot(@(x) cos(7*pi*(x-1)/7), [1,8], ':');
fplot(@(x) cos(8*pi*(x-1)/7), [1,8], '--');

A=[1 1 1 1 1 1 1 1]';
for k=2:8
A=[A arrayfun(@(x) cos(k*pi*(x-1)/7), [1:8]')];
end

A'*A;

%Classic GS
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n,n);
for j=1:n
v=A(:,j);
for i=1:j-1
R(i,j)=Q(:,i)'*A(:,j);
v=v-R(i,j)*Q(:,i);
end
R(j,j)=norm(v);
Q(:,j)=v/R(j,j);
end

figure(3);
for k=1:8
subplot(8,1,k);
stem(Q(:,k)); axis([0 9 -1 1]); axis off; hold on;
end

for k=1:8
subplot(8,1,k);
plot(Q(:,k));
end
n_inputU=norm(U,2); %2-norm of input matrix U
n_calculatedU=norm(Q,2); %2-norm of calculated matrix U


a=Q'*x;

[sorted_a2, indices] = sort(abs(a), 'descend');
a2 = zeros(size(a));
a2(indices(1:2)) = a(indices(1:2));

x2=Q*a2;
figure(1); stem(x2,'b*'); plot(x2,'b');

[sorted_a4, indices] = sort(abs(a), 'descend');
a4 = zeros(size(a));
a4(indices(1:4)) = a(indices(1:4));

x4=Q*a4;
figure(1); stem(x4,'kx'); plot(x4,'k');

x8=Q*a;
sqrt(sum((x-x8).^2)/sum(x.^2)); %relative error between x & x8
sqrt(sum((x-x2).^2)/sum(x.^2)); %relative error between x & x2
sqrt(sum((x-x4).^2)/sum(x.^2)); %relative error between x & x4

%Now I will load an actual audio file and compress it
S = load('CA_01_STRONG.mat');
y = S.him;
figure(4);
stem(y);hold on;plot(y);grid;
aa = dct(y);

[sorted_aa,ind] = sort(abs(aa), 'descend');
n = 1;
while norm(aa(ind(1:n)))/norm(aa) < 0.99
n=n+1;
end
new_aa=aa;
new_aa(ind(n+1:end))=0;
new_y=idct(new_aa);
figure(4); hold on; plot(new_y); legend('y','new_y');
sqrt(sum((y-new_y).^2)/sum(y.^2));