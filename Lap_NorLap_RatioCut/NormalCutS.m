function [Jfn, Cfn, kfn, dbn, Qn, accn,Vn, Dn] = NormalCutS(A,K,L,c,Partition)
%normalized cut taking
m = size(A,1);
if ~exist('K','var') && ~exist('L','var')&& ~exist('c','var')
    K = min(100, ceil(m/10));
    L = 2;
    c = 0;
end;

%d = sparse(diag(sum(abs(A))+1));
%normalized laplacian
La=laplacian(A);
%N = d^(-1)*A;
%d = sparse(diag(sum(abs(A))));
%N = d-A;
[V,D]=eigs(La,K,'SA');


dbn = zeros(K-1,3);
kfn = zeros(3,1);
Qn = zeros(K-1,1);
accn = zeros(K-1,1);
temp = [10000 0];
Cfn = cell(3,1);
Jfn = cell(3,1);


for i = 2:K
    [Jt, Ct, dbt] = AdjCl(V(:,1:i),L,c);
    dbn(i-1,:) = dbt;
    if dbn(i-1,1)<temp(1)
        temp(1) = dbn(i-1,1);
        Jfn{1} = Jt;
        Cfn{1} = Ct;
        kfn(1) = i;
    elseif temp(1) == dbn(i-1,1)
        %WARNING('Multiple Solution!)');
    end;
    
    P = idx2lgc(Jt);
    %Q(i-1) = Qfunction(A,P);
    Qn(i-1) = SignQfunction(A,P);
    if Qn(i-1)>temp(2)
        temp(2) = Qn(i-1);
        Jfn{2} = Jt;
        Cfn{2} = Ct;
        kfn(2) = i;
    elseif temp(2) == Qn(i-1)
        WARNING('Multiple Solution!)');
    end;
    
    if exist('Partition','var') && sum((size(P)-size(Partition)).^2)==0;
        accn(i-1) = PartitionAccuracy(P, Partition);
    end
end;

Vn=V(:,1:kfn(2));
Dn=D(:,1:kfn(2));

figure;plot(2:K,dbn(:,1),'-o','MarkerSize',6,'MarkerFaceColor','b');
xlabel('k');
ylabel('Davies-Buldin Index');
hold on;
plot([2,K],[1,1],'k--');
hold off;
figure;errorbar(2:K,dbn(:,1),dbn(:,2)-dbn(:,1),dbn(:,1)-dbn(:,3));
xlabel('k');
ylabel('Davies-Buldin Index');
xlim([2,K]);
hold on;
plot([2,K],[1,1],'k--');
hold off;

[qx,qy] = max(Qn(logical(dbn(:,1)<1)));
figure;plot(2:K,Qn,'-o','MarkerSize',6,'MarkerFaceColor','b');
xlabel('k');
ylabel('Modularity');
hold on;
plot(qy+1,qx,'ks','MarkerSize',12)
hold off;

d = diag(D);
[dx,dy]=sort(abs(d),'descend');
figure;plot(1:K,d(dy),'-o','MarkerSize',6,'MarkerFaceColor','b');
xlabel('k');
ylabel('Eigenvalues');
hold on;
plot([1,K],[0,0],'k--');
hold off;

