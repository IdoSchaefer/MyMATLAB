function [qk,tk,r,ierr]=mlanpro(A,n,kmax,r,qk,tk,ip1)
delta=sqrt(eps/kmax);
eta=10*eps^(3/4);
cgs=0;
elr=1;
p=rand(n,1)-0.5;
u=[];
anorm=[];
est_anorm=1;
eps1=sqrt(n)*eps/2;
np=0;nr=0;ierr=0;
%prepare for lanczos iteration
[s1,s2]=size(qk);
if s1==0
    qk=zeros(n,kmax);
    beta=zeros(kmax+1,1);alpha=zeros(kmax,1);
    q=zeros(n,1);beta(1)=norm(r);
    omega=zeros(kmax,1);omega_max=omega;omega_old=omega;
    omega(1)=0;force_reorth=0;
    j0=1;
else
    j=size(qk,2);
    qk=[qk zeros(n,kmax-j)];
    alpha=zeros(kmax+1,1);
    beta=zeros(kmax+1,1);
    alpha(1:j)=diag(tk);
    beta(2:j)=diag(tk,-1);
    q=qk(:,j);
    %reorthogonalize r
    beta(j+1)=norm(r);
    if j < kmax & beta(j+1)*delta < anorm*eps1
        fro=1;
        ierr=j;
    end
    int=1:j;
    [r, beta(j+1) ,rr]=mreorth(qk,r,beta(j+1),int,0.5,cgs);
    np=rr*j;
    nr=1;
    force_reorth=1;
    %compute gerscgorin bound on ||tk||_2 as sqrt(||tk'*tk||_1)
    if est_anorm
        anorm=sqrt(norm(tk'*tk,1));
    end
    omega=eps1*ones(kmax,1);
    omega_max=omega;
    omega_old=omega;
    j0=j+1;
end
    
%Initialize mu/nu - recurrence for monitoring loss of orthogonalization
% nu=zeros(kmax,1);mu=zeros(kmax,1);
% nu(1)=1;mu(1)=1;
% numax=zeros(kmax,1);mumax=zeros(kmax,1);
% j0=1;
%%%%%%%%%%%%%%%%%%%%%%%%%
fro=0;
%Perform Lanczos partial reorthogonalization
for j=j0:kmax
   q_old=q;
   if beta(j)==0
       q=r;
   else
       q=r/beta(j);
   end
   qk(:,j)=q;
   if ip1==0
       u=A*q;
   else
       u=feval(A,q);
   end
   r=u-beta(j)*q_old;
   alpha(j)=q'*r;
   r=r-alpha(j)*q;
   %extend local reorthogonalization
   if elr
       if j==1
           t1=0;
           for i=1:2
               t=q'*r;
               r=r-q*t;
               t1=t1+t;
           end
           alpha(j)=alpha(j)+t1;
       elseif j>1
           t1=q_old'*r;
           t2=q'*r;
           r=r-(q_old*t1+q*t2);
           if beta(j)~=0
               beta(j)=beta(j)+t1;
           end
           alpha(j)=alpha(j)+t2;
       end
   end
   beta(j+1)=sqrt(r'*r);%Quick and dirty estimate
   %Update Gersgorin estimate of ||T_k|| if required
   if est_anorm &beta(j+1)~=0
       anorm=update_gbound(anorm,alpha,beta,j);
   end
   %update omega reccurence
   if j > 1 & ~fro & beta(j+1) ~=0
       [omega ,omega_old]=update_omega(omega,omega_old,j,alpha,beta,eps1,anorm);
       omega_max(j)=max(abs(omega));
   end
   %reorthogonalize if required
   if j > 1 & (fro|force_reorth | omega_max(j)>delta)&beta(j+1)~=0
       if fro
           int=1:j;
       else
           if force_reorth==0
               force_reorth=1; % Do forced reorth to avoid spill-over from q_{j-1}
               int=compute_int(omega,j,delta,eta,0,0,0);
           else
               force_reorth=0;
           end
       end
       [r,beta(j+1),rr]=mreorth(qk,r,beta(j+1),int,0.5,cgs);
       omega(int)=eps1;
       np=np+rr*length(int(:));
       nr=nr+1;
   else
       beta(j+1)=norm(r); % compute norm accuracy
   end
   if j < kmax& beta(j+1) < n*anorm*eps
       %if beta is "small" we deflate by setting the off-diagonal of tk 
       %to 0 and attemp to restart with a basis for a new 
       % invariant subspace by replacing r with a random starting vector
       beta(j+1)=0;
       bailout=1;
       for attemp=1:3
           r=rand(n,1)-0.5;
           if ip1==0
               r=A*r;
           else
               r=feval(A,r);
           end
           nrm=sqrt(r'*r);
           int=1:j;
           [r,nrmnew,rr]=mreorth(qk,r,nrm,int,0.5,cgs);
           np=np+rr*length(int(:));
           nr=nr+1;
           if nrnew > 0
               %a vector orthogonali orthogonal to span{qk(1:j)} was found
               %Continue iteration
               bailout=0;
               break
           end
       end
       if bailout
           ierr=-j;
           break
       else
           r=r/nrmnew; % Continue with new normalized r as starting vector.
           force_reorth=1;
           omega(:)=eps1;
           if delta > 0
               fro=0;%Turn off full reorthogonalization
           end
       end
   elseif j < kmax & ~fro & beta(j+1)*delta < anorm*eps1
       fro=1;
       ierr=j;
   end
end
%Set up tridiagonal tk in sparse matrix data structure
tk=spdiags([[beta(2:j);0] alpha(1:j) beta(1:j)],-1:1,j,j);
qk=qk(:,1:j);
work=[nr np];