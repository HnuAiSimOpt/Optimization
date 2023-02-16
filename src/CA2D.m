function rb = CA2D(RB,K,oldK,U0,nfree,nCA)
Kr = RB'*oldK*RB;
invKr = inv(Kr);
deltk = K-oldK;
cab = zeros(nfree,nCA);
ROBdotInvKr = RB*invKr;
cab(:,1) = U0;
for i = 2:nCA
  vec1 = deltk*cab(:,i-1);
  vec2 = ROBdotInvKr*(RB'*vec1);
  cab(:,i) = -vec2;
end
[rb,~,~] = svd(cab,0);
end