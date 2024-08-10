fsnew = 100; 
tnew = 0:(1/fsnew):82;

% target
r = [20 33.33333 100 20];
ti = [2 4.5 5.5 8;
      3 4.5 5.5 7;
      4 4.5 5.5 6;
      2 4.5 5.5 8];

for j = 1:4
    
    Target = zeros(sum(tnew < 8),1);
    Target(tnew < ti(j,1)) = 0;
    Target(tnew > ti(j,1) & tnew < ti(j,2)) = r(j)*tnew(tnew > ti(j,1) & tnew < ti(j,2)) - ti(j,1)*r(j);
    Target(tnew >= ti(j,2) & tnew <= ti(j,3)) = 50;
    Target(tnew > ti(j,3) & tnew < ti(j,4)) = -r(j)*tnew(tnew > ti(j,3) & tnew < ti(j,4)) + ti(j,4)*r(j);

    if j == 4
        Target(tnew > ti(j,3) & tnew < 8) = 0;
    end
    
    rampTarget(j,:) = [repmat(Target,10,1); zeros(201,1)]';
end

% plot(tnew, rampTarget);

cd([codefolder, '\data\torque'])
save('RampTarget.mat')