function [Floe] = Weld_Floes(Floe,Nb,dhdt,Amax)
%%This function takes in two floe field and based upon a specified
%%probability function welds interacting floes together

id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

i = Nb+1; count = 0; FloeOld = Floe;
floenew = [];

for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
    FloeOld(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

%Set probability function
hvsd = @(x) [0.5*(x == 0) + (x > 0)];
ramp = @(frac) hvsd(frac)*frac;
SimpMin = @(A) 3*log10(A);

%Loop through all floes to determine which will weld
while i < length(Floe)
    
    if Floe(i).alive && ~isempty(Floe(i).potentialInteractions)
        c = unique(cat(1,Floe(i).potentialInteractions.floeNum));
        c = c(c<length(Floe));
        c = c(c>Nb);
        c = c(c>i);
        if ~isempty(c)
            OverlapA = area(intersect(Floe(i).poly,[Floe(c).poly]));
            [Areamax,k] = max(OverlapA);
            j = c(k);
            if Floe(j).alive && Areamax > 0
                
                frac = (Floe(i).area+Floe(j).area)/Amax;
                p = rand(1);
                
                %If probability is met then weld them together
                if p <dhdt*OverlapA/Floe(i).area %ramp((1-frac)^5*dhdt)
                    floe = FuseFloes(Floe(i),Floe(j));
                    count = count+1;
%                     if (Floe(i).mass+Floe(j).mass)*floe.area/(floe.mass*(Floe(i).area+Floe(j).area)) - 1 > 0.001
%                         xx = 1;
%                         xx(1) = [1 2];
%                     end
                    
                    if abs(floe.area/area(polyshape(floe.c_alpha'))-1)>1e-3
                        xx = 1;
                        xx(1) =[1 2];
                    end
                    
                    Floe(j).alive  = 0;
                    if length(floe.poly.Vertices) > SimpMin(floe.area)
                        floe2 = FloeSimplify(floe,false);
                    else
                        floe2 = floe;
                    end
                    for jj = 1:length(floe2)
                        if jj == 1
                            Floe(i) = floe2(jj);
                        else
                            floenew = [floenew floe2(jj)];
                        end
                        if abs(floe2(jj).area/area(polyshape(floe2(jj).c_alpha'))-1)>1e-3
                            xx = 1;
                            xx(1) =[1 2];
                        end
                    end
                    
                end
                
            end
        end
    end
    i = i+1;
end

for ii =1:length(Floe)
    Floe(ii).poly = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi Floe(ii).Yi]);
end

Floe = [Floe floenew];
live = cat(1,Floe.alive);
Floe(live == 0) = [];

Floe=rmfield(Floe,{'poly'});

warning('on',id)
warning('on',id3)

for ii = 1:length(Floe)
    if abs(Floe(ii).area/area(polyshape(Floe(ii).c_alpha'))-1)>1e-3
        xx = 1;
        xx(1) =[1 2];
    end
end

end

